import nanome
from nanome.util import Logs, Octree, Process
from nanome.api.structure import Complex
from nanome.util.stream import StreamCreationError
from nanome.util.enums import NotificationTypes, StreamType

from nanome._internal._structure._io._pdb.save import Options as PDBOptions
from nanome._internal._structure._bond import _Bond

from openff.toolkit.topology import Molecule
from openff.toolkit.topology import Topology
from openmmforcefields.generators import GAFFTemplateGenerator
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openmmforcefields.generators import SystemGenerator

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk.openmm.app.element import *
from pdbfixer.pdbfixer import PDBFixer, proteinResidues, dnaResidues, rnaResidues, _guessFileFormat

from rdkit.Chem import AllChem
from rdkit import Chem

import numpy as np
import tempfile
from collections import deque
from functools import partial
import os, sys, threading, io

PACKET_QUEUE_LEN = 20

SDFOPTIONS = Complex.io.SDFSaveOptions()

pdb_options = PDBOptions()
pdb_options.write_bonds = True

pdb_options_no_het = PDBOptions()
pdb_options_no_het.write_bonds = True
pdb_options_no_het.write_het_bonds = False
pdb_options_no_het.write_het_atoms = False

# hack to fix convert_to_frames killing atom indices:
_atom_shallow_copy = nanome._internal._structure._Atom._shallow_copy
def _atom_shallow_copy_fix(self, *args):
    atom = _atom_shallow_copy(self, *args)
    atom._index = self._index
    return atom
nanome._internal._structure._Atom._shallow_copy = _atom_shallow_copy_fix


class MinimizationProcess():
    def __init__(self, plugin):
        self.__plugin = plugin
        self._is_running = False
        self.__process_running = False
        self.__stream = None
        #TODO: change this to AdvancedSettings value
        #FF for openbabel: GeneralAmber (Gaff), Ghemical, MMFF94, MMFF94s, Universal
        #FF for OpenMM: GeneralAmber (Gaff), Amber (ff14SB), CHARMM, OpenForceField (openff & smirnoff)
        # https://github.com/openmm/openmmforcefields
        self.__forcefield_names = ["amber/protein.ff14SB.xml", "amber/tip3p_standard.xml", "amber/tip3p_HFE_multivalent.xml",\
        "charmm/toppar_all36_prot_model.xml", "gaff-2.11", "openff-1.3.0", "smirnoff99Frosst-1.1.0"]
        self._updateRate = 100

    def _get_bond_type(self, kind):
        if kind == _Bond.Kind.CovalentSingle:
            return Single
        if kind == _Bond.Kind.CovalentDouble:
            return Double
        if kind == _Bond.Kind.CovalentTriple:
            return Triple
        return None

    def delete_alternate_atoms(self, topology, positions):
        modeller = Modeller(topology, positions)
        delete_atoms = []
        for chain in topology.chains():
            for indexInChain, residue in enumerate(chain.residues()):
                atom_names = []
                for atom in residue.atoms():
                    if atom.name in atom_names:
                        delete_atoms.append(atom)
                    else:
                        atom_names.append(atom.name)

        modeller.delete(delete_atoms)
        return (modeller.getTopology(), modeller.getPositions())

    def pdbfixer_complexes(self, complex_list):
        fixed_complexes = []

        for complex in complex_list:
            for residue in complex.residues:
                atoms = residue._atoms
                for i in range(len(atoms) - 1, -1, -1):
                    if atoms[i].molecular.is_het == True:
                        del atoms[i]

        for complex in complex_list:

            temp_pdb = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
            complex.io.to_pdb(temp_pdb.name, pdb_options)

            fixer = PDBFixer(filename=temp_pdb.name)
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.removeHeterogens(False)
            fixer.addMissingHydrogens(7.0)

            (topology, positions) = self.delete_alternate_atoms(fixer.topology, fixer.positions)
            with open(temp_pdb.name, 'w') as pdb_file:
                PDBFile.writeFile(topology, positions, pdb_file)

            fixed_complex = nanome.structure.Complex.io.from_pdb(path=temp_pdb.name)
            fixed_complex.index = complex.index
            fixed_complex.position = complex.position
            fixed_complex.rotation = complex.rotation
            fixed_complex.molecular.name = complex.molecular.name
            fixed_complex.rendering.visible = True
            fixed_complexes.append(fixed_complex)

        return fixed_complexes

    def start_process(self, workspace, ff, steps, steepest):
        def on_stream_creation(stream, error):
            if error != StreamCreationError.NoError:
                Logs.error("Error while creating stream")
                return

            self.__stream = stream

            #Init openmm system

            self._custom_generator = None
            selected_forcefields = []
            forcefield_kwargs = { 'constraints' : HBonds, 'rigidWater' : True, 'removeCMMotion' : False, 'hydrogenMass' : 4*amu }

            for f in self.__forcefield_names:
                if f.startswith(ff):
                    selected_forcefields.append(f)

            if len(selected_forcefields) == 1:
                curff = selected_forcefields[0]
                if curff.startswith("gaff") or curff.startswith("openff") or curff.startswith("smirnoff"):
                    self._custom_generator = SystemGenerator(forcefields=["amber/protein.ff14SB.xml", "amber/tip3p_standard.xml"], small_molecule_forcefield=curff, forcefield_kwargs=forcefield_kwargs, cache='db.json')

            Logs.debug("Selected forcefields:", selected_forcefields)

            self.__integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

            if self._custom_generator == None: #Not a general forcefield 
                ommpdb = PDBFile(input_file.name)
                self.__forcefield = ForceField(*selected_forcefields)
                self.__openmm_system = self.__forcefield.createSystem(ommpdb.topology)
                (topology, positions) = self.delete_alternate_atoms(ommpdb.topology, ommpdb.positions)

                ommpdb.topology = topology
                ommpdb.positions = positions

                self.__openmm_topology = ommpdb.topology
                self.__openmm_positions = ommpdb.positions

            else: # General forcefields (Gaff/openff/smirnoff)
                # mols = Chem.SDMolSupplier(self.__het_sdf_path)
                # m = mols[0] #TODO: check if we have more than one molecule
                # mh = Chem.AddHs(m)
                # AllChem.EmbedMolecule(mh)  #generate a 3D conformation to avoid the errors before
                # Chem.rdmolops.AssignAtomChiralTagsFromStructure(mh)    #very important, else you'll repeat the errors from before

                # omm_mol = Molecule.from_rdkit(mh)
                # off_topology = omm_mol.to_topology()
                # self.__openmm_topology = off_topology.to_openmm()
                # self.__openmm_positions = omm_mol.conformers[0]._value
                # self.__openmm_system = self._custom_generator.create_system(self.__openmm_topology, molecules=omm_mol)
                from openff.toolkit.typing.engines.smirnoff import ForceField
                force_field = ForceField('openff_unconstrained-1.0.0.offxml')
                ligand_off_molecule = Molecule(self.__het_sdf_path)
                self.__openmm_topology = ligand_off_molecule.to_topology()
                self.__openmm_system = force_field.create_openmm_system(self.__openmm_topology)
                self.__openmm_positions = ligand_off_molecule.conformers[0]


            #By default, OpenMM picks the fastest platform
            self.__simulation = Simulation(self.__openmm_topology, self.__openmm_system, self.__integrator)

            #Assign zero mass to all __atom_to_restrain when creating the openmm system
            if self.__atom_to_restrain and len(self.__atom_to_restrain) > 0:
                #setup contraints to fix atoms in position
                self.__restraint = HarmonicBondForce()#used to fix an atom in one position
                self.__openmm_system.addForce(self.__restraint)
                nonbonded = [f for f in self.__openmm_system.getForces() if isinstance(f, NonbondedForce)][0]
                for i in self.__atom_to_restrain:
                    j = self.__openmm_system.addParticle(0)
                    nonbonded.addException(i, j, 0, 1, 0)
                    self.__restraint.addBond(i, j, 0*nanometers, 100*kilojoules_per_mole/nanometer**2)
                    self.__openmm_positions.append(self.__openmm_positions[i])
                Logs.debug("Restraining", len(self.__atom_to_restrain), "atoms")

            self.__simulation.context.setPositions(self.__openmm_positions)

            self.__process_running = True
            self._is_running = True

        input_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
        self.__output_lines = []
        self.__updates_done = {}
        self.__packet_id = 0
        self._steps = steps
        self._cur_step = 0
        self._lastEnergy = None

        self.__atom_to_restrain = []

        # fixed = self.pdbfixer_complexes(workspace.complexes)
        (saved_atoms, indices) = self.__save_atoms(input_file.name, workspace, workspace.complexes)
        if saved_atoms == None:
            self.__on_process_error("No atom selected")
            self.__plugin.send_notification(nanome.util.enums.NotificationTypes.warning, "No atom selected")
            self.stop_process()
            return

        self.__stream = self.__plugin.create_writing_stream(indices, StreamType.position, on_stream_creation)

    def stop_process(self):
        # if self.__process_running:
        #     self.__process.stop()
        if self.__stream is not None:
            self.__stream.destroy()
            self.__stream = None
        self._is_running = False
        self.__plugin.minimization_done()

    def update(self):
        if not self._is_running or \
                self.__packet_id > PACKET_QUEUE_LEN and \
                not self.__updates_done[self.__packet_id - PACKET_QUEUE_LEN]:
            return

        if not self.__process_running:
            Logs.debug('Minimization complete')
            self.stop_process()
        
        if self._is_running and self.__process_running:
            self.__simulation.minimizeEnergy(maxIterations=self._updateRate)
            energy = self.__simulation.context.getState(getEnergy=True).getPotentialEnergy()
            positions = self.__simulation.context.getState(getPositions=True).getPositions()
            self.minimization_result(positions, energy)

            self._cur_step += self._steps / self._updateRate 

            if self._lastEnergy is not None and math.isclose(self._lastEnergy, energy._value, abs_tol=1e-5):
                Logs.debug('Minimization converged')
                self.stop_process()

            if self._cur_step >= self._steps:
                Logs.debug('Minimization complete, maximum step')
                self.stop_process()

            self._lastEnergy = energy._value


    def __on_process_error(self, error):
        Logs.warning('Error in minimization process:')
        Logs.warning(error)


    # def __match_and_move(self, complex):
    #     positions = [0] * (len(self.__atom_position_index_by_serial) * 3)
    #     complex_absolute_to_relative = complex.get_workspace_to_complex_matrix()
    #     for atom in complex.atoms:
    #         if atom.serial in self.__atom_position_index_by_serial:
    #             (idx, complex) = self.__atom_position_index_by_serial[atom.serial]
    #             atom_relative_pos = complex_absolute_to_relative * atom.position
    #             positions[idx * 3] = atom_relative_pos.x
    #             positions[idx * 3 + 1] = atom_relative_pos.y
    #             positions[idx * 3 + 2] = atom_relative_pos.z
    #     if self.__stream == None:
    #         return
    #     self.__updates_done[self.__packet_id] = False
    #     self.__stream.update(positions, partial(self.__update_done, self.__packet_id))
    #     self.__packet_id += 1


    def __get_atom_symbol(self, name, atoms_nb):
        upper = name.upper()
        if upper.startswith('CL'):
            return chlorine
        elif upper.startswith('NA'):
            return sodium
        elif upper.startswith('MG'):
            return magnesium
        elif upper.startswith('BE'):
            return beryllium
        elif upper.startswith('LI'):
            return lithium
        elif upper.startswith('K'):
            return potassium
        elif upper.startswith('ZN'):
            return zinc
        elif (atoms_nb == 1 and upper.startswith('CA')):
            return calcium
        else:
            return Element.getBySymbol(upper[0])

    def __save_atoms(self, path, workspace, complexes):
        selected_atoms = Octree()
        count_selected_atoms = 0
        selected_atoms_per_complex = {}
        close_atoms_per_complex = {}
        fcomplexes = []

        #Add all atoms to an Octree
        for compl in complexes:
            fcomplex = compl.convert_to_frames()
            fcomplexes.append(fcomplex)
            complex_local_to_workspace_matrix = fcomplex.get_complex_to_workspace_matrix()
            selected_atoms_per_complex[fcomplex.name] = []
            for atom in fcomplex.atoms:
                if atom.selected:
                    atom_absolute_pos = complex_local_to_workspace_matrix * atom.position
                    selected_atoms.add(atom, atom_absolute_pos)
                    selected_atoms_per_complex[fcomplex.name].append(atom) 
                    count_selected_atoms += 1

        Logs.debug("Selected atoms:",count_selected_atoms)
        if count_selected_atoms < 2:
            return (None, None)

        # self.__atom_position_index_by_serial = dict()
        atom_by_index = dict()
        saved_atoms = []
        saved_atoms_indices = []
        atom_position_index = 0
        serial = 1
        found_atoms = []

        # Check for close atoms (7A radius) that are not in selection
        # These atoms will later be fixed in position when creating the OpenMM system
        for curComplex in fcomplexes:
            complex_local_to_workspace_matrix = curComplex.get_complex_to_workspace_matrix()
            molecule = curComplex._molecules[curComplex.current_frame]
            for atom in molecule.atoms:
                #check only not selected atoms
                if curComplex.name in selected_atoms_per_complex and atom not in selected_atoms_per_complex[curComplex.name]:
                    atom_absolute_pos = complex_local_to_workspace_matrix * atom.position
                    selected_atoms.get_near_append(atom_absolute_pos, 7, found_atoms, max_result_nb=1)
                    if len(found_atoms) > 0 and atom not in selected_atoms_per_complex[atom.complex]:
                        if not atom.complex.name in close_atoms_per_complex:
                            close_atoms_per_complex[atom.complex.name] = []
                        if not atom in close_atoms_per_complex[atom.complex.name]:
                            close_atoms_per_complex[atom.complex.name].append(atom)
                        if not atom in __atom_to_restrain:
                            self.__atom_to_restrain.add(atom)
                    found_atoms.clear()


        Logs.debug("Close to selection atoms:", len(self.__atom_to_restrain))

        #Create an empty complex to store hetero atoms
        het_complex = nanome.structure.Complex()
        het_molecule = nanome.structure.Molecule()
        het_chain = nanome.structure.Chain()
        het_residue = nanome.structure.Residue()
        het_complex.add_molecule(het_molecule)
        het_molecule.add_chain(het_chain)
        het_chain.add_residue(het_residue)
        het_bonds = []

        openmm_main_pdb = None
        PDBFile._loadNameReplacementTables()

        count_het = 0

        #Convert Nanome complex to OpenMM topologies, one per molecule
        #Only keep selected or close to selection atoms and not hetero atoms
        #Hetero atoms are stored in a Nanome complex
        for curComplex in fcomplexes:
            for molecule in curComplex.molecules:
                cur_topo = Topology()
                added_atoms = dict()
                containChains = False
                positions = []
                for chain in molecule.chains:
                    sim_chain = cur_topo.addChain()
                    containChains = True
                    for residue in chain.residues:
                        residueName = residue.name
                        if residueName in PDBFile._atomNameReplacements:
                            atomReplacements = PDBFile._atomNameReplacements[residueName]
                        else:
                            atomReplacements = {}
                        sim_residue = cur_topo.addResidue(residue.name, sim_chain)
                        for atom in residue.atoms:
                            #Atom is selected or close to the selection and not an het/ligand
                            if not atom.is_het and (atom in selected_atoms_per_complex[curComplex.name] or atom in close_atoms_per_complex[curComplex.name]):
                                symbol = self.__get_atom_symbol(atom.name, len(residue._atoms))
                                atom_name = atom.name
                                if atom_name in atomReplacements:
                                    atom_name = atomReplacements[atom_name]
                                sim_atom = cur_topo.addAtom(atom_name, symbol, sim_residue)
                                added_atoms[atom.index] = sim_atom
                                position = atom.position
                                positions.append(position.x * 0.1 * nanometer)
                                positions.append(position.y * 0.1 * nanometer)
                                positions.append(position.z * 0.1 * nanometer)
                                #positions.append(Vec3(position.x * 0.1 * nanometer, position.y * 0.1 * nanometer, position.z * 0.1 * nanometer))
                            if atom.is_het:
                                het_residue.add_atom(atom)
                                count_het+=1
                                for bond in atom.bonds:
                                    if bond not in het_bonds:
                                        het_bonds.append(bond)
                                        het_residue.add_bond(bond)
                if containChains:
                    for chain in molecule.chains:
                        for residue in chain.residues:
                            for bond in residue.bonds:
                                if bond.atom1.index in added_atoms and bond.atom2.index in added_atoms:
                                    atom1 = added_atoms[bond.atom1.index]
                                    atom2 = added_atoms[bond.atom2.index]
                                    btype = self._get_bond_type(bond.kind)
                                    cur_topo.addBond(atom1, atom2, btype)
                    n = int(len(positions)/3)
                    positions = np.reshape(positions, (n, 3))
                    if openmm_main_pdb == None:
                        openmm_main_pdb = Modeller(cur_topo, positions)
                    else:
                        openmm_main_pdb.add(cur_topo, positions)
                    Logs.debug(openmm_main_pdb.getPositions())
                
        Logs.debug("Hetero atoms count =", count_het)

        if count_het > 0:
            temp_sdf = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf')
            het_complex.io.to_sdf(temp_sdf.name, SDFOPTIONS)
            self.__het_sdf_path = temp_sdf.name
            Logs.debug("Wrote",count_het,"HET atoms to sdf file", temp_sdf.name)

        if openmm_main_pdb == None:
            Logs.error("Something went wrong processing minimization files")
            return (None, None)
        PDBFile.writeFile(openmm_main_pdb.topology, openmm_main_pdb.positions, open(path, "w"))
        Logs.debug("Wrote input file:", path)
        return (saved_atoms, saved_atoms_indices)

    def minimization_result(self, positions, energy):
        Logs.debug("Minimization step", energy)
        new_positions = []
        for position in positions:
            coords = [c._value * 10 for c in position]
            new_positions.extend(coords)

        if self.__stream == None:
            return
        self.__updates_done[self.__packet_id] = False
        self.__stream.update(new_positions, partial(self.on_result_processed, self.__packet_id))
        self.__packet_id += 1

    def on_result_processed(self, packetid):
        self.__updates_done[packetid] = True

        # if self._is_running:
            # self.__simulation.step(AdvancedSettings.instance.simulation_reporter_interval)
