import nanome
from nanome.util import Logs, Octree, Process
from nanome.api.structure import Complex
from nanome.util.stream import StreamCreationError
from nanome.util.enums import NotificationTypes, StreamType

from nanome._internal._structure._io._pdb.save import Options as PDBOptions

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk.openmm.app.element import *
from pdbfixer.pdbfixer import PDBFixer, proteinResidues, dnaResidues, rnaResidues, _guessFileFormat

import tempfile
from collections import deque
from functools import partial
import os, sys, threading, io

PACKET_QUEUE_LEN = 20

IS_WIN = sys.platform.startswith('win')
SDFOPTIONS = Complex.io.SDFSaveOptions()
SDFOPTIONS.write_bonds = True

pdb_options = PDBOptions()
pdb_options.write_bonds = True

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
        # self.__platform = Platform.getPlatformByName("CUDA")
        # self.__platform_properties = {'CudaPrecision': 'single'}
        self.__forcefield_names = ["amber99sb.xml"]
        self.__forcefield = ForceField(*self.__forcefield_names)
        #OpenMM minimizer needs one but does not use it
        self._updateRate = 100

    @staticmethod
    def get_bond_type(kind):
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
            ommpdb = PDBFile(input_file.name)

            self.__openmm_system = self.__forcefield.createSystem(ommpdb.topology)
            self.__integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

            (topology, positions) = self.delete_alternate_atoms(ommpdb.topology, ommpdb.positions)

            ommpdb.topology = topology
            ommpdb.positions = positions

            #self.__simulation = Simulation(ommpdb.topology, self.__openmm_system, self.__integrator, self.__platform, self.__platform_properties)
            #By default, OpenMM picks the fastest platform
            self.__simulation = Simulation(ommpdb.topology, self.__openmm_system, self.__integrator)

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
                    ommpdb.positions.append(ommpdb.positions[i])
                Logs.debug("Restraining", len(self.__atom_to_restrain), "atoms")

            self.__simulation.context.setPositions(ommpdb.positions)

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
        Logs.debug("Wrote input file:", input_file.name)

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
        # These atoms will be fixed in position when creating the OpenMM system
        for curComplex in fcomplexes:
            complex_local_to_workspace_matrix = curComplex.get_complex_to_workspace_matrix()
            molecule = curComplex._molecules[curComplex.current_frame]
            for atom in molecule.atoms:
                #check only not selected atoms
                if curComplex.name in selected_atoms_per_complex and atom not in selected_atoms_per_complex[curComplex.name]:
                    atom_absolute_pos = complex_local_to_workspace_matrix * atom.position
                    selected_atoms.get_near_append(atom_absolute_pos, 7, found_atoms, max_result_nb=1)
                    if len(found_atoms) > 0 and atom not in selected_atoms_per_complex[atom.complex]:
                        if not atom.complex in close_atoms_per_complex:
                            close_atoms_per_complex[atom.complex] = []
                        if not atom in close_atoms_per_complex[atom.complex]:
                            close_atoms_per_complex[atom.complex].append(atom)
                        if not atom in __atom_to_restrain:
                            self.__atom_to_restrain.add(atom)
                    found_atoms.clear()


        #Convert selected atoms to an openmm system with separated complexes
        openmm_main_pdb = None
        for curComplex in fcomplexes:
            atoms = selected_atoms_per_complex[curComplex.name]
            if len(atoms) > 0:
                # This cannot work for now, bug is getting fixed
                # pdb_options.only_save_these_atoms = atoms
                temp_pdb = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')

                curComplex.io.to_pdb(temp_pdb.name, pdb_options)
                temp_openmm_pdb = PDBFile(temp_pdb.name)
                if openmm_main_pdb == None:
                    openmm_main_pdb = Modeller(temp_openmm_pdb.topology, temp_openmm_pdb.positions)
                else:
                    openmm_main_pdb.add(temp_openmm_pdb.topology, temp_openmm_pdb.positions)
                for a in atoms:
                    saved_atoms_indices.append(a.index)

        if openmm_main_pdb == None:
            Logs.error("Something went wrong processing minimization files")
            return (None, None)
        PDBFile.writeFile(openmm_main_pdb.topology, openmm_main_pdb.positions, open(path, "w"))
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
