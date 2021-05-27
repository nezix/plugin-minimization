# Nanome - Minimization

Minimization arranges the selected molecule to try to find a geometry where inter-atomic forces are as close to zero as possible.

## Dependencies

[Docker](https://docs.docker.com/get-docker/)

When running outside of Docker:

This plugin requires `openmm`, `pdbfixer`, `openmmforcefields`.

## Behavior

This plugin receives a workspace from Nanome and finds not selected atoms that are close to the selection in a 7A radius (these atoms will be assigned a mass of 0 to fix them in space). It then converts these atoms into a PDB file with complexes separated in different models (thanks to the Modeller class). Then, the temporary PDB file is read to create an OpenMM system using one of the specified forcefields: Gaff, Amber, CHARMM, OpenFF, Smirnoff.
For general purpose forcefields (Gaff, OpenFF, Smirnoff), the plugin uses openmmforcefields and a template generator to create a system and minimize it.
The minimization process is done iteratively (OpenMM runs loops of 'updateRate' steps and sends the result to Nanome), until it converges (tolerance = 1e-5) or it reaches the number of specified steps (2500 by default). 

## Usage

To run Minimization in a Docker container:

```sh
$ cd docker
$ ./build.sh
$ ./deploy.sh -a <plugin_server_address> [optional args]
```

---

In Nanome:

- Activate Plugin
- Select atoms to minimize
- Click Run to immediately minimize or Advanced Settings to choose different options

## Development

To run Minimization with autoreload:

```sh
$ python3 -m pip install -r requirements.txt
$ python3 run.py -r -a <plugin_server_address> [optional args]
```

## License

MIT
