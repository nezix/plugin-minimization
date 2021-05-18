# Nanome - Minimization

Minimization arranges the selected molecule to try to find a geometry where inter-atomic forces are as close to zero as possible.

## Dependencies

[Docker](https://docs.docker.com/get-docker/)

When running outside of Docker:

This plugin requires `openmm`.

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
