import pathlib

from ._core import init as init_cpp


def init(mesh_filename: str, currents_filename: str):
    if pathlib.Path(mesh_filename).exists() == False:
        raise FileNotFoundError(
            "Input mesh file not found. The value was: {}".format(mesh_filename)
        )
    if pathlib.Path(currents_filename).exists() == False:
        raise FileNotFoundError(
            "Input currents file not found. The value was: {}".format(currents_filename)
        )
    return init_cpp(mesh_filename, currents_filename)


init.__doc__ = init_cpp.__doc__
