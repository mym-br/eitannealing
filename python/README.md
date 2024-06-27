# Electric Impedance Tomography Simulated Anneling Solver

## Instruction to compile pyeitsolver pybind11 Python library

1. If using vcpkg, create a environment variable CMAKE_TOOLCHAIN_FILE with the path to vcpkg.cmake. In Powershell, it would be:
   ```bash
   $ $env:CMAKE_TOOLCHAIN_FILE = "X:/path_to_vcpkg/scripts/buildsystems/vcpkg.cmake"
   ```
2. If the default CMake kit is not compatible with your Python, use the environment variable to set the CMake generator. In Powershell, it would be:
   ```bash
   $ $env:CMAKE_GENERATOR = "Visual Studio 17 2022"
   ```
3. Finally, to generate and install the Python library, run:
   ```bash
   $ pip install .
   ```

### Wheel generations

After performing optional steps 1-2 of the compilation process, execute:

```bash
$ pip install build scikit-build-core
$ python -m build --wheel
```
