# Electric Impedance Tomography Simulated Anneling Solver

## Instruction to compile pyeitsolver pybind11 Python library

1. If using vcpkg, create a environment variable CMAKE_TOOLCHAIN_FILE with the path to vcpkg.cmake. In Powershell, it would be:
   ```bash
   $ $env:CMAKE_TOOLCHAIN_FILE = "X:/path_to_vcpkg/scripts/buildsystems/vcpkg.cmake"
   ```
2. As the FindPython CMake module does not work well with the vcpkg and pybind11 combination, you need to manually set the Python paths CMake variables `PYTHON_INCLUDE_DIR` and `PYTHON_LIBRARY`. In Powershell, it would be, for example:
   ```bash
   $env:SKBUILD_CMAKE_ARGS = "-DPYTHON_INCLUDE_DIR=X:/path_to_python/include;-DPYTHON_LIBRARY=X:/path_to_python/libs/python3XX.lib"
   ```
3. If the default CMake kit is not compatible with your Python, use the environment variable to set the CMake generator. In Powershell, it would be:
   ```bash
   $ $env:CMAKE_GENERATOR = "Visual Studio 17 2022"
   ```
4. Other useful option can be set with:
- If the default platform is not compatible:
   ```bash
   $ $env:CMAKE_GENERATOR_PLATFORM = "x64"
   ```
- If you want to be able to debug it you could set the build type:
   ```bash
   $ $env:SKBUILD_CMAKE_BUILD_TYPE = "Debug"
   ```
- If you want to specify a build folder (for incremental build):
   ```bash
   $ $env:SKBUILD_BUILD_DIR = "build"
   ```
5. Finally, to generate and install the Python library, run:
   ```bash
   $ pip install .
   ```
   - use option `-v` for verbose output

## Wheel generation

1. After performing optional steps 1-4 of the compilation process, then execute:
   ```bash
   $ pip install build scikit-build-core
   $ python -m build --wheel
   ```
