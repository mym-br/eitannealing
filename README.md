# Electric Impedance Tomography Simulated Anneling Solver

## Instruction to compile pyeitgeo Python library

1. Download and extract Eigen 3.3.5 to folder `/eigen-3.3.5`
2. Install Python, Cython and the correct compiler for the installed Python version (on Windows see the list <https://wiki.python.org/moin/WindowsCompilers>)
3. To generate the Python library, run:
   ```bash
   $ python setup.py build_ext --inplace
   ```

### Wheel generations

After performing steps 1-2 of the compilation process, execute:

```bash
$ python -m build
```
