import numpy as np
import pyeitsolver

print(pyeitsolver.is_initialized())
print(pyeitsolver.init("a", "b"))
print(pyeitsolver.is_initialized())
print(pyeitsolver.solve_forward_problem([1, 2]))
print(pyeitsolver.solve_forward_problem(np.array([2, 3])))
print(pyeitsolver.solve_full_forward_problem(np.array([4, 5])))
