
# Use scipy to solve linprog problem(s)
#
# Run like this:
#
#     python ./utils/linprog.py
#
# This is helpful for comparison with numa linprog results

import numpy as np

## Import linprog like this to use the pip-installed scipy version
#from scipy.optimize import linprog

# Import like this to use manually-extracted scipy functions with added debug
# prints 
from _linprog import linprog

#===================================================================

c = [-1, 4]
#c = [4, -1]
A = [[-3, 1], [1, 2]]
b = [6, 4]

x0_bounds = (None, None)
x1_bounds = (-3, None)

#res = linprog(c, A_ub=A, b_ub=b)
res = linprog(c, A_ub=A, b_ub=b, bounds=[x0_bounds, x1_bounds])

#res = scipy.optimize.linprog(c, A_ub=A, b_ub=b, bounds=[x0_bounds, x1_bounds])

print("===============================================")
print("scipy doc example")
print("res = ", res)
print("fun = ", res["fun"])
print("x = ", res["x"])
print("msg = ", res["message"])
exit(0)
#print("fun = ", res.fun)
#print("x   = ", res.x)
#print("msg = ", res.message)

##===================================================================
#
#A_ub = None
#b_ub = None
#bounds = None
#
#c = [2.8, 6.3, 10.8, -2.8, -6.3, -10.8]
#A_eq = [[-1, -1, -1, 0, 0, 0],
#        [0, 0, 0, 1, 1, 1],
#        [1, 0, 0, 1, 0, 0],
#        [0, 1, 0, 0, 1, 0],
#        [0, 0, 1, 0, 0, 1]]
#b_eq = [-0.5, 0.4, 0.3, 0.3, 0.3]
#
#res = linprog(c, A_ub, b_ub, A_eq, b_eq, method = "simplex")
#
#print("===============================================")
#print("enzo example")
#print("res = ", res)
#print("fun = ", res["fun"])
#print("x = ", res["x"])  # expect  [0.1875 1.25]
#print("msg = ", res["message"])
#exit(0)
#
##===================================================================

#	! MATLAB example with bounds (all constraint types)
#
#	c = [-1.d0, -1.d0/3]  ! `f` in MATLAB's linprog()
#	a_ub = transpose(reshape([ &  ! `A` in MATLAB
#			1.0, 1.0, &
#			1.0, 0.25, &
#			1.0, -1.0, &
#			-0.25, -1.0, &
#			-1.0, -1.0, &
#			-1.0, 1.0  &
#		] &
#		, [2, 6] &
#	))
#	b_ub = [2, 1, 2, 1, -1, 2]
#	a_eq = reshape([1.0, 0.25] , [1, 2])
#	b_eq = [0.5]
#	lb = [-1.0, -0.5]
#	ub = [1.5, 1.25]
#
#	expect = [0.1875d0, 1.25d0]

c = [-1.0, -1.0/3]
a_ub = [
    [1.0, 1.0],
    [1.0, 0.25],
    [1.0, -1.0],
    [-0.25, -1.0],
    [-1.0, -1.0],
    [-1.0, 1.0],
]
b_ub = [2, 1, 2, 1, -1, 2]
a_eq = [[1.0, 0.25]]
b_eq = [0.5]
#lb = [-1.0, -0.5]
#ub = [1.5, 1.25]
x0_bounds = [-1.0, 1.5]
x1_bounds = [-0.5, 1.25]

res = linprog(c, A_ub=a_ub, b_ub=b_ub, 
              A_eq=a_eq, b_eq=b_eq, 
              bounds=[x0_bounds, x1_bounds]
)

print("===============================================")
print("MATLAB doc example")
print("res = ", res)
print("fun = ", res["fun"])
print("x = ", res["x"])  # expect  [0.1875 1.25]
print("msg = ", res["message"])
#exit(0)

#===================================================================

#  >>> objective = ('minimize', '4x_1 + 1x_2')
#  >>> constraints = ['3x_1 + 1x_2 = 3', '4x_1 + 3x_2 >= 6', '1x_1 + 2x_2 <= 4']
#  >>> Lp_system = Simplex(num_vars=2, constraints=constraints, objective_function=objective)
#  >>> print(Lp_system.solution)
#  {'x_2': Fraction(6, 5), 'x_1': Fraction(3, 5)}
#  
#  >>> print(Lp_system.optimize_val)
#  17/5

c = [4, 1]
Aeq = [[3, 1]]
beq = [3]
Ale = [[-4, -3], [1, 2]]
ble = [-6, 4]

#x0_bounds = (None, None)
#x1_bounds = (-3, None)

#res = linprog(c, A_ub=A, b_ub=b, bounds=[x0_bounds, x1_bounds])
res = linprog(c, A_ub=Ale, b_ub=ble, A_eq = Aeq, b_eq = beq)

print("===============================================")
print("khalibartan doc example")
#print("fun = ", res.fun)
print("x   = ", res.x)
print("msg = ", res.message)

#===================================================================

#def _linprog_simplex(c, c0, A, b, callback, postsolve_args,
#                     maxiter=1000, tol=1e-9, disp=False, bland=False,
#                     **unknown_options):
#    """
#    Minimize a linear objective function subject to linear equality and
#    non-negativity constraints using the two phase simplex method.
#    Linear programming is intended to solve problems of the following form:
#
#    Minimize::
#
#        c @ x
#
#    Subject to::
#
#        A @ x == b
#            x >= 0


# Standard form example from:
#
#     https://www3.nd.edu/%7Edgalvin1/30210/30210_F07/presentations/converting.pdf

# Note that "standard form" is quite different from "canonical form"

#c = [-1, -4, -1]
#A = [[-3, 1, 2], [1, 2, 5]]
#b = [6, 4]

c = np.array([8, 6, 0, 0])
A = np.array([
        [1, 1, -1, 0],
        [-0.05, 0.07, 0, 1],
])
b = np.array([
    1,
    0
])

x0b = (0, None)  # bounds
x1b = (0, None)
x2b = (0, None)
x3b = (0, None)

res = linprog(c, A_eq = A, b_eq = b, bounds = [x0b, x1b, x2b, x3b])

print("===============================================")
print("standard form example 1")
print("fun = ", res.fun)
print("x   = ", res.x)
print("msg = ", res.message)

#scipy.optimize._linprog_simplex(c)
#def _linprog_simplex(c, c0, A, b, callback, postsolve_args,
#                     maxiter=1000, tol=1e-9, disp=False, bland=False,
#                     **unknown_options):

#def _linprog_simplex(c, c0, A, b, callback, postsolve_args,
#                     maxiter=1000, tol=1e-9, disp=False, bland=False,
#                     **unknown_options):

c0 = 0
res = _linprog_simplex._linprog_simplex(c, c0, A, b, None, None)

print("===============================================")
print("standard form example 1 using _linprog_simplex()")
print("res = ", res)
#print("fun = ", res.fun)
#print("x   = ", res.x)
#print("msg = ", res.message)

#=============================

