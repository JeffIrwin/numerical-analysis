import numpy as np
from warnings import warn
#from ._optimize import OptimizeResult, OptimizeWarning, _check_unknown_options
#from ._linprog_util import _postsolve

def _pivot_col(T, tol=1e-9, bland=False):
    #print("T[-1, :-1] = ", T[-1, :-1])
    #print("T[-1, :-1] >= -tol = ", T[-1, :-1] >= -tol)
    ma = np.ma.masked_where(T[-1, :-1] >= -tol, T[-1, :-1], copy=False)
    #print("ma = ", ma)
    if ma.count() == 0:
        return False, np.nan
    if bland:
        # ma.mask is sometimes 0d
        return True, np.nonzero(np.logical_not(np.atleast_1d(ma.mask)))[0][0]
    #print("_pivot_col = ", np.ma.nonzero(ma == ma.min())[0][0])
    #exit(0)
    return True, np.ma.nonzero(ma == ma.min())[0][0]


def _pivot_row(T, basis, pivcol, phase, tol=1e-9, bland=False):
    if phase == 1:
        k = 2
    else:
        k = 1
    ma = np.ma.masked_where(T[:-k, pivcol] <= tol, T[:-k, pivcol], copy=False)
    if ma.count() == 0:
        return False, np.nan
    mb = np.ma.masked_where(T[:-k, pivcol] <= tol, T[:-k, -1], copy=False)
    q = mb / ma
    print("q = ", q)
    print("q == q.min() = ", q == q.min())
    min_rows = np.ma.nonzero(q == q.min())[0]
    print("min_rows = ", min_rows)
    if bland:
        return True, min_rows[np.argmin(np.take(basis, min_rows))]
    print("_pivot_row = ", min_rows[0])
    #exit(0)
    return True, min_rows[0]


def _apply_pivot(T, basis, pivrow, pivcol, tol=1e-9):
    basis[pivrow] = pivcol
    pivval = T[pivrow, pivcol]
    T[pivrow] = T[pivrow] / pivval
    for irow in range(T.shape[0]):
        if irow != pivrow:
            T[irow] = T[irow] - T[pivrow] * T[irow, pivcol]

    # The selected pivot should never lead to a pivot value less than the tol.
    if np.isclose(pivval, tol, atol=0, rtol=1e4):
        message = (
            f"The pivot operation produces a pivot value of:{pivval: .1e}, "
            "which is only slightly greater than the specified "
            f"tolerance{tol: .1e}. This may lead to issues regarding the "
            "numerical stability of the simplex method. "
            "Removing redundant constraints, changing the pivot strategy "
            "via Bland's rule or increasing the tolerance may "
            "help reduce the issue.")
        warn(message, OptimizeWarning, stacklevel=5)


def _solve_simplex(T, n, basis, callback, postsolve_args,
                   maxiter=1000, tol=1e-9, phase=2, bland=False, nit0=0,
                   ):
    nit = nit0
    status = 0
    message = ''
    complete = False

    if phase == 1:
        m = T.shape[1]-2
    elif phase == 2:
        m = T.shape[1]-1
    else:
        raise ValueError("Argument 'phase' to _solve_simplex must be 1 or 2")

    if phase == 2:
        # Check if any artificial variables are still in the basis.
        # If yes, check if any coefficients from this row and a column
        # corresponding to one of the non-artificial variable is non-zero.
        # If found, pivot at this term. If not, start phase 2.
        # Do this for all artificial variables in the basis.
        # Ref: "An Introduction to Linear Programming and Game Theory"
        # by Paul R. Thie, Gerard E. Keough, 3rd Ed,
        # Chapter 3.7 Redundant Systems (pag 102)
        for pivrow in [row for row in range(basis.size)
                       if basis[row] > T.shape[1] - 2]:
            non_zero_row = [col for col in range(T.shape[1] - 1)
                            if abs(T[pivrow, col]) > tol]
            if len(non_zero_row) > 0:
                pivcol = non_zero_row[0]
                _apply_pivot(T, basis, pivrow, pivcol, tol)
                nit += 1

    if len(basis[:m]) == 0:
        solution = np.empty(T.shape[1] - 1, dtype=np.float64)
    else:
        solution = np.empty(max(T.shape[1] - 1, max(basis[:m]) + 1),
                            dtype=np.float64)

    print("bland = ", bland)
    while not complete:
        print("nit = ", nit)

        print("T = ")
        print(T)

        # Find the pivot column
        pivcol_found, pivcol = _pivot_col(T, tol, bland)
        print("pivcol = ", pivcol)

        if not pivcol_found:
            pivcol = np.nan
            pivrow = np.nan
            status = 0
            complete = True
        else:
            # Find the pivot row
            pivrow_found, pivrow = _pivot_row(T, basis, pivcol, phase, tol, bland)
            print("pivrow = ", pivrow)
            if not pivrow_found:
                status = 3
                complete = True

        if not complete:
            if nit >= maxiter:
                # Iteration limit exceeded
                status = 1
                complete = True
            else:
                _apply_pivot(T, basis, pivrow, pivcol, tol)
                nit += 1
    return nit, status


def _linprog_simplex(c, c0, A, b, callback, postsolve_args,
                     maxiter=1000, tol=1e-9, disp=False, bland=False,
                     **unknown_options):
    #_check_unknown_options(unknown_options)

    status = 0
    messages = {0: "Optimization terminated successfully.",
                1: "Iteration limit reached.",
                2: "Optimization failed. Unable to find a feasible"
                   " starting point.",
                3: "Optimization failed. The problem appears to be unbounded.",
                4: "Optimization failed. Singular matrix encountered."}

    n, m = A.shape

    # All constraints must have b >= 0.
    is_negative_constraint = np.less(b, 0)
    A[is_negative_constraint] *= -1
    b[is_negative_constraint] *= -1

    # As all constraints are equality constraints the artificial variables
    # will also be basic variables.
    av = np.arange(n) + m
    basis = av.copy()

    # Format the phase one tableau by adding artificial variables and stacking
    # the constraints, the objective row and pseudo-objective row.
    row_constraints = np.hstack((A, np.eye(n), b[:, np.newaxis]))
    row_objective = np.hstack((c, np.zeros(n), c0))
    row_pseudo_objective = -row_constraints.sum(axis=0)
    row_pseudo_objective[av] = 0
    T = np.vstack((row_constraints, row_objective, row_pseudo_objective))

    np.set_printoptions(threshold=np.inf)  # Disable truncation

    print("T initial = ")
    print(T)

    nit1, status = _solve_simplex(T, n, basis, callback=callback,
                                  postsolve_args=postsolve_args,
                                  maxiter=maxiter, tol=tol, phase=1,
                                  bland=bland
                                  )

    print("T after phase 1= ")
    print(np.array2string(T))

    # if pseudo objective is zero, remove the last row from the tableau and
    # proceed to phase 2
    nit2 = nit1
    if abs(T[-1, -1]) < tol:
        # Remove the pseudo-objective row from the tableau
        T = T[:-1, :]
        # Remove the artificial variable columns from the tableau
        T = np.delete(T, av, 1)
    else:
        # Failure to find a feasible starting point
        status = 2
        messages[status] = (
            "Phase 1 of the simplex method failed to find a feasible "
            "solution. The pseudo-objective function evaluates to "
            f"{abs(T[-1, -1]):.1e} "
            f"which exceeds the required tolerance of {tol} for a solution to be "
            "considered 'close enough' to zero to be a basic solution. "
            "Consider increasing the tolerance to be greater than "
            f"{abs(T[-1, -1]):.1e}. "
            "If this tolerance is unacceptably large the problem may be "
            "infeasible."
        )

    if status == 0:
        # Phase 2
        nit2, status = _solve_simplex(T, n, basis, callback=callback,
                                      postsolve_args=postsolve_args,
                                      maxiter=maxiter, tol=tol, phase=2,
                                      bland=bland, nit0=nit1
                                      )

    solution = np.zeros(n + m)
    solution[basis[:n]] = T[:n, -1]
    x = solution[:m]

    return x, status, messages[status], int(nit2)
