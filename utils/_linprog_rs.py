# Author: Matt Haberland

import numpy as np
from numpy.linalg import LinAlgError

from scipy.linalg import solve
#from ._optimize import _check_unknown_options
from _bglu_dense import LU
from _bglu_dense import BGLU as BGLU
from _linprog_util import _postsolve
#from ._optimize import OptimizeResult


def _phase_one(A, b, x0, callback, postsolve_args, maxiter, tol, disp,
               maxupdate, mast, pivot):
    print("Starting _phase_one()")
    m, n = A.shape
    status = 0

    # generate auxiliary problem to get initial BFS
    A, b, c, basis, x, status = _generate_auxiliary_problem(A, b, x0, tol)

    if status == 6:
        residual = c.dot(x)
        iter_k = 0
        return x, basis, A, b, residual, status, iter_k

    # solve auxiliary problem
    phase_one_n = n
    iter_k = 0
    x, basis, status, iter_k = _phase_two(c, A, x, basis, callback,
                                          postsolve_args,
                                          maxiter, tol, disp,
                                          maxupdate, mast, pivot,
                                          iter_k, phase_one_n)

    # check for infeasibility
    residual = c.dot(x)
    if status == 0 and residual > tol:
        status = 2

    # drive artificial variables out of basis
    # TODO: test redundant row removal better
    # TODO: make solve more efficient with BGLU? This could take a while.
    keep_rows = np.ones(m, dtype=bool)
    for basis_column in basis[basis >= n]:
        B = A[:, basis]
        try:
            basis_finder = np.abs(solve(B, A))  # inefficient
            pertinent_row = np.argmax(basis_finder[:, basis_column])
            eligible_columns = np.ones(n, dtype=bool)
            eligible_columns[basis[basis < n]] = 0
            eligible_column_indices = np.where(eligible_columns)[0]
            index = np.argmax(basis_finder[:, :n]
                              [pertinent_row, eligible_columns])
            new_basis_column = eligible_column_indices[index]
            if basis_finder[pertinent_row, new_basis_column] < tol:
                keep_rows[pertinent_row] = False
            else:
                basis[basis == basis_column] = new_basis_column
        except LinAlgError:
            status = 4

    # form solution to original problem
    A = A[keep_rows, :n]
    basis = basis[keep_rows]
    x = x[:n]
    m = A.shape[0]
    return x, basis, A, b, residual, status, iter_k


def _get_more_basis_columns(A, basis):
    print("Starting _get_more_basis_columns()")
    m, n = A.shape

    # options for inclusion are those that aren't already in the basis
    a = np.arange(m+n)
    bl = np.zeros(len(a), dtype=bool)
    bl[basis] = 1
    options = a[~bl]
    options = options[options < n]  # and they have to be non-artificial
    print("options = ", options)

    # form basis matrix
    B = np.zeros((m, m))
    B[:, 0:len(basis)] = A[:, basis]
    print("B = ")
    print(B)

    if (basis.size > 0 and
            np.linalg.matrix_rank(B[:, :len(basis)]) < len(basis)):
        raise Exception("Basis has dependent columns")

    rank = 0  # just enter the loop
    for i in range(n):  # somewhat arbitrary, but we need another way out
        # permute the options, and take as many as needed
        new_basis = np.random.permutation(options)[:m-len(basis)]
        B[:, len(basis):] = A[:, new_basis]  # update the basis matrix
        rank = np.linalg.matrix_rank(B)      # check the rank
        if rank == m:
            break

    return np.concatenate((basis, new_basis))


def _generate_auxiliary_problem(A, b, x0, tol):
    print("Starting _generate_auxiliary_problem()")
    status = 0
    m, n = A.shape

    if x0 is not None:
        x = x0
    else:
        x = np.zeros(n)

    r = b - A@x  # residual; this must be all zeros for feasibility
    print("r = ", r)

    A[r < 0] = -A[r < 0]  # express problem with RHS positive for trivial BFS
    b[r < 0] = -b[r < 0]  # to the auxiliary problem
    r[r < 0] *= -1
    print("A after neg = ")
    print(A)
    print("r = ", r)

    # Rows which we will need to find a trivial way to zero.
    # This should just be the rows where there is a nonzero residual.
    # But then we would not necessarily have a column singleton in every row.
    # This makes it difficult to find an initial basis.
    if x0 is None:
        nonzero_constraints = np.arange(m)
    else:
        nonzero_constraints = np.where(r > tol)[0]
    print("nonzero_constraints = ", nonzero_constraints)

    # these are (at least some of) the initial basis columns
    basis = np.where(np.abs(x) > tol)[0]
    print("basis = ")
    print(basis)

    if len(nonzero_constraints) == 0 and len(basis) <= m:  # already a BFS
        c = np.zeros(n)
        basis = _get_more_basis_columns(A, basis)
        print("basis after getting more = ")
        print(basis)
        return A, b, c, basis, x, status
    elif (len(nonzero_constraints) > m - len(basis) or
          np.any(x < 0)):  # can't get trivial BFS
        c = np.zeros(n)
        status = 6
        print("status = 6")
        return A, b, c, basis, x, status

    # chooses existing columns appropriate for inclusion in initial basis
    cols, rows = _select_singleton_columns(A, r)
    print("rows = ", rows)
    print("cols = ", cols)
    #exit(0)

    # find the rows we need to zero that we _can_ zero with column singletons
    i_tofix = np.isin(rows, nonzero_constraints)
    print("rows = ", rows)
    print("nonzero_constraints = ", nonzero_constraints)
    print("i_tofix = ", i_tofix)

    # these columns can't already be in the basis, though
    # we are going to add them to the basis and change the corresponding x val
    i_notinbasis = np.logical_not(np.isin(cols, basis))
    i_fix_without_aux = np.logical_and(i_tofix, i_notinbasis)
    rows = rows[i_fix_without_aux]
    cols = cols[i_fix_without_aux]

    print("i_notinbasis = ", i_notinbasis)

    # indices of the rows we can only zero with auxiliary variable
    # these rows will get a one in each auxiliary column
    arows = nonzero_constraints[np.logical_not(
                                np.isin(nonzero_constraints, rows))]
    n_aux = len(arows)
    acols = n + np.arange(n_aux)          # indices of auxiliary columns

    print("arows = ", arows)

    basis_ng = np.concatenate((cols, acols))   # basis columns not from guess
    basis_ng_rows = np.concatenate((rows, arows))  # rows we need to zero

    # add auxiliary singleton columns
    A = np.hstack((A, np.zeros((m, n_aux))))
    A[arows, acols] = 1

    # generate initial BFS
    x = np.concatenate((x, np.zeros(n_aux)))
    x[basis_ng] = r[basis_ng_rows]/A[basis_ng_rows, basis_ng]
    print("x = ", x)

    # generate costs to minimize infeasibility
    c = np.zeros(n_aux + n)
    c[acols] = 1
    print("c = ", c)

    # basis columns correspond with nonzeros in guess, those with column
    # singletons we used to zero remaining constraints, and any additional
    # columns to get a full set (m columns)
    basis = np.concatenate((basis, basis_ng))
    print("basis after ng cat = ")
    print(basis)
    basis = _get_more_basis_columns(A, basis)  # add columns as needed
    print("basis after getting more = ")
    print(basis)

    print("Ending _generate_auxiliary_problem")
    print("=" * 60)

    return A, b, c, basis, x, status


def _select_singleton_columns(A, b):
    print("Starting _select_singleton_columns()")
    # find indices of all singleton columns and corresponding row indices
    column_indices = np.nonzero(np.sum(np.abs(A) != 0, axis=0) == 1)[0]
    print("column_indices = ", column_indices)

    columns = A[:, column_indices]          # array of singleton columns
    print("columns = ")
    print(columns)

    row_indices = np.zeros(len(column_indices), dtype=int)
    nonzero_rows, nonzero_columns = np.nonzero(columns)
    row_indices[nonzero_columns] = nonzero_rows   # corresponding row indices

    print("nonzero_rows    = ", nonzero_rows)
    print("nonzero_columns = ", nonzero_columns)
    print("row_indices     = ", row_indices)

    print("A slice = ")
    print(A[row_indices, column_indices])
    print("b slice = ")
    print(b[row_indices])
    print("mat = ")
    print(A[row_indices, column_indices]*b[row_indices])

    # keep only singletons with entries that have same sign as RHS
    # this is necessary because all elements of BFS must be non-negative
    same_sign = A[row_indices, column_indices]*b[row_indices] >= 0
    column_indices = column_indices[same_sign][::-1]
    row_indices = row_indices[same_sign][::-1]

    print("same_sign = ", same_sign)
    print("column_indices (reversed) = ", column_indices)
    print("row_indices = ", row_indices)

    # Reversing the order so that steps below select rightmost columns
    # for initial basis, which will tend to be slack variables. (If the
    # guess corresponds with a basic feasible solution but a constraint
    # is not satisfied with the corresponding slack variable zero, the slack
    # variable must be basic.)

    # for each row, keep rightmost singleton column with an entry in that row
    unique_row_indices, first_columns = np.unique(row_indices,
                                                  return_index=True)
    return column_indices[first_columns], unique_row_indices

def _select_enter_pivot(c_hat, bl, a, rule="bland", tol=1e-12):
    print("Starting _select_enter_pivot()")
    if rule.lower() == "mrc":  # index with minimum reduced cost
        return a[~bl][np.argmin(c_hat)]
    else:  # smallest index w/ negative reduced cost
        return a[~bl][c_hat < -tol][0]

def _phase_two(c, A, x, b, callback, postsolve_args, maxiter, tol, disp,
               maxupdate, mast, pivot, iteration=0, phase_one_n=None):
    print("Starting _phase_two()")
    m, n = A.shape
    status = 0
    a = np.arange(n)                    # indices of columns of A
    ab = np.arange(m)                   # indices of columns of B

    # BGLU works but I might want to port LU to Fortran instead first
    if False and maxupdate:
        # basis matrix factorization object; similar to B = A[:, b]
        B = BGLU(A, b, maxupdate, mast)
    else:
        B = LU(A, b)
    print("b = ", b)
    print("A = ")
    print(A)
    print("B.B = ")
    print(B.B)

    for iteration in range(iteration, maxiter):

        #if disp or callback is not None:
        #    _display_and_callback(phase_one_n, x, postsolve_args, status,
        #                          iteration, disp, callback)

        bl = np.zeros(len(a), dtype=bool)
        bl[b] = 1

        xb = x[b]       # basic variables
        cb = c[b]       # basic costs

        try:
            v = B.solve(cb, transposed=True)    # similar to v = solve(B.T, cb)
        except LinAlgError:
            status = 4
            break

        # TODO: cythonize?
        c_hat = c - v.dot(A)    # reduced cost
        c_hat = c_hat[~bl]
        # Above is much faster than:
        # N = A[:, ~bl]                 # slow!
        # c_hat = c[~bl] - v.T.dot(N)
        # Can we perform the multiplication only on the nonbasic columns?

        if np.all(c_hat >= -tol):  # all reduced costs positive -> terminate
            break

        j = _select_enter_pivot(c_hat, bl, a, rule=pivot, tol=tol)
        u = B.solve(A[:, j])        # similar to u = solve(B, A[:, j])

        i = u > tol                 # if none of the u are positive, unbounded
        if not np.any(i):
            status = 3
            break

        th = xb[i]/u[i]
        l = np.argmin(th)           # implicitly selects smallest subscript
        th_star = th[l]             # step size

        x[b] = x[b] - th_star*u     # take step
        x[j] = th_star
        B.update(ab[i][l], j)       # modify basis
        b = B.b                     # similar to b[ab[i][l]] =

    else:
        # If the end of the for loop is reached (without a break statement),
        # then another step has been taken, so the iteration counter should
        # increment, info should be displayed, and callback should be called.
        iteration += 1
        status = 1
        #if disp or callback is not None:
        #    _display_and_callback(phase_one_n, x, postsolve_args, status,
        #                          iteration, disp, callback)

    return x, b, status, iteration


def _linprog_rs(c, c0, A, b, x0, callback, postsolve_args,
                maxiter=5000, tol=1e-12, disp=False,
                maxupdate=10, mast=False, pivot="mrc",
                **unknown_options):
    print("Starting _linprog_rs()")
    #_check_unknown_options(unknown_options)

    messages = ["Optimization terminated successfully.",
                "Iteration limit reached.",
                "The problem appears infeasible, as the phase one auxiliary "
                "problem terminated successfully with a residual of {0:.1e}, "
                "greater than the tolerance {1} required for the solution to "
                "be considered feasible. Consider increasing the tolerance to "
                "be greater than {0:.1e}. If this tolerance is unacceptably "
                "large, the problem is likely infeasible.",
                "The problem is unbounded, as the simplex algorithm found "
                "a basic feasible solution from which there is a direction "
                "with negative reduced cost in which all decision variables "
                "increase.",
                "Numerical difficulties encountered; consider trying "
                "method='interior-point'.",
                "Problems with no constraints are trivially solved; please "
                "turn presolve on.",
                "The guess x0 cannot be converted to a basic feasible "
                "solution. "
                ]

    if A.size == 0:  # address test_unbounded_below_no_presolve_corrected
        return np.zeros(c.shape), 5, messages[5], 0

    x, basis, A, b, residual, status, iteration = (
        _phase_one(A, b, x0, callback, postsolve_args,
                   maxiter, tol, disp, maxupdate, mast, pivot))

    if status == 0:
        x, basis, status, iteration = _phase_two(c, A, x, basis, callback,
                                                 postsolve_args,
                                                 maxiter, tol, disp,
                                                 maxupdate, mast, pivot,
                                                 iteration)

    return x, status, messages[status].format(residual, tol), iteration
