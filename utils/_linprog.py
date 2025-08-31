import numpy as np

#from ._optimize import OptimizeResult, OptimizeWarning
from warnings import warn
#from ._linprog_highs import _linprog_highs
#from ._linprog_ip import _linprog_ip

#from ._linprog_simplex import _linprog_simplex
from _linprog_simplex import _linprog_simplex

#from ._linprog_rs import _linprog_rs
#from ._linprog_doc import (_linprog_highs_doc, _linprog_ip_doc,  # noqa: F401
#                           _linprog_rs_doc, _linprog_simplex_doc,
#                           _linprog_highs_ipm_doc, _linprog_highs_ds_doc)

#from ._linprog_util import (
from _linprog_util import (
    _parse_linprog, _presolve, _get_Abc, _LPProblem, _autoscale,
    _postsolve, _check_result, _display_summary)
from copy import deepcopy

__all__ = ['linprog', 'linprog_verbose_callback', 'linprog_terse_callback']

__docformat__ = "restructuredtext en"

LINPROG_METHODS = [
    'simplex', 'revised simplex', 'interior-point', 'highs', 'highs-ds', 'highs-ipm'
]


def linprog_verbose_callback(res):
    x = res['x']
    fun = res['fun']
    phase = res['phase']
    status = res['status']
    nit = res['nit']
    message = res['message']
    complete = res['complete']

    saved_printoptions = np.get_printoptions()
    np.set_printoptions(linewidth=500,
                        formatter={'float': lambda x: f"{x: 12.4f}"})
    if status:
        print('--------- Simplex Early Exit -------\n')
        print(f'The simplex method exited early with status {status:d}')
        print(message)
    elif complete:
        print('--------- Simplex Complete --------\n')
        print(f'Iterations required: {nit}')
    else:
        print(f'--------- Iteration {nit:d}  ---------\n')

    if nit > 0:
        if phase == 1:
            print('Current Pseudo-Objective Value:')
        else:
            print('Current Objective Value:')
        print('f = ', fun)
        print()
        print('Current Solution Vector:')
        print('x = ', x)
        print()

    np.set_printoptions(**saved_printoptions)


def linprog_terse_callback(res):
    nit = res['nit']
    x = res['x']

    if nit == 0:
        print("Iter:   X:")
    print(f"{nit: <5d}   ", end="")
    print(x)


def linprog(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None,
            bounds=(0, None),
            #method='highs',
            method='simplex',
            callback=None,
            options=None, x0=None, integrality=None):
    meth = method.lower()
    methods = {"highs", "highs-ds", "highs-ipm",
               "simplex", "revised simplex", "interior-point"}

    if meth not in methods:
        raise ValueError(f"Unknown solver '{method}'")

    #if x0 is not None and meth != "revised simplex":
    #    warning_message = "x0 is used only when method is 'revised simplex'. "
    #    warn(warning_message, OptimizeWarning, stacklevel=2)

    #if np.any(integrality) and not meth == "highs":
    #    integrality = None
    #    warning_message = ("Only `method='highs'` supports integer "
    #                       "constraints. Ignoring `integrality`.")
    #    warn(warning_message, OptimizeWarning, stacklevel=2)
    #elif np.any(integrality):
    #    integrality = np.broadcast_to(integrality, np.shape(c))
    #else:
    integrality = None

    print("bounds = ", bounds)

    lp = _LPProblem(c, A_ub, b_ub, A_eq, b_eq, bounds, x0, integrality)
    lp, solver_options = _parse_linprog(lp, options, meth)
    tol = solver_options.get('tol', 1e-9)

    ## Give unmodified problem to HiGHS
    #if meth.startswith('highs'):
    #    if callback is not None:
    #        raise NotImplementedError("HiGHS solvers do not support the "
    #                                  "callback interface.")
    #    highs_solvers = {'highs-ipm': 'ipm', 'highs-ds': 'simplex',
    #                     'highs': None}

    #    sol = _linprog_highs(lp, solver=highs_solvers[meth],
    #                         **solver_options)
    #    sol['status'], sol['message'] = (
    #        _check_result(sol['x'], sol['fun'], sol['status'], sol['slack'],
    #                      sol['con'], lp.bounds, tol, sol['message'],
    #                      integrality))
    #    sol['success'] = sol['status'] == 0
    #    return OptimizeResult(sol)


    #warn(f"`method='{meth}'` is deprecated and will be removed in SciPy "
    #     "1.11.0. Please use one of the HiGHS solvers (e.g. "
    #     "`method='highs'`) in new code.", DeprecationWarning, stacklevel=2)

    iteration = 0
    complete = False  # will become True if solved in presolve
    undo = []

    # Keep the original arrays to calculate slack/residuals for original
    # problem.
    lp_o = deepcopy(lp)

    # Solve trivial problem, eliminate variables, tighten bounds, etc.
    rr_method = solver_options.pop('rr_method', None)  # need to pop these;
    rr = solver_options.pop('rr', True)  # they're not passed to methods
    c0 = 0  # we might get a constant term in the objective
    if solver_options.pop('presolve', True):
        (lp, c0, x, undo, complete, status, message) = _presolve(lp, rr,
                                                                 rr_method,
                                                                 tol)

    C, b_scale = 1, 1  # for trivial unscaling if autoscale is not used
    postsolve_args = (lp_o._replace(bounds=lp.bounds), undo, C, b_scale)

    if not complete:
        A, b, c, c0, x0 = _get_Abc(lp, c0)
        if solver_options.pop('autoscale', False):
            A, b, c, x0, C, b_scale = _autoscale(A, b, c, x0)
            postsolve_args = postsolve_args[:-2] + (C, b_scale)

        if meth == 'simplex':
            x, status, message, iteration = _linprog_simplex(
                c, c0=c0, A=A, b=b, callback=callback,
                postsolve_args=postsolve_args, **solver_options)

            print("x with slack = ", x)

        #elif meth == 'interior-point':
        #    x, status, message, iteration = _linprog_ip(
        #        c, c0=c0, A=A, b=b, callback=callback,
        #        postsolve_args=postsolve_args, **solver_options)
        #elif meth == 'revised simplex':
        #    x, status, message, iteration = _linprog_rs(
        #        c, c0=c0, A=A, b=b, x0=x0, callback=callback,
        #        postsolve_args=postsolve_args, **solver_options)

    # Eliminate artificial variables, re-introduce presolved variables, etc.
    disp = solver_options.get('disp', False)

    x, fun, slack, con = _postsolve(x, postsolve_args, complete)

    status, message = _check_result(x, fun, status, slack, con, lp_o.bounds,
                                    tol, message, integrality)

    if disp:
        _display_summary(message, status, fun, iteration)

    sol = {
        'x': x,
        'fun': fun,
        'slack': slack,
        'con': con,
        'status': status,
        'message': message,
        'nit': iteration,
        'success': status == 0}

    #return OptimizeResult(sol)
    return sol
