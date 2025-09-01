
#include "panic.F90"

!> Module for linear programming, i.e. linear optimization
module numa__linprog

	!use numa__core
	use numa__blarg

	implicit none

	integer, parameter :: &
		LINPROG_SIMPLEX = 1, &
		LINPROG_REVISED_SIMPLEX = 2

	private :: linprog_solve_simplex

contains

!===============================================================================

subroutine linprog_get_abc(c, a_ub, b_ub, a_eq, b_eq, lbs, ubs, a, b)
	! Convert a general linprog problem to the standard form taken by
	! `linprog_std()`
	!
	! This is based on scipy:
	!
	!     https://github.com/scipy/scipy/blob/5f112b048b552e8788027dec7e95fb112daeeec1/scipy/optimize/_linprog_util.py#L1030

	use ieee_arithmetic
	use numa__blarg
	use numa__utils

	double precision, allocatable, intent(inout) :: c(:)
	double precision, allocatable, intent(inout) :: a_ub(:,:)
	double precision, allocatable, intent(inout) :: b_ub(:)
	double precision, intent(inout) :: a_eq(:,:)
	double precision, intent(inout) :: b_eq(:)  ! could be intent in
	double precision, intent(inout) :: lbs(:)
	double precision, intent(inout) :: ubs(:)   ! could be intent in

	double precision, allocatable, intent(out) :: a(:,:)
	double precision, allocatable, intent(out) :: b(:)

	!********
	double precision, allocatable :: ub_newub(:), a1(:,:), a2(:,:), lb_shift(:)
	integer :: m_ub, n_ub, n_bounds, n_free
	integer, allocatable :: i_nolb(:), i_newub(:), shape_(:), i_free(:), &
		i_shift(:)
	logical, allocatable :: lb_none(:), ub_none(:), lb_some(:), ub_some(:), &
		l_nolb_someub(:), l_free(:)

	!print *, "starting linprog_get_abc()"
	!print *, "a_ub = "
	!print "("//to_str(size(a_ub,2))//"es14.4)", transpose(a_ub)

	m_ub = size(a_ub, 1)
	n_ub = size(a_ub, 2)
	!print *, "m_ub, n_ub = ", m_ub, n_ub

	lb_none = lbs == -inf()
	ub_none = ubs == inf()
	ub_some = .not. ub_none

	!print *, "lb_none = ", lb_none
	!print *, "ub_none = ", ub_none

	! Unbounded below: substitute xi = -xi' (unbounded above)
	l_nolb_someub = lb_none .and. ub_some
	i_nolb = mask_to_index(l_nolb_someub)
	!print *, "l_nolb_someub = ", l_nolb_someub
	!print *, "i_nolb = ", i_nolb

	lbs(i_nolb) = -ubs(i_nolb)
	!print *, "lbs = ", lbs
	!print *, "ubs = ", ubs
	!print *, "c initial = ", c

	lb_none = lbs == -inf()
	ub_none = ubs == inf()
	lb_some = .not. lb_none
	ub_some = .not. ub_none
	c(i_nolb) = -c(i_nolb)
	!print *, "c after *-1 = ", c

	a_ub(:, i_nolb) = -a_ub(:, i_nolb)
	a_eq(:, i_nolb) = -a_eq(:, i_nolb)
	!print *, "a_eq = ", a_eq

	! Upper bound: add inequality constraint
	i_newub = mask_to_index(ub_some)
	ub_newub = ubs(i_newub)
	!print *, "i_newub = ", i_newub
	!print *, "ub_newub = ", ub_newub

	n_bounds = size(i_newub)
	if (n_bounds > 0) then
		! This condition seems unnecessary, but doesn't hurt

		shape_ = [n_bounds, size(a_ub, 2)]
		!print *, "shape_ = ", shape_

		! a_ub and b_ub change size here, thus they need to be allocatable
		a_ub = vstack(a_ub, zeros(shape_(1), shape_(2)))

		a_ub(m_ub+1:, i_newub) = eye(size(i_newub))

		!print *, "size(b_ub) = ", size(b_ub)
		b_ub = [b_ub, zeros(n_bounds)]
		!print *, "size(b_ub) = ", size(b_ub)

		!print *, "n_bounds = ", n_bounds
		!print *, "b_ub = ", b_ub
		!print *, "m_ub = ", m_ub

		b_ub(m_ub+1:) = ub_newub

	end if
	!call print_mat(a_ub, "a_ub after n bound > 0 = ")
	!print *, "b_ub = ", b_ub

	a1 = vstack(a_ub, a_eq)
	!print *, "a1 = "
	!print "("//to_str(size(a1,2))//"es14.4)", transpose(a1)

	!print *, "b_ub = ", b_ub
	!print *, "b_eq = ", b_eq
	!print *, "size(b_eq) = ", size(b_eq)
	!print *, "allocated(b_eq) = ", allocated(b_eq)

	b = [b_ub, b_eq]
	c = [c, zeros(size(a_ub, 1))]
	!print *, "b = ", b
	!print *, "c after cat = ", c

	! Unbounded: substitute xi = xi+ + xi-
	l_free = lb_none .and. ub_none
	i_free = mask_to_index(l_free)
	n_free = size(i_free)
	c = [c, zeros(n_free)]
	!print *, "c = ", c

	a1 = hstack(a1(:, :n_ub), -a1(:, i_free))
	!print *, "i_free = ", i_free

	c(n_ub+1: n_ub+n_free) = -c(i_free)

	!print *, "a1 = "
	!print "("//to_str(size(a1,2))//"es14.4)", transpose(a1)
	!print *, "n_ub, n_free = ", n_ub, n_free
	!print *, "c after i_free = ", c

	! Add slack variables
	a2 = vstack(eye(size(a_ub,1)), zeros(size(a_eq, 1), size(a_ub, 1)))

	a = hstack(a1, a2)
	!print *, "a = "
	!print "("//to_str(size(a,2))//"es14.4)", transpose(a)

	! Lower bound: substitute xi = xi' + lb
	i_shift = mask_to_index(lb_some)
	lb_shift = lbs(i_shift)
	!print *, "lb_some = ", lb_some
	!print *, "lb_shift = ", lb_shift

	!print *, "b = ", b
	!print *, "size(a) = ", size(a)
	b = b - sum(matmul(a(:, i_shift), diag(lb_shift)), 2)

	!print *, "****************"
	!print *, "linprog_get_abc() return values:"
	!print *, "a = "
	!print "("//to_str(size(a,2))//"es14.4)", transpose(a)
	!print *, "b = ", b
	!print *, "c = ", c
	!stop
	!print *, "****************"

end subroutine linprog_get_abc

!===============================================================================

function linprog &
	( &
		c, a_ub, b_ub, a_eq, b_eq, lb, ub, tol, iters, method, fval, iostat &
	) result(x)

	! Solve a general linear programming problem:
	!
	! Minimize:
	!
	!     dot_product(c, x)
	!
	! Subject to:
	!
	!     matmul(a_ub, x) <= b_ub
	!     matmul(a_eq, x) == b_eq  (optional)
	!     lb <= x <= ub            (lb default 0, ub default infinity)
	!
	! A wrapper fn could be added which avoids modifying a_ub, b_ub, a_eq, and
	! b_eq, like inv() vs invert(), by doing extra copying

	use ieee_arithmetic

	use numa__blarg
	use numa__utils

	double precision, allocatable, intent(in) :: c(:)
	double precision, allocatable, intent(inout) :: a_ub(:,:), b_ub(:)
	double precision, optional, allocatable, intent(inout) :: a_eq(:,:), b_eq(:)
	double precision, optional, allocatable, intent(in) :: lb(:), ub(:)

	double precision, optional, intent(in) :: tol
	integer, optional, intent(in) :: iters
	integer, optional, intent(in) :: method
	double precision, optional, intent(out) :: fval
	integer, optional, intent(out) :: iostat
	double precision, allocatable :: x(:)
	!********

	character(len = :), allocatable :: msg
	double precision :: lbi, ubi, tol_
	double precision, allocatable :: a(:,:), b(:), a_eq_(:,:), b_eq_(:), &
		lb_(:), ub_(:), c_(:), lbs0(:), ubs0(:)
	integer :: i, nx, n_unbounded, io, iters_, method_

	if (present(iostat)) iostat = 0

	tol_ = 1.d-10
	if (present(tol)) tol_ = tol

	iters_ = 1000
	if (present(iters)) iters_ = iters

	method_ = LINPROG_SIMPLEX
	if (present(method)) method_ = method

	x = [0]

	!print *, repeat("=", 60)
	!print *, "starting linprog()"

	if (size(c) /= size(a_ub,2)) then
		msg = "size(c) does not match size(a_ub,2) in linprog()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if
	if (present(a_eq)) then
		if (size(c) /= size(a_eq,2)) then
			msg = "size(c) does not match size(a_eq,2) in linprog()"
			call PANIC(msg, present(iostat))
			iostat = 2
			return
		end if
	end if
	if (present(lb)) then
		if (size(c) /= size(lb)) then
			msg = "size(c) does not match size(lb) in linprog()"
			call PANIC(msg, present(iostat))
			iostat = 3
			return
		end if
	end if
	if (present(ub)) then
		if (size(c) /= size(ub)) then
			msg = "size(c) does not match size(ub) in linprog()"
			call PANIC(msg, present(iostat))
			iostat = 4
			return
		end if
	end if
	if (size(a_ub,1) /= size(b_ub)) then
		msg = "size(a_ub,1) does not match size(b_ub) in linprog()"
		call PANIC(msg, present(iostat))
		iostat = 5
		return
	end if
	if (present(a_eq) .and. present(b_eq)) then
		if (size(a_eq,1) /= size(b_eq)) then
			msg = "size(a_eq,1) does not match size(b_eq) in linprog()"
			call PANIC(msg, present(iostat))
			iostat = 6
			return
		end if
	end if
	if (present(a_eq) .neqv. present(b_eq)) then
		msg = "`a_eq` must be present if and only if `b_eq` is present in linprog()"
		call PANIC(msg, present(iostat))
		iostat = 7
		return
	end if

	if (present(a_eq)) then
		a_eq_ = a_eq
	else
		allocate(a_eq_(0,0))
	end if
	if (present(b_eq)) then
		!print *, "b_eq present"
		b_eq_ = b_eq
	else
		!print *, "b_eq not present"
		allocate(b_eq_(0))
	end if

	if (present(lb)) then
		lb_ = lb
	else
		lb_ = zeros(size(c))
	end if
	if (present(ub)) then
		ub_ = ub
	else
		allocate(ub_(size(c)))
		ub_ = inf()
	end if
	!print *, "ub_ = ", ub_

	nx = size(c)

	! Make modifiable copies
	c_ = c

	! Copy backups
	lbs0 = lb_
	ubs0 = ub_

	call linprog_get_abc(c_, a_ub, b_ub, a_eq_, b_eq_, lb_, ub_, a, b)

	select case (method_)
	case (LINPROG_SIMPLEX)
		x = linprog_std(c_, a, b, tol_, iters_, io)
		if (io /= 0) then
			msg = "linprog_std() failed in linprog()"
			call PANIC(msg, present(iostat))
			iostat = 8
			return
		end if

	case (LINPROG_REVISED_SIMPLEX)
		!print *, "NOT IMPLEMENTED"
		x = linprog_rs(c_, a, b, tol_, iters_, io)
		if (io /= 0) then
			msg = "linprog_rs() failed in linprog()"
			call PANIC(msg, present(iostat))
			iostat = 10
			return
		end if

	case default
		msg = "bad method in linprog()"
		call PANIC(msg, present(iostat))
		iostat = 9
		return

	end select
	!print *, "x with slack = ", x

	!********
	! Post-solve

	! Undo variable substitutions of linprog_get_abc()
	n_unbounded = 0
	do i = 1, nx
		lbi = lbs0(i)
		ubi = ubs0(i)
		!print *, "lbi, ubi = ", lbi, ubi

		if (lbi == -inf() .and. ubi == inf()) then
			!print *, "unbounded"
			n_unbounded = n_unbounded + 1
			x(i) = x(i) - x(nx + n_unbounded)
		else
			if (lbi == -inf()) then
				x(i) = ubi - x(i)
			else
				x(i) = x(i) + lbi
			end if
		end if
	end do
	!print *, "nx = ", nx

	! Remove slack vars from x solution
	x = x(:nx)

	if (present(fval)) fval = dot_product(x, c)

end function linprog

!===============================================================================

function linprog_rs(c, a, b, tol, iters, iostat) result(x)
	! Revised simplex method
	!
	! Minimize:
	!
	!     dot_product(c, x)
	!
	! Subject to:
	!
	!     matmul(a, x) == b
	!     all(x >= 0)
	!
	! Note that *standard* form is different from *canonical* form

	use numa__blarg
	use numa__utils

	double precision, intent(in) :: c(:)
	double precision, intent(inout) :: a(:,:)
	double precision, intent(inout) :: b(:)

	double precision, optional, intent(in) :: tol
	integer, optional, intent(in) :: iters
	integer, optional, intent(out) :: iostat
	double precision, allocatable :: x(:)
	!********

	character(len = :), allocatable :: msg
	double precision, parameter :: c0 = 0  ! opt arg in scipy but never changed here
	double precision :: tol_
	double precision, allocatable :: r(:)
	!double precision, allocatable :: row_constraints(:,:), t(:,:)
	!double precision, allocatable :: row_objective(:)
	!double precision, allocatable :: row_pseudo_objective(:), solution(:)
	integer :: i, m, n, io, iters_, i1(1)
	!integer, allocatable :: av(:), basis(:)
	integer, allocatable :: basis(:), cols(:), rows(:), ineg(:), &
		nonzero_constraints(:)
	logical, allocatable :: l_tofix(:)

	if (present(iostat)) iostat = 0
	x = [0]

	!********
	! Generate auxiliary problem

	!def _generate_auxiliary_problem(A, b, x0, tol):
	!    print("Starting _generate_auxiliary_problem()")
	!    status = 0
	!    m, n = A.shape

	m = size(a, 1)
	n = size(a, 2)

	!    if x0 is not None:
	!        x = x0
	!    else:
	!        x = np.zeros(n)
	x = zeros(n)

	!    r = b - A@x  # residual; this must be all zeros for feasibility
	r = b - matmul(a, x)
	print *, "r = ", r

	!    A[r < 0] = -A[r < 0]  # express problem with RHS positive for trivial BFS
	!    b[r < 0] = -b[r < 0]  # to the auxiliary problem
	!    r[r < 0] *= -1

	!call print_mat(a, "a before neg = ")
	ineg = mask_to_index(r < 0)
	a(ineg, :) = -a(ineg, :)
	b(ineg) = -b(ineg)
	r(ineg) = -r(ineg)
	call print_mat(a, "a after neg = ")

	!    # Rows which we will need to find a trivial way to zero.
	!    # This should just be the rows where there is a nonzero residual.
	!    # But then we would not necessarily have a column singleton in every row.
	!    # This makes it difficult to find an initial basis.
	!    if x0 is None:
	!        nonzero_constraints = np.arange(m)
	!    else:
	!        nonzero_constraints = np.where(r > tol)[0]

	nonzero_constraints = [(i, i = 1, m)]
	print *, "nonzero_constraints = ", nonzero_constraints

	!    # these are (at least some of) the initial basis columns
	!    basis = np.where(np.abs(x) > tol)[0]

	!! TODO: not sure if `basis` is supposed to be all >tol locations or just
	!! the first
	!basis = mask_to_index(abs(x) > tol)
	i1 = findloc(abs(x) > tol, .true.)
	if (i1(1) > 0) then
		basis = [i1(1)]
	else
		allocate(basis(0))
	end if
	print *, "basis = ", basis

	!    if len(nonzero_constraints) == 0 and len(basis) <= m:  # already a BFS
	!        c = np.zeros(n)
	!        basis = _get_more_basis_columns(A, basis)
	!        return A, b, c, basis, x, status
	!    elif (len(nonzero_constraints) > m - len(basis) or
	!          np.any(x < 0)):  # can't get trivial BFS
	!        c = np.zeros(n)
	!        status = 6
	!        return A, b, c, basis, x, status

	if (size(nonzero_constraints) == 0 .and. size(basis) <= m) then
		print *, "TODO: already a basic feasible solution"
		! TODO
		stop
	else if (size(nonzero_constraints) > m - size(basis) .or. &
		any(x < 0)) then
		print *, "TODO: can't get trivial basis feasible solution"
		! TODO
		stop
	end if


	!********
	! Select singleton columns

	!    # chooses existing columns appropriate for inclusion in initial basis
	!    cols, rows = _select_singleton_columns(A, r)

	block

	double precision, allocatable :: columns(:,:), anz(:)
	integer, allocatable :: column_indices(:), row_indices(:), idx2(:,:), &
		nonzero_rows(:), nonzero_cols(:), same_sign(:), &
		unique_row_indices(:), first_columns(:)

	!def _select_singleton_columns(A, b):
	!    print("Starting _select_singleton_columns()")

	! Find indices of all singleton columns and corresponding row indices

	!    # find indices of all singleton columns and corresponding row indices
	!    column_indices = np.nonzero(np.sum(np.abs(A) != 0, axis=0) == 1)[0]

	column_indices = mask_to_index( &
		count(abs(a) /= 0, dim = 1) == 1 &
	)
	print *, "column_indices = ", column_indices

	!    columns = A[:, column_indices]          # array of singleton columns

	columns = a(:, column_indices)  ! array of singleton columns
	call print_mat(columns, "columns = ")

	!    row_indices = np.zeros(len(column_indices), dtype=int)
	!    nonzero_rows, nonzero_columns = np.nonzero(columns)
	!    row_indices[nonzero_columns] = nonzero_rows   # corresponding row indices

	row_indices = zeros_i32(size(column_indices))
	idx2 = mask_to_index(columns /= 0)
	nonzero_rows = idx2(1,:)
	nonzero_cols = idx2(2,:)
	call print_mat(idx2, "idx2 = ")
	print *, "nonzero_rows = ", nonzero_rows
	print *, "nonzero_cols = ", nonzero_cols

	row_indices(nonzero_cols) = nonzero_rows
	print *, "row_indices = ", row_indices

	!    # keep only singletons with entries that have same sign as RHS
	!    # this is necessary because all elements of BFS must be non-negative
	!    same_sign = A[row_indices, column_indices]*b[row_indices] >= 0
	!    column_indices = column_indices[same_sign][::-1]
	!    row_indices = row_indices[same_sign][::-1]

	!call print_mat(a(row_indices, column_indices), "a(rows, cols) = ")
	allocate(anz( size(idx2, 2) ))
	do i = 1, size(idx2, 2)
		anz(i) = columns(idx2(1,i), idx2(2,i))
	end do
	print *, "anz = ", anz

	! Keep only singletons with entries that have the same sign as RHS
	same_sign = mask_to_index( &
		anz * b(row_indices) >= 0 &
		!matmul(a(row_indices, column_indices), diag(b(row_indices))) >= 0 &
	)
	print *, "same_sign = ", same_sign

	print *, "column_indices slice = ", column_indices(same_sign)
	column_indices = reverse(column_indices(same_sign))
	print *, "column_indices reversed = ", column_indices

	row_indices = reverse(row_indices(same_sign))
	print *, "row_indices = ", row_indices

	!    # Reversing the order so that steps below select rightmost columns
	!    # for initial basis, which will tend to be slack variables. (If the
	!    # guess corresponds with a basic feasible solution but a constraint
	!    # is not satisfied with the corresponding slack variable zero, the slack
	!    # variable must be basic.)

	!    # for each row, keep rightmost singleton column with an entry in that row
	!    unique_row_indices, first_columns = np.unique(row_indices,
	!                                                  return_index=True)
	!    return column_indices[first_columns], unique_row_indices
	!    cols, rows = _select_singleton_columns(A, r)

	! For each row, keep the rightmost singleton column_indices with an entry in
	! that row
	call lp_unique(row_indices, unique_row_indices, first_columns)

	! End select singleton columns
	cols = column_indices(first_columns)
	rows = unique_row_indices

	! Scipy actually reverses cols and rows compared to what I have, but I don't
	! think that actually matters.  Probably because np.unique() does an
	! expensive and maybe unnecessary sorting operation

	end block
	!********

	print *, "rows = ", rows
	print *, "cols = ", cols

	!    # find the rows we need to zero that we _can_ zero with column singletons
	!    i_tofix = np.isin(rows, nonzero_constraints)

	l_tofix = is_in_vec(rows, nonzero_constraints)
	print *, "l_tofix = ", l_tofix

	!    # these columns can't already be in the basis, though
	!    # we are going to add them to the basis and change the corresponding x val
	!    i_notinbasis = np.logical_not(np.isin(cols, basis))
	!    i_fix_without_aux = np.logical_and(i_tofix, i_notinbasis)
	!    rows = rows[i_fix_without_aux]
	!    cols = cols[i_fix_without_aux]

	!    # indices of the rows we can only zero with auxiliary variable
	!    # these rows will get a one in each auxiliary column
	!    arows = nonzero_constraints[np.logical_not(
	!                                np.isin(nonzero_constraints, rows))]
	!    n_aux = len(arows)
	!    acols = n + np.arange(n_aux)          # indices of auxiliary columns

	!    basis_ng = np.concatenate((cols, acols))   # basis columns not from guess
	!    basis_ng_rows = np.concatenate((rows, arows))  # rows we need to zero

	!    # add auxiliary singleton columns
	!    A = np.hstack((A, np.zeros((m, n_aux))))
	!    A[arows, acols] = 1

	!    # generate initial BFS
	!    x = np.concatenate((x, np.zeros(n_aux)))
	!    x[basis_ng] = r[basis_ng_rows]/A[basis_ng_rows, basis_ng]

	!    # generate costs to minimize infeasibility
	!    c = np.zeros(n_aux + n)
	!    c[acols] = 1

	!    # basis columns correspond with nonzeros in guess, those with column
	!    # singletons we used to zero remaining constraints, and any additional
	!    # columns to get a full set (m columns)
	!    basis = np.concatenate((basis, basis_ng))
	!    basis = _get_more_basis_columns(A, basis)  # add columns as needed

	!    return A, b, c, basis, x, status


	! End generate auxiliary problem
	!********

end function linprog_rs

!===============================================================================

function is_in_vec(needles, haystack)
	integer, intent(in) :: needles(:), haystack(:)
	logical, allocatable :: is_in_vec(:)

	integer :: i, nn, nh

	nn = size(needles)
	nh = size(haystack)
	allocate(is_in_vec(nn))
	do i = 1, nn

		! TODO: make this not O(n^2).  Maybe sort both? Or hash map or
		! pigeonhole.  Pigeonhole might be effective and not explode memory
		! usage for the revised simplex application -- inputs are expected to be
		! in the range of matrix sizes
		!
		! In the generate case, pigeonhole could use a *lot* of memory for a
		! wide range of possibly haystack values

		is_in_vec(i) = any(haystack == needles(i))
	end do

end function is_in_vec

!===============================================================================

!call lp_unique(row_indices, unique_row_indices, first_columns)
subroutine lp_unique(vec, vals, idxs)
	integer, intent(in) :: vec(:)
	integer, allocatable, intent(out) :: vals(:)
	integer, allocatable, intent(out) :: idxs(:)
	!********
	integer :: i, n, nu

	n = size(vec)

	! Assume descending.  From my first revised simplex example, this may be ok
	!print *, "vec = ", vec
	if (any(vec(1:n-1) - vec(2:n) < 0)) then
		print *, "Error: vec is not sorted descending in lp_unique()"
		stop
	end if

	allocate(vals(n), idxs(n))

	vals(1) = vec(1)
	idxs(1) = 1
	nu = 1
	do i = 2, n
		if (vec(i) /= vec(i-1)) then
			nu = nu + 1
			vals(nu) = vec(i)
			idxs(nu) = i
		end if
	end do

end subroutine lp_unique

!===============================================================================

function linprog_std(c, a, b, tol, iters, iostat) result(x)
	! Solve a linear programming problem in standard form (?):
	!
	! Minimize:
	!
	!     dot_product(c, x)
	!
	! Subject to:
	!
	!     matmul(a, x) == b
	!     all(x >= 0)
	!
	! Note that *standard* form is different from *canonical* form
	!
	! This implementation is based on `_linprog_simplex()` from scipy:
	!
	!     https://github.com/scipy/scipy/blob/7d48c99615028935614943007fe61ce361dddebf/scipy/optimize/_linprog_simplex.py#L438C5-L438C21

	use numa__blarg
	use numa__utils

	double precision, intent(in) :: c(:)
	double precision, intent(inout) :: a(:,:)
	double precision, intent(inout) :: b(:)

	double precision, optional, intent(in) :: tol
	integer, optional, intent(in) :: iters
	integer, optional, intent(out) :: iostat
	double precision, allocatable :: x(:)
	!********

	character(len = :), allocatable :: msg
	double precision, parameter :: c0 = 0  ! opt arg in scipy but never changed here
	double precision :: tol_
	double precision, allocatable :: row_constraints(:,:), t(:,:)
	double precision, allocatable :: row_objective(:)
	double precision, allocatable :: row_pseudo_objective(:), solution(:)
	integer :: i, m, n, io, iters_
	integer, allocatable :: av(:), basis(:)

	if (present(iostat)) iostat = 0
	x = [0]

	tol_ = 1.d-10
	if (present(tol)) tol_ = tol

	iters_ = 1000
	if (present(iters)) iters_ = iters

	if (size(c) /= size(a,2)) then
		msg = "size(c) does not match size(a,2) in linprog_std()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if
	if (size(a,1) /= size(b)) then
		msg = "size(a,1) does not match size(b) in linprog_std()"
		call PANIC(msg, present(iostat))
		iostat = 5
		return
	end if

	n = size(a, 1)
	m = size(a, 2)
	!print *, "n, m = ", n, m

	! All constraints must have b >= 0
	do i = 1, size(b)
		if (b(i) >= 0) cycle
		a(i,:) = -a(i,:)
		b(i) = -b(i)
	end do

	! As all constraints are equality constraints, the artificial variables `av`
	! will also be initial basic variables
	av = [(i, i = 1, n)] + m
	basis = av

	!print *, "av    = ", av
	!print *, "basis = ", basis

	! Format the phase one tableau by adding artificial variables and stacking
	! the constraints, the objective row and pseudo-objective row
	row_constraints = hstack(hstack(a, eye(n)), b)

	!print *, "row_constraints = "
	!print "("//to_str(size(row_constraints,2))//"es14.4)", transpose(row_constraints)

	row_objective = [c, zeros(n), c0]
	!print *, "row_objective = "
	!print "(*(es14.4))", row_objective

	row_pseudo_objective = -sum(row_constraints, 1)
	row_pseudo_objective(av) = 0
	!print *, "row_pseudo_objective = "
	!print "(*(es14.4))", row_pseudo_objective

	t = vstack(vstack(row_constraints, row_objective), row_pseudo_objective)

	!print *, "t initial = "
	!print "("//to_str(size(t,2))//"es14.4)", transpose(t)

	call linprog_solve_simplex(t, basis, tol_, iters_, phase = 1, iostat = io)
	if (io /= 0) then
		msg = "linprog_solve_simplex() phase 1 failed in linprog_std()"
		call PANIC(msg, present(iostat))
		iostat = 9
		return
	end if

	!print *, "t after phase 1 = "
	!print "("//to_str(size(t,2))//"es14.4)", transpose(t)

	!********

	if (abs(t(size(t,1), size(t,2))) >= tol_) then
		msg = "solution is infeasible in linprog_std()"
		call PANIC(msg, present(iostat))
		iostat = 8
		return
	end if

	! Remove the pseudo-objective row from the tableau
	t = t(1: size(t,1) - 1, :)

	!print *, "t after removing pseudo-objective = "
	!print "("//to_str(size(t,2))//"es14.4)", transpose(t)

	! Remove the artificial variables columns from the tableau
	call rm_cols(t, av)

	!print *, "t after removing artificial vars = "
	!print "("//to_str(size(t,2))//"es14.4)", transpose(t)

	!********

	call linprog_solve_simplex(t, basis, tol_, iters_, phase = 2, iostat = io)
	if (io /= 0) then
		msg = "linprog_solve_simplex() phase 2 failed in linprog_std()"
		call PANIC(msg, present(iostat))
		iostat = 9
		return
	end if

	solution = zeros(n + m)
	solution(basis(:n)) = t(:n, size(t,2))
	x = solution(:m)

end function linprog_std

!********

subroutine rm_cols(a, indices)
	! Remove columns at given `indices` from matrix `a`.  Shrink `a` on output
	use numa__utils
	double precision, allocatable, intent(inout) :: a(:,:)
	integer, intent(in) :: indices(:)
	!********
	logical, allocatable :: mask(:)

	! Use a temp array to track which indices to keep.  Alternatively you could
	! assume that `indices` are sorted (possibly with verification)
	allocate(mask(size(a,2)))
	mask = .false.
	mask(indices) = .true.
	mask = .not. mask

	!print *, "indices = ", indices
	!print *, "mask    = ", mask
	!print *, "mask_to_index = ", mask_to_index(mask)
	!call print_mat(a, "before removing cols")

	a = a(:, mask_to_index(mask))

	!call print_mat(a, "after  removing cols")

end subroutine rm_cols

!********

subroutine linprog_solve_simplex(t, basis, tol, iters, phase, iostat)
	! This is based on `_solve_simplex()` from scipy:
	!
	!     https://github.com/scipy/scipy/blob/7d48c99615028935614943007fe61ce361dddebf/scipy/optimize/_linprog_simplex.py#L232
	!
	use numa__utils
	double precision, intent(inout) :: t(:,:)
	integer, intent(inout) :: basis(:)
	double precision, intent(in) :: tol
	integer, intent(in) :: iters
	integer, intent(in) :: phase
	integer, optional, intent(out) :: iostat
	!********
	character(len = :), allocatable :: msg
	double precision :: pv
	double precision, allocatable :: q(:)
	integer :: i, m, nb, mb, nc, nr, ntrunc, kc, kr, i1(1), iter, ir, ic
	logical :: complete

	if (present(iostat)) iostat = 0

	!print *, "starting linprog_solve_simplex(), phase = ", to_str(phase)

	if (phase == 1) then
		m = size(t, 2) - 2
		ntrunc = 2
	else
		m = size(t, 2) - 1
		ntrunc = 1
	end if
	!print *, "m = ", m

	if (phase /= 1) then
		! Check if any artificial variables are still in the basis ...

		do ir = 1, size(basis)
			if (basis(ir) <= size(t,2) - 2) cycle

			do ic = 1, size(t,2) - 1
				if (abs(t(ir, ic)) <= tol) cycle

				!print *, "phase 2 pre-pivot"
				!stop

				! Apply pivot
				kr = ir
				kc = ic

				basis(kr) = kc
				pv = t(kr, kc)
				t(kr, :) = t(kr, :) / pv
				do i = 1, size(t, 1)
					if (i == kr) cycle
					t(i,:) = t(i,:) - t(kr,:) * t(i, kc)
				end do

				! Pivot on first non-cycled `ic` only, but then iterate to next
				! row `ir`
				exit

			end do
		end do

	end if

	nb = size(basis)
	mb = min(nb, m)

	nr = size(t, 1)
	nc = size(t, 2)

	complete = .false.
	iter = 0
	do while (.not. complete)
		iter = iter + 1

		!print *, "iter"
		!print *, "t = "
		!print "("//to_str(size(t,2))//"es14.4)", transpose(t)

		! Find the pivot column
		!
		! Beware numpy.masked_where() works exactly the opposite as sensible
		! Fortran mask :explode:
		!
		! Also this minloc requires -standard-semantics for Intel compilers
		i1 = minloc(t(nr, :nc-1), t(nr, :nc-1) < -tol .and. t(nr, :nc-1) /= 0)
		kc = i1(1)

		!print *, "c = ", t(nr, :nc-1)
		!print *, "m = ", t(nr, :nc-1) < -tol
		!print *, "i1 = ", i1
		!print *, "pivcol = ", kc

		if (kc == 0) then
			! Successful end of iteration
			complete = .true.
			cycle
		end if

		q = t(:nr-ntrunc, nc) / t(:nr-ntrunc, kc)
		!print *, "q = ", q
		i1 = minloc(q, t(:nr-ntrunc, kc) > tol)

		kr = i1(1)
		!print *, "pivrow = ", kr
		if (kr == 0) then
			msg = "pivot row not found in linprog_solve_simplex()"
			call PANIC(msg, present(iostat))
			iostat = 1
			return
		end if

		! Apply pivot
		basis(kr) = kc
		pv = t(kr, kc)
		t(kr, :) = t(kr, :) / pv
		do i = 1, size(t, 1)
			if (i == kr) cycle
			t(i,:) = t(i,:) - t(kr,:) * t(i, kc)
		end do
		if (abs(pv) <= tol) then
			write(*,*) YELLOW // "Warning" // COLOR_RESET // &
				": pivot value is nearly 0. The linprog problem may be ill-conditioned"
		end if

		!print *, ""
		!stop
		if (iter >= iters) then
			msg = "max iterations reached in linprog_solve_simplex()"
			call PANIC(msg, present(iostat))
			iostat = 2
			return
		end if
	end do

end subroutine linprog_solve_simplex

!===============================================================================

end module numa__linprog

