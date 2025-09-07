
#include "panic.F90"

!> Module for linear programming, i.e. linear optimization
module numa__linprog

	use numa__blarg
	use numa__linalg
	use numa__utils
	use numa__rand

	implicit none

	integer, parameter :: &
		LINPROG_SIMPLEX_METH = 1, &
		LINPROG_REVISED_SIMPLEX_METH = 2

	private :: lp_solve_simplex, rs_phase_two, lp_get_abc, &
		lp_get_more_cols, lp_unique, lp_is_in_vec, rs_gen_aux

contains

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

	method_ = LINPROG_SIMPLEX_METH
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

	call lp_get_abc(c_, a_ub, b_ub, a_eq_, b_eq_, lb_, ub_, a, b)

	select case (method_)
	case (LINPROG_SIMPLEX_METH)
		x = linprog_simplex(c_, a, b, tol_, iters_, io)
		if (io /= 0) then
			msg = "linprog_simplex() failed in linprog()"
			call PANIC(msg, present(iostat))
			iostat = 8
			return
		end if

	case (LINPROG_REVISED_SIMPLEX_METH)
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

	! Undo variable substitutions of lp_get_abc()
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

subroutine lp_get_abc(c, a_ub, b_ub, a_eq, b_eq, lbs, ubs, a, b)
	! Convert a general linprog problem to the standard form taken by
	! `linprog_simplex()`
	!
	! This is based on scipy:
	!
	!     https://github.com/scipy/scipy/blob/5f112b048b552e8788027dec7e95fb112daeeec1/scipy/optimize/_linprog_util.py#L1030

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

	!print *, "starting lp_get_abc()"
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
	!print *, "lp_get_abc() return values:"
	!print *, "a = "
	!print "("//to_str(size(a,2))//"es14.4)", transpose(a)
	!print *, "b = ", b
	!print *, "c = ", c
	!stop
	!print *, "****************"

end subroutine lp_get_abc

!===============================================================================

subroutine rs_gen_aux(c, a, b, x, basis, tol, iostat)
	! Generate auxiliary problem for revised simplex
	double precision, allocatable, intent(inout) :: c(:)
	double precision, allocatable, intent(inout) :: a(:,:)
	double precision, intent(inout) :: b(:)
	double precision, allocatable, intent(inout) :: x(:)
	integer, allocatable, intent(out) :: basis(:)
	double precision, intent(in) :: tol
	integer, intent(out) :: iostat

	!********

	character(len = :), allocatable :: msg
	double precision, allocatable :: r(:)
	integer :: i, m, n, n_aux, io
	integer, allocatable :: cols(:), rows(:), ineg(:), &
		nonzero_constraints(:), i_fix_without_aux(:), arows(:), acols(:), &
		basis_ng(:), basis_ng_rows(:)
	logical, allocatable :: l_tofix(:), l_notinbasis(:), l_fix_without_aux(:)

	iostat = 0

	m = size(a, 1)
	n = size(a, 2)

	x = zeros(n)
	r = b - matmul(a, x)
	!print *, "r = ", r

	!call print_mat(a, "a before neg = ")
	ineg = mask_to_index(r < 0)
	a(ineg, :) = -a(ineg, :)
	b(ineg) = -b(ineg)
	r(ineg) = -r(ineg)
	!call print_mat(a, "a after neg = ")

	nonzero_constraints = range_i32(m)
	!print *, "nonzero_constraints = ", nonzero_constraints

	basis = mask_to_index(abs(x) > tol)
	!print *, "basis = ", basis

	if (size(nonzero_constraints) == 0 .and. size(basis) <= m) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": untested execution path. Auxiliary problem is already a basic" &
			// " feasible solution in revised simplex"

		c = zeros(n)
		basis = lp_get_more_cols(a, basis, io)
		if (io /= 0) then
			msg = "lp_get_more_cols() failed in revised simplex"
			call PANIC(msg, .true.)
			iostat = 11
			return
		end if
		return

	else if (size(nonzero_constraints) > m - size(basis) .or. &
		any(x < 0)) then

		msg = "problem is infeasible for trivial basis in linprog_rs()"
		call PANIC(msg, .true.)
		iostat = 1
		return
	end if

	!********
	! Select singleton columns
	block

	double precision, allocatable :: columns(:,:), anz(:)
	integer, allocatable :: column_indices(:), row_indices(:), idx2(:,:), &
		nonzero_rows(:), nonzero_cols(:), same_sign(:), &
		unique_row_indices(:), first_columns(:)

	! Find indices of all singleton columns and corresponding row indices
	column_indices = mask_to_index(count(abs(a) /= 0, dim = 1) == 1)
	!print *, "column_indices = ", column_indices

	columns = a(:, column_indices)  ! array of singleton columns
	!call print_mat(columns, "columns = ")

	row_indices = zeros_i32(size(column_indices))
	idx2 = mask_to_index(columns /= 0)
	nonzero_rows = idx2(1,:)
	nonzero_cols = idx2(2,:)
	!call print_mat(idx2, "idx2 = ")
	!print *, "nonzero_rows = ", nonzero_rows
	!print *, "nonzero_cols = ", nonzero_cols

	row_indices(nonzero_cols) = nonzero_rows
	!print *, "row_indices = ", row_indices
	!call print_mat(a(row_indices, column_indices), "a(rows, cols) = ")

	allocate(anz( size(idx2, 2) ))
	do i = 1, size(idx2, 2)
		anz(i) = columns(idx2(1,i), idx2(2,i))
	end do
	!print *, "anz = ", anz

	! Keep only singletons with entries that have the same sign as RHS
	same_sign = mask_to_index(anz * b(row_indices) >= 0)
	column_indices = reverse(column_indices(same_sign))
	row_indices = reverse(row_indices(same_sign))

	!print *, "same_sign = ", same_sign
	!!print *, "column_indices slice = ", column_indices(same_sign)
	!print *, "column_indices reversed = ", column_indices
	!print *, "row_indices = ", row_indices

	! For each row, keep the rightmost singleton column_indices with an entry in
	! that row
	call lp_unique(row_indices, unique_row_indices, first_columns, io)
	if (io /= 0) then
		msg = "lp_unique() failed in revised simplex"
		call PANIC(msg, .true.)
		iostat = 10
		return
	end if

	! End select singleton columns
	cols = column_indices(first_columns)
	rows = unique_row_indices

	! Scipy actually reverses cols and rows compared to what I have, but I don't
	! think that actually matters.  Probably because np.unique() does an
	! expensive and maybe unnecessary sorting operation

	end block
	!********

	!print *, "rows = ", rows
	!print *, "cols = ", cols

	l_tofix = lp_is_in_vec(rows, nonzero_constraints)
	!print *, "l_tofix = ", l_tofix

	l_notinbasis = .not. lp_is_in_vec(cols, basis)
	l_fix_without_aux = l_tofix .and. l_notinbasis
	i_fix_without_aux = mask_to_index(l_fix_without_aux)
	rows = rows(i_fix_without_aux)
	cols = cols(i_fix_without_aux)

	!print *, "l_notinbasis = ", l_notinbasis
	!print *, "l_fix_without_aux = ", l_fix_without_aux

	arows = nonzero_constraints( &
		mask_to_index(.not. lp_is_in_vec(nonzero_constraints, rows)) &
	)
	n_aux = size(arows)
	acols = n + range_i32(n_aux)

	!print *, "arows = ", arows
	!print *, "acols = ", acols

	basis_ng = [cols, acols]
	basis_ng_rows = [rows, arows]
	!print *, "basis_ng = ", basis_ng
	!print *, "basis_ng_rows = ", basis_ng_rows

	! Add auxiliary singleton columns
	a = hstack(a, zeros(m, n_aux))
	do i = 1, size(arows)
		a(arows(i), acols(i)) = 1
	end do

	! Generate initial basic feasible solution (BFS)
	x = [x, zeros(n_aux)]
	do i = 1, size(basis_ng)
		x(basis_ng(i)) = r(basis_ng_rows(i)) / a(basis_ng_rows(i), basis_ng(i))
	end do
	!print *, "x bfs = ", x

	! Generate costs to minimize infeasibility
	c = zeros(n_aux + n)
	c(acols) = 1
	!print *, "c = ", c

	basis = [basis, basis_ng]
	!print *, "basis after ng cat       = ", basis
	basis = lp_get_more_cols(a, basis, io)
	if (io /= 0) then
		msg = "lp_get_more_cols() failed in revised simplex"
		call PANIC(msg, .true.)
		iostat = 11
		return
	end if
	!print *, "basis after getting more = ", basis
	!print *, "x = ", x

	!print *, "ending rs_gen_aux()"
	!print *, "basis = ", basis
	!print *, "c = ", c
	!print *, repeat("=", 60)

end subroutine rs_gen_aux

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

	double precision, allocatable, intent(inout) :: c(:)
	double precision, allocatable, intent(inout) :: a(:,:)
	double precision, intent(inout) :: b(:)

	double precision, optional, intent(in) :: tol
	integer, optional, intent(in) :: iters
	integer, optional, intent(out) :: iostat
	double precision, allocatable :: x(:)
	!********

	character(len = :), allocatable :: msg
	double precision :: tol_, residual, dmax
	double precision, allocatable :: bb(:,:), basis_finder(:,:), c0(:)
	integer :: i, m, n, iters_, i1(1), ib, basis_column, &
		pertinent_row, new_basis_column, io
	integer, allocatable :: basis(:)
	logical, allocatable :: keep_rows(:), eligible_columns(:)

	if (present(iostat)) iostat = 0
	x = [0]

	tol_ = 1.d-10
	if (present(tol)) tol_ = tol

	iters_ = 1000
	if (present(iters)) iters_ = iters

	! Backup before overwriting with aux problem
	c0 = c
	!print *, "c = ", c
	!stop

	m = size(a, 1)
	n = size(a, 2)

	call rs_gen_aux(c, a, b, x, basis, tol_, io)
	!print *, "c = ", c
	!print *, "basis = ", basis
	!print *, "size a = ", size(a,1), size(a,2)
	!print *, "x = ", x
	if (io /= 0) then
		msg = "auxiliary generation failed in revised simplex"
		call PANIC(msg, present(iostat))
		iostat = 12
		return
	end if

	!********

	! Note: this is just phase 2 to solve the auxiliary problem inside phase 1.
	! We still have to do a bunch of other stuff and then the actual phase 2
	! later
	call rs_phase_two(c, a, x, basis, iters_, tol_, io)
	if (io /= 0) then
		msg = "auxiliary phase 2 failed in revised simplex"
		call PANIC(msg, present(iostat))
		iostat = 9
		return
	end if

	! Check for infeasibility
	residual = dot_product(c, x)
	if (residual > tol_) then
		msg = "auxiliary problem is infeasible in revised simplex"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	! Drive artificial variables out of basis
	keep_rows = trues(m)
	do ib = 1, size(basis)
		if (basis(ib) <= n) cycle

		basis_column = basis(ib)
		bb = a(:, basis)
		basis_finder = abs(invmul(bb, a, io))
		if (io /= 0) then
			msg = "invmul() failed in revised simplex"
			call PANIC(msg, present(iostat))
			iostat = 1
			return
		end if

		i1 = maxloc(basis_finder(:, basis_column))
		pertinent_row = i1(1)

		eligible_columns = trues(n)
		eligible_columns(basis(mask_to_index(basis < n))) = .false.

		new_basis_column = 0
		dmax = -huge(dmax)
		do i = 1, n
			if (.not. eligible_columns(i)) cycle
			if (basis_finder(pertinent_row, i) > dmax) then
				dmax = basis_finder(pertinent_row, i)
				new_basis_column = i
			end if
		end do
		!print *, "new_basis_column = ", new_basis_column

		if (basis_finder(pertinent_row, new_basis_column) < tol_) then
			keep_rows(pertinent_row) = .false.
		else
			basis(mask_to_index(basis == basis_column)) = new_basis_column
		end if

	end do
	!print *, "keep_rows = ", keep_rows

	! Form solution to original problem
	a = a(mask_to_index(keep_rows), :n)
	basis = basis(mask_to_index(keep_rows))
	x = x(:n)

	! Phase 2, for real this time
	c = c0  ! restore
	call rs_phase_two(c, a, x, basis, iters_, tol_, io)
	if (io /= 0) then
		msg = "phase 2 failed in revised simplex"
		call PANIC(msg, present(iostat))
		iostat = 9
		return
	end if

end function linprog_rs

!===============================================================================

subroutine rs_phase_two(c, aa, x, b, maxiter, tol, iostat)

	double precision, intent(inout) :: c(:), aa(:,:), x(:)
	integer, intent(inout) :: b(:)  ! basis
	integer, intent(in) :: maxiter
	double precision, intent(in) :: tol
	integer, intent(out) :: iostat
	!********
	character(len = :), allocatable :: msg
	double precision :: th_star
	double precision, allocatable :: bb(:,:), xb(:), cb(:), v(:), c_hat(:), &
		u(:), th(:)
	integer :: i, j, l, i1(1), m, n, iteration, io
	integer, allocatable :: a(:), ab(:), ipos(:)
	logical :: converged
	logical, allocatable :: bl(:)

	iostat = 0

	!print *, repeat("=", 60)
	!print *, "Starting rs_phase_two()"
	!print *, "b = ", b
	!print *, "x = ", x

	m = size(aa, 1)
	n = size(aa, 2)
	a = range_i32(n)
	ab = range_i32(m)

	bb = aa(:, b)
	!call print_mat(aa, "aa = ")
	!call print_mat(bb, "bb = ")

	converged = .false.
	do iteration = 1, maxiter

		bl = falses(size(a))
		bl(b) = .true.

		xb = x(b)
		cb = c(b)

		v = invmul(transpose(bb), cb, io)
		if (io /= 0) then
			msg = "invmul() failed in revised simplex phase 2"
			call PANIC(msg, .true.)
			iostat = 1
			return
		end if

		c_hat = c - matmul(v, aa)

		!print *, "c = ", c
		!print *, "cb = ", cb
		!print *, "v = ", v
		!print *, "c_hat = ", c_hat
		!!print *, "bl = ", bl
		!print *, "c_hat (unmasked) = ", c_hat
		!print *, "size(bl)    = ", size(bl)
		!print *, "size(c_hat) = ", size(c_hat)

		if (all(c_hat >= -tol .or. bl)) then
			!print *, "converged"
			converged = .true.
			exit
		end if

		! Select enter pivot `j`
		j = 0
		do i = 1, size(a)
			if (bl(i)) cycle
			if (c_hat(i) >= -tol) cycle
			j = i
			exit
		end do
		!print *, "j = ", j

		u = invmul(bb, aa(:,j), io)
		if (io /= 0) then
			msg = "invmul() failed in revised simplex phase 2"
			call PANIC(msg, .true.)
			iostat = 1
			return
		end if
		!print *, "u = ", u

		! If none of `u` are positive, unbounded
		ipos = mask_to_index(u > tol)
		if (size(ipos) <= 0) then
			msg = "problem is unbounded in revised simplex"
			call PANIC(msg, .true.)
			iostat = 1
			return
		end if

		th = xb(ipos) / u(ipos)
		i1 = minloc(th)
		l = i1(1)
		th_star = th(l)

		x(b) = x(b) - th_star * u
		x(j) = th_star

		!print *, "u > tol = ", u > tol
		!print *, "th = ", th
		!print *, "th_star = ", th_star

		! Update (modify) basis and basis matrix
		i = ab(ipos(l))
		b(i: m-1) = b(i+1: m)
		b(size(b)) = j
		bb = aa(:, b)
		!print *, "b = ", b

	end do
	if (.not. converged) then
		msg = "revised simplex phase 2 did not converge"
		call PANIC(msg, .true.)
		iostat = 1
		return
	end if

end subroutine rs_phase_two

!===============================================================================

function lp_get_more_cols(aa, basis, iostat) result(res)
	! Get more basis columns
	double precision, intent(in) :: aa(:,:)
	integer, intent(in) :: basis(:)
	integer, allocatable :: res(:)
	integer, intent(out) :: iostat
	!********
	character(len = :), allocatable :: msg
	double precision, allocatable :: b(:,:)
	integer :: i, m, n, rank_
	integer, allocatable :: a(:), options(:), perm(:), new_basis(:)
	logical, allocatable :: bl(:)

	m = size(aa, 1)
	n = size(aa, 2)

	a = range_i32(m+n)
	bl = falses(m+n)
	bl(basis) = .true.
	options = a(mask_to_index(.not. bl))
	options = options(mask_to_index(options <= n))
	!print *, "options = ", options

	! Form basis matrix
	b = zeros(m, m)
	b(:, 1: size(basis)) = aa(:, basis)
	!call print_mat(b, "b = ")

	if (size(basis) > 0 .and. &
		.not. is_full_rank(b(:, :size(basis)))) then

		msg = "basis has dependent columns"
		call PANIC(msg, .true.)
		iostat = 1
		return
	end if
	!print *, "ok"

	rank_ = 0
	do i = 1, n
		perm = rand_perm(size(options))
		new_basis = options(perm(: m - size(basis)))

		!print *, "new_basis = ", new_basis
		if (size(new_basis) > 0) then
			b(:, size(basis):) = aa(:, new_basis)
		end if

		!!if (i > 1) stop  ! unreachable in tests
		if (is_full_rank(b)) exit
	end do

	res = [basis, new_basis]

end function lp_get_more_cols

!===============================================================================

function lp_is_in_vec(needles, haystack)
	integer, intent(in) :: needles(:), haystack(:)
	logical, allocatable :: lp_is_in_vec(:)
	!********
	integer :: i, nn, nh, hmin, hmax
	logical, allocatable :: in_haystack(:)

	nn = size(needles)
	nh = size(haystack)

	if (nh <= 0) then
		!print *, "haystack empty"
		lp_is_in_vec = falses(nn)
		return
	end if

	! Pigeonhole array for O(nn + nh)
	hmin = minval(haystack)
	hmax = maxval(haystack)
	allocate(in_haystack(hmin: hmax))
	in_haystack = .false.
	in_haystack(haystack) = .true.

	allocate(lp_is_in_vec(nn))

	do i = 1, nn

		! There are several possible strategies to make this not O(n^2).  Maybe
		! sort both? Or hash map or pigeonhole.  Pigeonhole might be effective
		! and not explode memory usage for the revised simplex application --
		! inputs are expected to be in the range of matrix sizes
		!
		! In the general case, pigeonhole could use a *lot* of memory for a
		! wide range of possibly haystack values

		!! O(n^2)
		!lp_is_in_vec(i) = any(haystack == needles(i))

		! O(nn + nh)
		if (hmin <= needles(i) .and. needles(i) <= hmax) then
			lp_is_in_vec(i) = in_haystack(needles(i))
		else
			lp_is_in_vec(i) = .false.
		end if

	end do

end function lp_is_in_vec

!===============================================================================

subroutine lp_unique(vec, vals, idxs, iostat)
	integer, intent(in) :: vec(:)
	integer, allocatable, intent(out) :: vals(:)
	integer, allocatable, intent(out) :: idxs(:)
	integer, intent(out) :: iostat
	!********
	character(len = :), allocatable :: msg
	integer :: i, n, nu

	iostat = 0
	n = size(vec)

	! Assume vec is descending.  From my test cases, this seems ok
	!print *, "vec = ", vec
	if (any(vec(1:n-1) - vec(2:n) < 0)) then
		msg = "vec is not sorted descending in lp_unique()"
		call PANIC(msg, .true.)
		iostat = 1
		return
	end if

	allocate(vals(n), idxs(n))
	if (n <= 0) return

	vals(1) = vec(1)
	idxs(1) = 1
	nu = 1
	do i = 2, n
		if (vec(i) /= vec(i-1)) then
			nu = nu + 1
			vals(nu) = vec(i)
			idxs(nu) = i
		!else
		!	! Unreachable.  Why does scipy call unique() at all?
		!	print *, "non-unique elem"
		!	stop
		end if
	end do

	!! Not necessary, but helpful for comparison with scipy
	!vals = reverse(vals)
	!idxs = reverse(idxs)

end subroutine lp_unique

!===============================================================================

function linprog_simplex(c, a, b, tol, iters, iostat) result(x)
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
		msg = "size(c) does not match size(a,2) in linprog_simplex()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if
	if (size(a,1) /= size(b)) then
		msg = "size(a,1) does not match size(b) in linprog_simplex()"
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
	av = m + range_i32(n)
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

	call lp_solve_simplex(t, basis, tol_, iters_, phase = 1, iostat = io)
	if (io /= 0) then
		msg = "phase 1 failed in linprog_simplex()"
		call PANIC(msg, present(iostat))
		iostat = 9
		return
	end if

	!print *, "t after phase 1 = "
	!print "("//to_str(size(t,2))//"es14.4)", transpose(t)

	!********

	if (abs(t(size(t,1), size(t,2))) >= tol_) then
		msg = "problem is infeasible in linprog_simplex()"
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

	call lp_solve_simplex(t, basis, tol_, iters_, phase = 2, iostat = io)
	if (io /= 0) then
		msg = "phase 2 failed in linprog_simplex()"
		call PANIC(msg, present(iostat))
		iostat = 9
		return
	end if

	solution = zeros(n + m)
	solution(basis(:n)) = t(:n, size(t,2))
	x = solution(:m)

end function linprog_simplex

!********

subroutine rm_cols(a, indices)
	! Remove columns at given `indices` from matrix `a`.  Shrink `a` on output
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

subroutine lp_solve_simplex(t, basis, tol, iters, phase, iostat)
	! This is based on `_solve_simplex()` from scipy:
	!
	!     https://github.com/scipy/scipy/blob/7d48c99615028935614943007fe61ce361dddebf/scipy/optimize/_linprog_simplex.py#L232
	!
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

	!print *, "starting lp_solve_simplex(), phase = ", to_str(phase)

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
			msg = "pivot row not found in lp_solve_simplex()"
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
			msg = "max iterations reached in lp_solve_simplex()"
			call PANIC(msg, present(iostat))
			iostat = 2
			return
		end if
	end do

end subroutine lp_solve_simplex

!===============================================================================

end module numa__linprog

