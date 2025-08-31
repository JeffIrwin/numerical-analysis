
#include "panic.F90"

!> Module for linear programming, i.e. linear optimization
module numa__linprog

	!use numa__core
	use numa__blarg

	implicit none

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

function linprog(c, a_ub, b_ub, a_eq, b_eq, lb, ub, tol, iters, fval, iostat) result(x)
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
	double precision, optional, intent(out) :: fval
	integer, optional, intent(out) :: iostat
	double precision, allocatable :: x(:)
	!********

	character(len = :), allocatable :: msg
	double precision :: lbi, ubi, tol_
	double precision, allocatable :: a(:,:), b(:), a_eq_(:,:), b_eq_(:), &
		lb_(:), ub_(:), c_(:), lbs0(:), ubs0(:)
	integer :: i, nx, n_unbounded, io, iters_

	if (present(iostat)) iostat = 0

	tol_ = 1.d-10
	if (present(tol)) tol_ = tol

	iters_ = 1000
	if (present(iters)) iters_ = iters

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

	x = linprog_std(c_, a, b, tol_, iters_, io)
	if (io /= 0) then
		msg = "linprog_std() failed in linprog()"
		call PANIC(msg, present(iostat))
		iostat = 8
		return
	end if
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
	double precision :: pivval
	double precision, allocatable :: q(:)
	integer :: i, k, m, nb, mb, nc, nr, pivcol, pivrow, i1(1), iter, ir, ic
	logical :: complete

	if (present(iostat)) iostat = 0

	! TODO: rename overly long variables

	!print *, "starting linprog_solve_simplex(), phase = ", to_str(phase)

	if (phase == 1) then
		m = size(t, 2) - 2
		! TODO: rename `k` so I don't accidentally change it by using it as a
		! loop var e.g.
		k = 2
	else
		m = size(t, 2) - 1
		k = 1
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
				pivrow = ir
				pivcol = ic

				basis(pivrow) = pivcol
				pivval = t(pivrow, pivcol)
				t(pivrow, :) = t(pivrow, :) / pivval
				do i = 1, size(t, 1)
					if (i == pivrow) cycle
					t(i,:) = t(i,:) - t(pivrow,:) * t(i, pivcol)
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
		pivcol = i1(1)

		!print *, "c = ", t(nr, :nc-1)
		!print *, "m = ", t(nr, :nc-1) < -tol
		!print *, "i1 = ", i1
		!print *, "pivcol = ", pivcol

		if (pivcol == 0) then
			! Successful end of iteration
			complete = .true.
			cycle
		end if

		q = t(:nr-k, nc) / t(:nr-k, pivcol)
		!print *, "q = ", q
		i1 = minloc(q, t(:nr-k, pivcol) > tol)

		pivrow = i1(1)
		!print *, "pivrow = ", pivrow
		if (pivrow == 0) then
			msg = "pivot row not found in linprog_solve_simplex()"
			call PANIC(msg, present(iostat))
			iostat = 1
			return
		end if

		! Apply pivot
		basis(pivrow) = pivcol
		pivval = t(pivrow, pivcol)
		t(pivrow, :) = t(pivrow, :) / pivval
		do i = 1, size(t, 1)
			if (i == pivrow) cycle
			t(i,:) = t(i,:) - t(pivrow,:) * t(i, pivcol)
		end do
		! TODO: check if pivval is too close to tol
		!
		! "The pivot operation produces a pivot value of ..."

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

