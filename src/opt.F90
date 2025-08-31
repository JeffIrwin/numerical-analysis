
#include "panic.F90"

!> Module for optimization.  See module numa__linprog for linear programming
module numa__opt

	use numa__core

	implicit none

contains

!===============================================================================

function nelder_mead(f, x0, x_tol, iters) result(x)
	! Solve an optimization problem by finding `x` in order to minimize 
	! f(x)

	use numa__utils
	procedure(fn_vec_f64_to_f64) :: f
	double precision, intent(in) :: x0(:)
	double precision, optional, intent(in) :: x_tol
	integer, optional, intent(in) :: iters

	double precision, allocatable :: x(:)
	!********

	double precision :: fr, fe, fc, x_tol_
	double precision, parameter :: alpha_ = 1, gamma_ = 2, rho_ = 0.5, sigma_ = 0.5
	double precision, allocatable :: xs(:,:), fs(:), xo(:), xr(:), &
		xe(:), xc(:)

	integer :: i, nb, nb1, iter, iters_
	integer, allocatable :: idx(:)

	logical :: converged

	x_tol_ = 1.d-3
	iters_ = 1000
	if (present(x_tol)) x_tol_ = x_tol
	if (present(iters)) iters_ = iters
	nb = size(x0)
	nb1 = nb + 1  ! number of simplex points

	! Initial simplex
	allocate(xs(nb, nb1))
	do i = 1, nb
		xs(:,i) = x0
		xs(i,i) = xs(i,i) + 1
	end do
	xs(:,nb1) = x0

	! Evaluate fn on initial simplex
	allocate(fs(nb1))
	do i = 1, nb1
		fs(i) = f(xs(:,i))
	end do
	!print *, "fs init = ", fs

	converged = .false.
	do iter = 1, iters_

		! Sort
		call sortidx_f64_1(fs, idx)
		!print *, "idx = ", idx

		fs = fs(idx)
		xs = xs(:, idx)
		!print *, "fs = ", fs

		!print *, "xs = "
		!print "(2es16.6)", xs

		! Note: this is a bad criterion because 1 and nb1 may not be
		! representative of the simplex size
		if (norm2(xs(:,1) - xs(:,nb1)) < x_tol_) then
			!print *, "iter = ", iter
			converged = .true.
			exit
		end if

		! Centroid
		xo = sum(xs(:, 1: nb), dim = 2) / nb
		!print *, "xo = ", xo

		! Reflect
		xr = xo + alpha_ * (xo - xs(:,nb1))
		fr = f(xr)

		if (fs(1) <= fr .and. fr < fs(nb)) then
			! Replace
			fs(nb1) = fr
			xs(:,nb1) = xr
			cycle
		end if

		if (fr < fs(1)) then
			! Expand
			xe = xo + gamma_ * (xr - xo)
			fe = f(xe)

			if (fe < fr) then
				fs(nb1) = fe
				xs(:,nb1) = xe
			else
				fs(nb1) = fr
				xs(:,nb1) = xr
			end if
			cycle
		end if

		if (fr < fs(nb1)) then
			! Contract outward
			xc = xo + rho_ * (xr - xo)
			fc = f(xc)

			if (fc < fr) then
				fs(nb1) = fc
				xs(:,nb1) = xc
				cycle
			end if

		else
			! Contract inward
			xc = xo + rho_ * (xs(:,nb1) - xo)
			fc = f(xc)

			if (fc < fs(nb1)) then
				fs(nb1) = fc
				xs(:,nb1) = xc
				cycle
			end if

		end if

		! Shrink
		do i = 2, nb1
			xs(:,i) = xs(:,1) + sigma_ * (xs(:,i) - xs(:,1))
			fs(i) = f(xs(:,i))
		end do

	end do

	if (.not. converged) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": nelder_mead() has not converged"
	end if
	x = xs(:,1)

end function nelder_mead

!===============================================================================

end module numa__opt
