
#include "panic.F90"

!> Module for data fitting
module numa__fit

	use numa__core
	use numa__linalg

	implicit none

	interface polyval
		procedure :: polyval_sca
		procedure :: polyval_vec
	end interface

contains

!===============================================================================

function polyfit_lu(x, y, n, iostat) result(p)
	! Fit a polynomial using LU decomposition.  The function polyfit() should be
	! used instead, which avoid unnecessary matmul's, using QR decomposition
	! instead
	!
	! Polynomial coefficients `p` are in ascending powers, unlike MATLAB's polyfit()
	use numa__utils
	double precision, intent(in) :: x(:), y(:)
	integer, intent(in) :: n
	integer, optional, intent(out) :: iostat

	double precision, allocatable :: p(:)
	!********

	character(len = :), allocatable :: msg
	double precision, allocatable :: xx(:,:), xtx(:,:)
	integer :: i, nx, io

	if (present(iostat)) iostat = 0

	nx = size(x)

	if (nx /= size(y)) then
		msg = "size(x) does not match size(y) in polyfit_lu()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if
	if (n+1 > nx) then
		msg = "polynomial degree is too high for size of data in polyfit_lu()"
		call PANIC(msg, present(iostat))
		iostat = 2
		return
	end if
	if (n < 0) then
		msg = "polynomial degree is negative in polyfit_lu()"
		call PANIC(msg, present(iostat))
		iostat = 3
		return
	end if

	allocate(xx(nx, n+1))
	xx(:,1) = 1
	do i = 2, n+1
		xx(:,i) = x * xx(:, i-1)
	end do

	xtx = matmul(transpose(xx), xx)
	p   = matmul(transpose(xx), y)

	call lu_invmul(xtx, p, io)
	if (io /= 0) then
		msg = "lu_invmul() failed in polyfit_lu()"
		call PANIC(msg, present(iostat))
		iostat = 4
		return
	end if
	!print *, "p = ", p
	!print *, "y = ", y

end function polyfit_lu

!===============================================================================

function polyfit(x, y, n, iostat) result(p)
	! Polynomial coefficients `p` are in ascending powers, unlike MATLAB's polyfit()
	use numa__utils
	double precision, intent(in) :: x(:), y(:)
	integer, intent(in) :: n
	integer, optional, intent(out) :: iostat

	double precision, allocatable :: p(:)
	!********

	character(len = :), allocatable :: msg
	double precision, allocatable :: xx(:,:)
	integer :: i, nx, io

	if (present(iostat)) iostat = 0

	nx = size(x)

	if (nx /= size(y)) then
		msg = "size(x) does not match size(y) in polyfit()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if
	if (n+1 > nx) then
		msg = "polynomial degree is too high for size of data in polyfit()"
		call PANIC(msg, present(iostat))
		iostat = 2
		return
	end if
	if (n < 0) then
		msg = "polynomial degree is negative in polyfit()"
		call PANIC(msg, present(iostat))
		iostat = 3
		return
	end if

	allocate(xx(nx, n+1))
	xx(:,1) = 1
	do i = 2, n+1
		xx(:,i) = x * xx(:, i-1)
	end do
	!print *, "xx = "
	!!print "("//to_str(nx)//"es16.6)", xx
	!print "("//to_str(n+1)//"es16.6)", transpose(xx)

	p = qr_solve(xx, y, allow_rect = .true., iostat = io)
	if (io /= 0) then
		msg = "qr_solve() failed in polyfit()"
		call PANIC(msg, present(iostat))
		iostat = 4
		return
	end if

end function polyfit

!===============================================================================

function polyder(p) result(dp)
	! Return the polynomial coefficients `dp` for the derivative of polynomial `p`
	double precision, intent(in) :: p(:)
	double precision, allocatable :: dp(:)
	!********
	integer :: np
	np = size(p)

	dp = p(2:) * range_i32(1, np-1)
	!allocate(dp(np-1))
	!do i = 1, np-1
	!	dp(i) = p(i+1) * i
	!end do
end function polyder

!===============================================================================

function polyval_vec(p, x) result(y)
	! Evaluate a polynomial `p` with coefficients given in ascending powers
	double precision, intent(in) :: p(:), x(:)
	double precision, allocatable :: y(:)
	!********
	double precision, allocatable :: xpow(:)
	integer :: j, nx

	nx = size(x)
	allocate(y(nx))

	y = p(1)
	xpow = x
	do j = 2, size(p)
		y = y + p(j) * xpow
		xpow = xpow * x
	end do

end function polyval_vec

function polyval_sca(p, x) result(y)
	double precision, intent(in) :: p(:), x
	double precision, allocatable :: y
	!********
	double precision :: y1(1)
	y1 = polyval_vec(p, [x])
	y = y1(1)
end function polyval_sca
!===============================================================================

function gauss_newton(x, y, f, df, beta0, iters, iostat) result(beta)
	! Use the Gauss-Newton algorithm to find parameters `beta` to fit a function
	! `f` to data `x` and `y` with an initial guess `beta0`.  The function has a
	! gradient `df` == df/dx
	!
	! This requires knowledge of the analytic derivative of f.  If the
	! derivative is not easily known, a derivative-free optimization algorithm
	! could be used instead like nelder_mead_fit()

	use numa__utils
	double precision, intent(in) :: x(:), y(:)
	procedure(fn_f64_params_to_f64) :: f
	procedure(fn_f64_params_to_vec_f64) :: df
	double precision, intent(in) :: beta0(:)
	integer, intent(in) :: iters
	integer, optional, intent(out) :: iostat

	double precision, allocatable :: beta(:)
	!********

	character(len = :), allocatable :: msg
	double precision, allocatable :: res(:), jac(:,:), delta(:)
	integer :: i, nx, nb, iter, io

	if (present(iostat)) iostat = 0

	nx = size(x)
	nb = size(beta0)
	allocate(res(nx))
	allocate(jac(nx, nb))

	beta = beta0
	do iter = 1, iters

		do i = 1, nx
			! Evaluate residual
			res(i) = y(i) - f(x(i), beta)
		end do
		!print *, "res = ", res

		do i = 1, nx
			! Evaluate Jacobian
			jac(i,:) = -df(x(i), beta)
		end do
		!print *, "jac = "
		!print "(2es16.6)", transpose(jac)

		delta = qr_solve(jac, res, allow_rect = .true., iostat = io)
		if (io /= 0) then
			msg = "qr_solve() failed in gauss_newton()"
			call PANIC(msg, present(iostat))
			iostat = 1
			return
		end if

		beta = beta - delta

	end do

end function gauss_newton

!===============================================================================

function nelder_mead_fit(x, y, f, beta0, beta_tol, iters) result(beta)
	! Fit data `x` and `y` by finding `beta` in order to minimize 
	! norm2(y - f(x, beta))
	!
	! This is 1D in terms of `x` and `y` for now, i.e. `f` takes a scalar `x`
	! and returns a scalar `y`, but it could probably be generalized to
	! multidimensional fns, and more easily so than gauss_newton().  The
	! parameters `beta` of course can have any size/dimension
	!
	! Maybe this should have an f_tol option instead of or in addition too
	! beta_tol

	use numa__utils
	double precision, intent(in) :: x(:), y(:)
	procedure(fn_f64_params_to_f64) :: f
	double precision, intent(in) :: beta0(:)
	double precision, optional, intent(in) :: beta_tol
	integer, optional, intent(in) :: iters

	double precision, allocatable :: beta(:)
	!********

	double precision :: fr, fe, fc, beta_tol_
	double precision, parameter :: alpha_ = 1, gamma_ = 2, rho_ = 0.5, sigma_ = 0.5
	double precision, allocatable :: bs(:,:), fs(:), bo(:), br(:), &
		be(:), bc(:)

	integer :: i, nx, nb, nb1, iter, iters_
	integer, allocatable :: idx(:)

	logical :: converged

	beta_tol_ = 1.d-3
	iters_ = 1000
	if (present(beta_tol)) beta_tol_ = beta_tol
	if (present(iters)) iters_ = iters

	nx = size(x)
	nb = size(beta0)
	nb1 = nb + 1  ! number of simplex points

	! Initial simplex
	allocate(bs(nb, nb1))
	do i = 1, nb
		bs(:,i) = beta0
		bs(i,i) = bs(i,i) + 1
	end do
	bs(:,nb1) = beta0

	! Evaluate fn on initial simplex
	allocate(fs(nb1))
	do i = 1, nb1
		fs(i) = nm_eval_res(bs(:,i))
	end do
	!print *, "fs init = ", fs

	converged = .false.
	do iter = 1, iters_

		! Sort
		call sortidx_f64_1(fs, idx)
		!print *, "idx = ", idx

		fs = fs(idx)
		bs = bs(:, idx)
		!print *, "fs = ", fs

		!print *, "bs = "
		!print "(2es16.6)", bs

		! Note: this is a bad criterion because 1 and nb1 may not be
		! representative of the simplex size
		if (norm2(bs(:,1) - bs(:,nb1)) < beta_tol_) then
			!print *, "iter = ", iter
			converged = .true.
			exit
		end if

		! Centroid
		bo = sum(bs(:, 1: nb), dim = 2) / nb
		!print *, "bo = ", bo

		! Reflect
		br = bo + alpha_ * (bo - bs(:,nb1))
		fr = nm_eval_res(br)

		if (fs(1) <= fr .and. fr < fs(nb)) then
			! Replace
			fs(nb1) = fr
			bs(:,nb1) = br
			cycle
		end if

		if (fr < fs(1)) then
			! Expand
			be = bo + gamma_ * (br - bo)
			fe = nm_eval_res(be)

			if (fe < fr) then
				fs(nb1) = fe
				bs(:,nb1) = be
			else
				fs(nb1) = fr
				bs(:,nb1) = br
			end if
			cycle
		end if

		if (fr < fs(nb1)) then
			! Contract outward
			bc = bo + rho_ * (br - bo)
			fc = nm_eval_res(bc)

			if (fc < fr) then
				fs(nb1) = fc
				bs(:,nb1) = bc
				cycle
			end if

		else
			! Contract inward
			bc = bo + rho_ * (bs(:,nb1) - bo)
			fc = nm_eval_res(bc)

			if (fc < fs(nb1)) then
				fs(nb1) = fc
				bs(:,nb1) = bc
				cycle
			end if

		end if

		! Shrink
		do i = 2, nb1
			bs(:,i) = bs(:,1) + sigma_ * (bs(:,i) - bs(:,1))
			fs(i) = nm_eval_res(bs(:,i))
		end do

	end do

	if (.not. converged) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": nelder_mead_fit() has not converged"
	end if
	beta = bs(:,1)

	!--------------------------------
	contains

		double precision function nm_eval_res(beta_) result(res)
			! Evaluate the residual and sum its squares
			double precision, intent(in) :: beta_(:)
			integer :: k
			res = 0
			do k = 1, nx
				res = res + (y(k) - f(x(k), beta_)) ** 2
			end do
		end function nm_eval_res

end function nelder_mead_fit

!===============================================================================

end module numa__fit

