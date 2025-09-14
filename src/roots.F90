
!===============================================================================

#include "panic.F90"

!> Module for high-level linear algebra.  Some other more basic routines can be
!> found in the module numa__blarg
module numa__roots

	use numa__core
	use numa__linalg
	use numa__utils

	implicit none

	interface newton_raphson
		procedure :: newton_raphson_1d
		procedure :: newton_raphson_nd
	end interface

contains

!===============================================================================

double precision function newton_raphson_1d(f, df, x0, maxiters, tol, iostat) result(x)
	! Find the root of a scalar function `f` with derivative `df` using Newton-Raphson
	procedure(fn_f64_to_f64) :: f, df
	double precision, optional, intent(in) :: x0
	integer, optional, intent(in) :: maxiters
	double precision, optional, intent(in) :: tol
	integer, optional, intent(out) :: iostat
	!********
	character(len = :), allocatable :: msg
	integer :: maxiters_, iter
	double precision :: tol_, fx, dfx, x0_
	logical :: converged

	if (present(iostat)) iostat = 0

	maxiters_ = 15
	if (present(maxiters)) maxiters_ = maxiters

	tol_ = 1.d-10
	if (present(tol)) tol_ = tol

	! TODO: maybe 1 should be default?  Derivative is sometimes 0 at 1 for many
	! fns.  Same for newton_raphson_nd()
	x0_ = 0.d0
	if (present(x0)) x0_ = x0

	x = x0_
	converged = .false.
	do iter = 1, maxiters_

		fx = f(x)
		if (abs(fx) < tol_) then
			converged = .true.
			exit
		end if

		!! Maybe add verbose option?
		!print *, "iter, fx = ", iter, fx

		dfx = df(x)
		x = x - fx / dfx

	end do

	if (.not. converged) then
		! One last fn eval
		fx = f(x)
		converged = abs(fx) < tol_
	end if
	if (.not. converged) then
		msg = "newton_raphson_1d() did not converge"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

end function newton_raphson_1d

!===============================================================================

function newton_raphson_nd(f, df, x0, nx, maxiters, tol, iostat) result(x)
	! Find the root of a scalar function using Newton-Raphson
	!
	! df(i,j) is the derivative of fx(i) wrt x(j) for fx = f(x)
	procedure(fn_vec_f64_to_vec_f64) :: f
	procedure(fn_vec_f64_to_mat_f64) :: df
	double precision, optional, intent(in) :: x0(:)
	integer, optional, intent(in) :: nx
	integer, optional, intent(in) :: maxiters
	double precision, optional, intent(in) :: tol
	integer, optional, intent(out) :: iostat

	double precision, allocatable :: x(:)
	!********
	character(len = :), allocatable :: msg
	integer :: maxiters_, iter, io
	double precision :: tol_
	double precision, allocatable :: fx(:), dfx(:,:), x0_(:)
	logical :: converged

	if (present(iostat)) iostat = 0

	maxiters_ = 15
	if (present(maxiters)) maxiters_ = maxiters

	tol_ = 1.d-10
	if (present(tol)) tol_ = tol

	if (.not. present(x0) .and. .not. present(nx)) then
		! This is a bit of a weird Fortran problem where usually everything can
		! be allocatable and size agnostic, but we can't know the size of the
		! problem before we call the fn `f`, and we also can't call the function
		! to determine its output size without having a properly allocated input
		! variable
		msg = "either initial guess `x0` or system size `nx` is required as an argument for newton_raphson_nd()"
		call PANIC(msg, present(iostat))
		iostat = 2
		return
	end if

	if (present(x0)) then
		x0_ = x0
	else
		x0_ = zeros(nx)
	end if

	x = x0_
	converged = .false.
	do iter = 1, maxiters_

		fx = f(x)
		if (norm2(fx) < tol_) then
			converged = .true.
			exit
		end if

		!! Maybe add verbose option?
		!print *, "iter, norm(fx) = ", iter, norm2(fx)

		dfx = df(x)

		call lu_invmul(dfx, fx, io)
		if (io /= 0) then
			msg = "lu_invmul() failed in newton_raphson_nd()"
			call PANIC(msg, present(iostat))
			iostat = 3
			return
		end if

		x = x - fx

	end do

	if (.not. converged) then
		! One last fn eval
		fx = f(x)
		converged = norm2(fx) < tol_
	end if
	if (.not. converged) then
		msg = "newton_raphson_nd() did not converge"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

end function newton_raphson_nd

!===============================================================================

end module numa__roots

