
!===============================================================================

#include "panic.F90"

!> Module for high-level linear algebra.  Some other more basic routines can be
!> found in the module numa__blarg
module numa__roots

	use numa__core
	!use numa__linalg
	use numa__utils

	implicit none

	interface newton_raphson
		procedure :: newton_raphson_1d
	end interface

contains

!===============================================================================

double precision function newton_raphson_1d(f, df, x0, maxiters, tol, iostat) result(x)

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

end module numa__roots

