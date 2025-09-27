
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

	interface golden_search
		procedure :: golden_search_1d
		!procedure :: golden_search_nd
	end interface

contains

!===============================================================================

double precision function fzero(f, xmin, xmax, tol, iters)
	! Source:  https://github.com/dr-nikolai/FMM/blob/master/fmm/zeroin.f
	procedure(fn_f64_to_f64) :: f
	double precision :: xmin, xmax, tol
	integer, optional, intent(out) :: iters
	! This function subprogram is a slightly  modified  translation  of
	! the algol procedure zero  given in richard brent, algorithms for
	! minimization without derivatives, prentice-hall, inc. (1973).

	double precision, parameter :: eps = epsilon(eps)
	double precision :: a, b, c, d, e, fa, fb, fc, tol1, xm, p, q, r, s
	integer :: iter_
	logical :: do_begin, do_bisect

	! initialization
	a = xmin
	b = xmax
	fa = f(a)
	fb = f(b)
	do_begin = .true.
	iter_ = 0

	do
		iter_ = iter_ + 1
		!print *, "iter = ", iter_

		if (do_begin) then
			! `do_begin` could probably have a better name
			c = a
			fc = fa
			d = b - a
			e = d
		end if

		if (abs(fc) < abs(fb)) then
			a = b
			b = c
			c = a
			fa = fb
			fb = fc
			fc = fa
		end if

		! convergence test
		tol1 = 2.d0*eps*abs(b) + 0.5d0*tol
		xm = 0.5d0 * (c - b)

		if (abs(xm) <= tol1) exit
		if (fb == 0.d0) exit

		! is bisection necessary?
		do_bisect = .true.
		if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then

			! is quadratic interpolation possible?
			if (a /= c) then
				! inverse quadratic interpolation
				!print *, "quadratic interpolation"
				q = fa/fc
				r = fb/fc
				s = fb/fa
				p = s*(2.d0*xm*q*(q - r) - (b - a)*(r - 1.d0))
				q = (q - 1.d0)*(r - 1.d0)*(s - 1.d0)
			else
				! linear interpolation
				!print *, "linear interpolation"
				s = fb/fa
				p = 2.d0*xm*s
				q = 1.d0 - s
			end if

			! adjust signs
			if (p > 0.d0) q = -q
			p = abs(p)

			! is interpolation acceptable?
			do_bisect = 2.d0*p >= 3.d0*xm*q - abs(tol1*q) .or. p >= abs(0.5d0*e*q)
			if (.not. do_bisect) then
				e = d
				d = p/q
			end if
		end if

		if (do_bisect) then
			! bisection
			!print *, "bisect"
			d = xm
			e = d
		end if

		! complete step
		a = b
		fa = fb
		if (abs(d) > tol1) b = b + d
		if (abs(d) <= tol1) b = b + sign(tol1, xm)
		fb = f(b)

		do_begin = fb*(fc/abs(fc)) > 0

	end do

	fzero = b
	if (present(iters)) iters = iter_

end function fzero

!===============================================================================

double precision function bisect_root(f, xmin, xmax, maxiters, tol, iostat) result(x)
	! Find the root of a scalar function `f` with between `xmin` and `xmax`
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax
	integer, optional, intent(in) :: maxiters
	double precision, optional, intent(in) :: tol
	integer, optional, intent(out) :: iostat
	!********
	character(len = :), allocatable :: msg
	integer :: maxiters_, iter
	double precision :: tol_, fx, xlo, xhi, flo, fhi
	logical :: converged

	if (present(iostat)) iostat = 0

	if (xmin >= xmax) then
		msg = "xmin must be less than xmax in bisect_root()"
		call PANIC(msg, present(iostat))
		iostat = 2
		return
	end if

	maxiters_ = 25
	if (present(maxiters)) maxiters_ = maxiters

	tol_ = 1.d-9
	if (present(tol)) tol_ = tol

	xlo = xmin
	xhi = xmax
	flo = f(xlo)
	fhi = f(xhi)

	if (sign_(flo) == sign_(fhi)) then
		msg = "f(xmin) and f(xmax) have the same sign in bisect_root()"
		call PANIC(msg, present(iostat))
		iostat = 3
		return
	end if

	converged = .false.
	do iter = 1, maxiters_
		x = 0.5d0 * (xlo + xhi)
		fx = f(x)
		if (abs(fx) < tol_) then
			converged = .true.
			exit
		end if

		!print *, "iter, xlo, x, xhi = ", iter, xlo, x, xhi

		if (sign_(flo) == sign_(fx)) then
			flo = fx
			xlo = x
		else
			fhi = fx
			xhi = x
		end if
	end do

	if (.not. converged) then
		msg = "bisect_root() did not converge"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

end function bisect_root

!===============================================================================

double precision function newton_raphson_1d(f, df, x0, maxiters, tol, iostat) result(x)
	! Find the root of a scalar function `f` with derivative `df` using Newton-Raphson
	!
	! TODO: this tol is a y tol. Maybe take an x tol like fzero, but how can it
	! be checked? Either difference between x iterations, or estimate
	! xtol = ytol / (dy/dx) (or conservatively check both).  Same for bisect
	! (although checking is easier there).  X tol makes more sense because the
	! root-finding fn returns an X value, not Y
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

function broyden(f, df, x0, nx, maxiters, tol, iostat) result(x)
	! Find the root of a scalar function using Broyden's method
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
	double precision :: tol_, alpha
	double precision, allocatable :: fx(:), dfx_inv(:,:), x0_(:), dx(:), &
		fx0(:), dfx(:), xdir(:)
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
		msg = "either initial guess `x0` or system size `nx` is required as an argument for broyden()"
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

	!dfx_inv = inv(df(x))
	!dfx_inv = gauss_jordan(df(x))
	dfx_inv = df(x)
	!call invert(dfx_inv, io)
	call gauss_jordan(dfx_inv, iostat = io)  ! TODO: inv wrapper instead of gauss_jordan direct?
	if (io /= 0) then
		msg = "gauss_jordan() failed in broyden()"
		call PANIC(msg, present(iostat))
		iostat = 3
		return
	end if

	fx = f(x)
	do iter = 1, maxiters_

		if (norm2(fx) < tol_) then
			converged = .true.
			exit
		end if
		print *, "iter, norm(fx) = ", iter, norm2(fx)

		xdir = matmul(dfx_inv, fx)

		alpha = line_search(f, x, xdir, 0.d0, 2.d0)
		!alpha = line_search(f, x, xdir, 0.d0, 1.d0)
		!!alpha = line_search(f, x, xdir, -0.1d0, 2.d0)
		!!alpha = line_search(f, x, xdir, -1.0d0, 3.d0)

		print *, "alpha = ", alpha

		!if (alpha == 0) then
		if (norm2(f(x - alpha * xdir)) >= norm2(f(x))) then
			dfx_inv = df(x)
			call gauss_jordan(dfx_inv, iostat = io)  ! TODO: inv wrapper instead of gauss_jordan direct?
			if (io /= 0) then
				msg = "gauss_jordan() failed in broyden()"
				call PANIC(msg, present(iostat))
				iostat = 3
				return
			end if
			alpha = 1.d0
			xdir = matmul(dfx_inv, fx)
		end if

		dx = alpha * xdir

		x = x - dx
		fx0 = fx
		fx = f(x)
		dfx = fx - fx0

		! TODO: optimize repeated mat ops
		dfx_inv = dfx_inv + &
			outer_product((dx - matmul(dfx_inv, dfx)) / &
			dot_product(dx, matmul(dfx_inv, dfx)), &
			matmul(dx, dfx_inv))

	end do

	if (.not. converged) then
		! One last fn eval
		fx = f(x)
		converged = norm2(fx) < tol_
	end if
	if (.not. converged) then
		msg = "broyden() did not converge"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

end function broyden

!===============================================================================

double precision function line_search(f, x, dx, amin, amax) result(alpha)
	! Find amin <= alpha <= amax, such that ||f(x - alpha * dx)|| is minimized
	procedure(fn_vec_f64_to_vec_f64) :: f
	double precision, intent(in) :: x(:), dx(:)
	double precision, intent(in) :: amin, amax

	double precision :: a, b, c, z, fb, fc, fa, fz
	integer :: i

	a = amin
	z = amax
	fa = norm2(f(x - a * dx))
	fz = norm2(f(x - z * dx))
	do i = 1, 10

		! TODO: golden section

		! Ternary search with equal-length intervals
		b = a + 1.d0 * (z - a) / 3.d0
		c = a + 2.d0 * (z - a) / 3.d0

		!! Ternary search close to midpoint
		!b = a + 0.49d0 * (z - a)
		!c = a + 0.51d0 * (z - a)

		fb = norm2(f(x - b * dx))
		fc = norm2(f(x - c * dx))

		if (fb <= fc) then
			z = c
			alpha = b
		else
			a = b
			alpha = c
		end if

	end do

	if (fa < fb .and. fa < fc) alpha = amin
	if (fz < fb .and. fz < fc) alpha = amax

end function line_search

!===============================================================================

double precision function golden_search_1d(f, xmin, xmax, xtol) result(xopt)
	! Find xopt in range [xmin, xmax] to minimize f(x)
	!
	! The golden section search is a line search that reduces function
	! evaluations, because the interior points of subintervals are shared from
	! one iteration to the next
	!
	! TODO: add ND version (like line_search()). Then probably delete
	! line_search()
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax
	double precision, optional, intent(in) :: xtol
	!********
	double precision, parameter :: INVPHI = (sqrt(5.d0) - 1.d0) / 2.d0
	double precision :: a, b, c, d, fa, fb, fc, fd, fxmin, fxmax, fopt

	! TODO: panic if xmin >= xmax

	! a < b < c < d
	a = xmin
	d = xmax
	b =  d - (d-a) * INVPHI
	c = a + (d-a) * INVPHI

	fa = f(a)
	fb = f(b)
	fc = f(c)
	fd = f(d)
	fxmin = fa
	fxmax = fd

	do while (d - a > xtol)
		if (fb <= fc) then
			xopt = b
			fopt = fb
			d = c
			c = b
			b =  d - (d-a) * INVPHI

			fd = fc
			fc = fb
			fb = f(b)
		else
			xopt = c
			fopt = fc
			a = b
			b = c
			c = a + (d-a) * INVPHI

			fa = fb
			fb = fc
			fc = f(c)
		end if
	end do

	if (fxmin < fopt) xopt = xmin
	if (fxmax < fopt) xopt = xmax

end function golden_search_1d

!===============================================================================

end module numa__roots

