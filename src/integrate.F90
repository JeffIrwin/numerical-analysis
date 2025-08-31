
#include "panic.F90"

!> Module for integration and quadrature
module numa__integrate

	use numa__core

	implicit none

	private :: simpson_adapt_aux, gk15_aux, gk15i_aux

	! Tables of abscissas and weights are taken from Jacob Williams' quadpack,
	! BSD-3 license:
	!
	!     https://github.com/jacobwilliams/quadpack/blob/702abfd5f0acbdb51439695334347a4b3c0dc87a/src/quadpack_generic.F90#L5209
	!
	! See 3p/LICENSES for details
	integer, parameter :: NG_GK15 = 7, NK_GK15 = 15
	double precision, parameter :: WXG_GK15(*,*) = reshape([ &
		1.29484966168869693270611432679082018329d-1,  9.49107912342758524526189684047851262401d-1, &
		1.29484966168869693270611432679082018329d-1, -9.49107912342758524526189684047851262401d-1, &
		2.79705391489276667901467771423779582487d-1,  7.41531185599394439863864773280788407074d-1, &
		2.79705391489276667901467771423779582487d-1, -7.41531185599394439863864773280788407074d-1, &
		3.81830050505118944950369775488975133878d-1,  4.05845151377397166906606412076961463347d-1, &
		3.81830050505118944950369775488975133878d-1, -4.05845151377397166906606412076961463347d-1, &
		4.17959183673469387755102040816326530612d-1,  0.00000000000000000000000000000000000000d0   &
	], [2, NG_GK15])
	double precision, parameter :: WXK_GK15(*,*) = reshape([ &
		6.30920926299785532907006631892042866651d-2,  9.49107912342758524526189684047851262401d-1, &
		6.30920926299785532907006631892042866651d-2, -9.49107912342758524526189684047851262401d-1, &
		1.40653259715525918745189590510237920400d-1,  7.41531185599394439863864773280788407074d-1, &
		1.40653259715525918745189590510237920400d-1, -7.41531185599394439863864773280788407074d-1, &
		1.90350578064785409913256402421013682826d-1,  4.05845151377397166906606412076961463347d-1, &
		1.90350578064785409913256402421013682826d-1, -4.05845151377397166906606412076961463347d-1, &
		2.09482141084727828012999174891714263698d-1,  0.00000000000000000000000000000000000000d0 , &
		2.29353220105292249637320080589695919936d-2,  9.91455371120812639206854697526328516642d-1, &
		2.29353220105292249637320080589695919936d-2, -9.91455371120812639206854697526328516642d-1, &
		1.04790010322250183839876322541518017444d-1,  8.64864423359769072789712788640926201211d-1, &
		1.04790010322250183839876322541518017444d-1, -8.64864423359769072789712788640926201211d-1, &
		1.69004726639267902826583426598550284106d-1,  5.86087235467691130294144838258729598437d-1, &
		1.69004726639267902826583426598550284106d-1, -5.86087235467691130294144838258729598437d-1, &
		2.04432940075298892414161999234649084717d-1,  2.07784955007898467600689403773244913480d-1, &
		2.04432940075298892414161999234649084717d-1, -2.07784955007898467600689403773244913480d-1  &
	], [2, NK_GK15])

contains

!===============================================================================

double precision function simpson_integrator_vals(yi, xmin, xmax) result(area)
	! Integrate an array of function values `yi` from xmin to xmax using
	! Simpson's 1/3 and 3/8 rules, depending on whether the size of `yi` is even
	! or odd

	! You could make other value-array-based integrators for higher-order
	! methods like Milne's or Weddle's rules, but more cases for the number of
	! points in the last subintervals would be required

	double precision, intent(in) :: yi(:)
	double precision, intent(in) :: xmin, xmax

	!********

	integer :: i, n

	area = 0.d0

	n = size(yi)
	if (mod(n, 2) == 1) then
		! Pure 1/3 rule
		do i = 1, n-1, 2
			area = area + 1 * yi(i+0)
			area = area + 4 * yi(i+1)
			area = area + 1 * yi(i+2)
		end do

		area = area * (xmax - xmin) / (n - 1) / 3.d0

	else
		! Use the 1/3 rule for most points and the 3/8 rule for the last 4
		! points
		do i = 1, n-4, 2
			area = area + 1 * yi(i+0)
			area = area + 4 * yi(i+1)
			area = area + 1 * yi(i+2)
		end do
		area = area * (xmax - xmin) / (n - 1) / 3.d0

		area = area &
			+ (1*yi(n-3) + 3*yi(n-2) + 3*yi(n-1) + 1*yi(n-0)) &
			* 3 * (xmax - xmin) / (n - 1) / 8

	end if

end function simpson_integrator_vals

!===============================================================================

! TODO: make trapezoid rule integrators, both fn and val versions

!===============================================================================

double precision function simpson_13_integrator(f, xmin, xmax, dx) result(area)
	! Integrate function `f` from xmin to xmax with step size dx using Simpson's
	! 1/3 rule
	!
	! Error is O(dx**5)

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax, dx

	!********

	integer, parameter :: coefs(*) = [1, 4, 1]

	area = newton_cotes_integrator(f, xmin, xmax, dx, coefs)

end function simpson_13_integrator

!===============================================================================

double precision function simpson_38_integrator(f, xmin, xmax, dx) result(area)
	! Integrate function `f` from xmin to xmax with step size dx using Simpson's
	! 3/8 rule
	!
	! Error is O(dx**5), same order as simpson_13_integrator()

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax, dx

	!********

	integer, parameter :: coefs(*) = [1, 3, 3, 1]

	area = newton_cotes_integrator(f, xmin, xmax, dx, coefs)

end function simpson_38_integrator

!===============================================================================

double precision function milne_integrator(f, xmin, xmax, dx) result(area)
	! Integrate function `f` from xmin to xmax with step size dx using Milne's
	! rule
	!
	! Error is O(dx**7)
	!
	! Note: the text describes another integrator with one more coefficient and
	! error order O(dx**7) not implemented here

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax, dx

	!********

	! Wikipedia calls this one Boole's rule and has something entirely
	! difference for Milne's rule :shrug:
	integer, parameter :: coefs(*) = [7, 32, 12, 32, 7]

	area = newton_cotes_integrator(f, xmin, xmax, dx, coefs)

end function milne_integrator

!===============================================================================

double precision function weddle_integrator(f, xmin, xmax, dx) result(area)
	! Integrate function `f` from xmin to xmax with step size dx using Weddle's
	! rule
	!
	! Error is O(dx**9)

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax, dx

	!********

	integer, parameter :: coefs(*) = [41, 216, 27, 272, 27, 216, 41]

	area = newton_cotes_integrator(f, xmin, xmax, dx, coefs)

end function weddle_integrator

!===============================================================================

double precision function newton_cotes_integrator(f, xmin, xmax, dx, coefs) result(area)
	! Integrate function `f` from xmin to xmax with step size dx using given
	! Newton-Cotes coefficients `coefs`

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax, dx
	integer, intent(in) :: coefs(:)

	!********

	double precision :: dxl, darea, x, f0

	integer :: i, j, j0, n, nc

	nc = size(coefs)
	n = max(1, int((xmax - xmin) / dx))
	dxl = (xmax - xmin) / n

	area = 0.d0
	j0 = 0
	do i = 0, n - 1

		!print *, "i = ", i
		darea = 0.d0
		if (j0 > 0) darea = coefs(nc) * f0
		do j = j0, nc - 1
			x = xmin + i * (xmax - xmin) / n + j * dxl / (nc - 1)
			!print *, "x = ", x
			f0 = f(x)
			darea = darea + coefs(j+1) * f0
		end do
		!print *, "darea = ", darea
		!print *, ""
		area = area + darea

		! Start point of next subinterval is the end point of this one
		j0 = 1

	end do

	area = area * dxl / sum(coefs)

end function newton_cotes_integrator

!===============================================================================

double precision function romberg_integrator_fixed(f, xmin, xmax, n) result(area)
	! Integrate `f` from xmin to xmax using Romberg integration with `n` levels
	! of extrapolation
	!
	! Extrapolation methods estimate the error term of lower-order methods at
	! multiple resolutions

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax
	integer, intent(in) :: n
	!********

	double precision :: dx, s, q
	double precision, allocatable :: t(:)
	integer :: i, k, nt

	allocate(t(n + 1))  ! triangular tableau

	dx = xmax - xmin
	nt = 1  ! number of intervals at level k
	t(1) = 0.5d0 * dx * (f(xmin) + f(xmax))

	do k = 1, n

		s = 0.d0
		dx = 0.5d0 * dx
		nt = 2 * nt
		q = 1
		do i = 1, nt-1, 2
			s = s + f(xmin + i * dx)
		end do
		t(k+1) = 0.5d0 * t(k) + s * dx
		!print *, "tk = ", t(k+1)

		do i = k-1, 0, -1
			q = q * 4
			t(i+1) = t(i+2) + (t(i+2) - t(i+1)) / (q - 1)
			!print *, "ti = ", t(i+1)
		end do
		!print *, "t(1) = ", t(1)  ! this is the best estimate so far after k levels
	end do
	area = t(1)

	!print *, "t = "
	!print "(es22.12)", t

end function romberg_integrator_fixed

!===============================================================================

double precision function romberg_integrator(f, xmin, xmax, tol, max_levels) &
		result(area)
	! Integrate `f` from xmin to xmax using Romberg integration with tolerance
	! `tol`
	!
	! tol could be optional with a default

	use numa__utils
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax
	double precision, intent(in) :: tol
	integer, optional, intent(in) :: max_levels
	!********

	double precision :: dx, s, q, area0
	double precision, allocatable :: t(:)
	integer :: i, k, nt, n
	logical :: converged

	n = 16
	if (present(max_levels)) n = max_levels

	allocate(t(n + 1))  ! triangular tableau

	dx = xmax - xmin
	nt = 1  ! number of intervals at level k
	t(1) = 0.5d0 * dx * (f(xmin) + f(xmax))

	converged = .false.
	do k = 1, n
		area0 = t(1)

		s = 0.d0
		dx = 0.5d0 * dx
		nt = 2 * nt
		q = 1
		do i = 1, nt-1, 2
			s = s + f(xmin + i * dx)
		end do
		t(k+1) = 0.5d0 * t(k) + s * dx
		!print *, "tk = ", t(k+1)

		do i = k-1, 0, -1
			q = q * 4
			t(i+1) = t(i+2) + (t(i+2) - t(i+1)) / (q - 1)
			!print *, "ti = ", t(i+1)
		end do
		!print *, "t(1) = ", t(1)  ! this is the best estimate so far after k levels
		!print *, "diff = ", t(1) - area0

		converged = abs(t(1) - area0) < tol
		if (converged) exit

	end do
	area = t(1)

	if (.not. converged) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": Romberg integrator did not converge under tolerance "// &
			to_str(tol)//" after "//to_str(n)//" iterations"
	end if

	!print *, "t = "
	!print "(es22.12)", t

end function romberg_integrator

!===============================================================================

double precision function gauss2_single(f, xmin, xmax) result(area)
	! Integrate `f` from xmin to xmax using 2-point Gaussian integration with
	! only a single interval

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax
	!********

	double precision, parameter :: &
		x1 = sqrt(3.d0) / 3.d0
	double precision, parameter :: &
		w1 = 1.d0

	area = gauss_general_single(f, xmin, xmax, &
		[ w1, w1], &
		[-x1, x1] &
	)

end function gauss2_single

!===============================================================================

double precision function gauss3_single(f, xmin, xmax) result(area)
	! Integrate `f` from xmin to xmax using 3-point Gaussian integration with
	! only a single interval

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax
	!********

	double precision, parameter :: &
		x1 = 0.7745966692414834d0, &
		x2 = 0.d0
	double precision, parameter :: &
		w1 = 5.d0 / 9, &
		w2 = 8.d0 / 9

	area = gauss_general_single(f, xmin, xmax, &
		[ w1, w1, w2], &
		[-x1, x1, x2] &
	)

end function gauss3_single

!===============================================================================

double precision function gauss4_single(f, xmin, xmax) result(area)
	! Integrate `f` from xmin to xmax using 4-point Gaussian integration with
	! only a single interval

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax
	!********

	double precision, parameter :: &
		x1 = 0.8611363115940526d0, &
		x2 = 0.3399810435848563d0
	double precision, parameter :: &
		w1 = 0.3478548451374538d0, &
		w2 = 0.6521451548625461d0

	area = gauss_general_single(f, xmin, xmax, &
		[ w1, w1,  w2, w2], &
		[-x1, x1, -x2, x2] &
	)

end function gauss4_single

!===============================================================================

double precision function gauss5_single(f, xmin, xmax) result(area)
	! Integrate `f` from xmin to xmax using 5-point Gaussian integration with
	! only a single interval

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax
	!********

	double precision, parameter :: &
		x1 = 0.9061798459386640d0, &
		x2 = 0.5384693101056831d0, &
		x0 = 0.d0
	double precision, parameter :: &
		w1 = 0.2369268850561891d0, &
		w2 = 0.4786286704993665d0, &
		w0 = 0.5688888888888889d0

	area = gauss_general_single(f, xmin, xmax, &
		[ w1, w1, w0,  w2, w2], &
		[-x1, x1, x0, -x2, x2] &
	)

end function gauss5_single

!===============================================================================

double precision function gauss7_single(f, xmin, xmax) result(area)
	! Integrate `f` from xmin to xmax using 7-point Gaussian integration with
	! only a single interval

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax
	!********

	double precision, parameter :: wx(*,*) = reshape([ &
		1.29484966168869693270611432679082018329d-1,  9.49107912342758524526189684047851262401d-1, &
		1.29484966168869693270611432679082018329d-1, -9.49107912342758524526189684047851262401d-1, &
		2.79705391489276667901467771423779582487d-1,  7.41531185599394439863864773280788407074d-1, &
		2.79705391489276667901467771423779582487d-1, -7.41531185599394439863864773280788407074d-1, &
		3.81830050505118944950369775488975133878d-1,  4.05845151377397166906606412076961463347d-1, &
		3.81830050505118944950369775488975133878d-1, -4.05845151377397166906606412076961463347d-1, &
		4.17959183673469387755102040816326530612d-1,  0.00000000000000000000000000000000000000d0   &
	], [2, 7])

	area = gauss_general_single(f, xmin, xmax, &
		wx(1,:), &
		wx(2,:)  &
	)

end function gauss7_single

!===============================================================================

recursive double precision function gk15_aux &
	( &
		f, xmin, xmax, tol, whole, n, &
		neval, is_eps_underflow, is_max_level &
	) &
	result(area)

	! Recursive core.  See gk15_adaptive_integrator() for the user-facing
	! wrapper

	use numa__utils
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax, tol, whole
	integer, intent(in) :: n
	integer, intent(inout) :: neval
	logical, intent(inout) :: is_eps_underflow, is_max_level
	!********

	double precision :: fx, areag, areak, mid, dx2, diff

	integer :: i

	dx2 = 0.5d0 * (xmax - xmin)
	mid = 0.5d0 * (xmin + xmax)

	if (tol/2 == tol .or. xmin == mid) then
		is_eps_underflow = .true.
		area = whole
		return
	end if

	areag = 0.d0
	areak = 0.d0
	do i = 1, NG_GK15
		fx = f(mid + WXK_GK15(2, i) * dx2)
		areag = areag + WXG_GK15(1, i) * fx
		areak = areak + WXK_GK15(1, i) * fx
	end do
	do i = NG_GK15+1, NK_GK15
		fx = f(mid + WXK_GK15(2, i) * dx2)
		areak = areak + WXK_GK15(1, i) * fx
	end do
	neval = neval + NK_GK15
	areag = areag * dx2
	areak = areak * dx2

	!print *, "areag = ", areag
	!print *, "areak = ", areak

	diff = abs(areak - areag)
	if (n <= 0 .and. diff >= tol) is_max_level = .true.
	if (diff < tol .or. n <= 0) then
		area = areak
		return
	end if

	! Here, quadpack does some crazy stuff with pushing interval errors onto a
	! minheap to prioritize which subinterval to evaluate next.  Maybe there's a
	! good reason for this that I don't understand, or maybe it's a limitation
	! to do as best as possible with a limitted number of intervals and fn
	! evals, statically-allocated work arrays, and an avoidance of recursion
	! (except in quadpack's adaptive Simpson and Lobatto integrators)
	!
	! But this method is so simple.  I don't see the need to prioritize
	! subdivision ordering.  If the diff/error is over tolerance, you're going
	! to evaluate both halves eventually anyway unless you give up on staying
	! under tol

	area = &
		gk15_aux(f, xmin, mid, tol/2, areak/2, n-1, neval, is_eps_underflow, is_max_level) + &
		gk15_aux(f, mid, xmax, tol/2, areak/2, n-1, neval, is_eps_underflow, is_max_level)

end function gk15_aux

!===============================================================================

recursive double precision function gk15_adaptive_integrator &
	( &
		f, xmin, xmax, tol, max_levels, iostat &
	) &
	result(area)

	! Integrate `f` from xmin to xmax using the adaptive Gauss-Kronrod method,
	! with 7-point Gauss and 15-point Kronrod

	use numa__utils
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax, tol
	integer, optional, intent(in) :: max_levels
	integer, optional, intent(out) :: iostat
	!********

	character(len = :), allocatable :: msg
	integer :: n, neval
	logical :: is_eps_underflow, is_max_level

	if (present(iostat)) iostat = 0
	if (xmax <= xmin) then
		msg = "xmax is less than xmin in gk15_adaptive_integrator()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if
	if (tol <= 0) then
		msg = "tolerance is 0 in gk15_adaptive_integrator()"
		call PANIC(msg, present(iostat))
		iostat = 2
		return
	end if

	n = 16
	if (present(max_levels)) n = max_levels

	neval = 0
	is_eps_underflow = .false.
	is_max_level = .false.
	area = gk15_aux(f, xmin, xmax, tol, 0.d0, n, neval, is_eps_underflow, is_max_level)

	!print *, "neval gk15 = ", neval

	if (is_eps_underflow) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": gk15_adaptive_integrator() reached numeric tolerance" &
			//" or bound underflow"
	end if
	if (is_max_level) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": gk15_adaptive_integrator() reached max recursion level"
	end if

end function gk15_adaptive_integrator

!===============================================================================

recursive double precision function gk15i_aux &
	( &
		f, xmin, xmax, tol, whole, n, &
		neval, is_eps_underflow, is_max_level &
	) &
	result(area)

	! Recursive core.  See gk15i_adaptive_integrator() for the user-facing
	! wrapper

	use numa__utils
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax, tol, whole
	integer, intent(in) :: n
	integer, intent(inout) :: neval
	logical, intent(inout) :: is_eps_underflow, is_max_level
	!********

	double precision :: x, fx, areag, areak, mid, dx2, diff

	integer :: i

	dx2 = 0.5d0 * (xmax - xmin)
	mid = 0.5d0 * (xmin + xmax)

	if (tol/2 == tol .or. xmin == mid) then
		is_eps_underflow = .true.
		area = whole
		return
	end if

	areag = 0.d0
	areak = 0.d0
	do i = 1, NG_GK15
		x = mid + WXK_GK15(2, i) * dx2
		fx = f(1.d0/x) / x ** 2
		areag = areag + WXG_GK15(1, i) * fx
		areak = areak + WXK_GK15(1, i) * fx
	end do
	do i = NG_GK15+1, NK_GK15
		x = mid + WXK_GK15(2, i) * dx2
		fx = f(1.d0/x) / x ** 2
		areak = areak + WXK_GK15(1, i) * fx
	end do
	neval = neval + NK_GK15
	areag = areag * dx2
	areak = areak * dx2

	!print *, "areag = ", areag
	!print *, "areak = ", areak

	diff = abs(areak - areag)
	if (n <= 0 .and. diff >= tol) is_max_level = .true.
	if (diff < tol .or. n <= 0) then
		area = areak
		return
	end if

	area = &
		gk15i_aux(f, xmin, mid, tol/2, areak/2, n-1, neval, is_eps_underflow, is_max_level) + &
		gk15i_aux(f, mid, xmax, tol/2, areak/2, n-1, neval, is_eps_underflow, is_max_level)

end function gk15i_aux

!===============================================================================

recursive double precision function gk15i_adaptive_integrator &
	( &
		f, xmin, tol, max_levels, iostat &
	) &
	result(area)

	! Integrate `f` from xmin to infinity using the adaptive Gauss-Kronrod method,
	! with 7-point Gauss and 15-point Kronrod
	!
	! This makes use of the following change of variables:
	!
	!     Integrate f(x) dx from xmin to infinity
	!     ==
	!     Integrate 1/t**2 * f(1/t) dt from 0 to 1/xmin
	!
	! See page 184 of Stoer and Bulirsch
	!
	! These are technically "improper" integrals, although improper integrals
	! don't necessarily have infinit bounds, they can also have finite bounds
	! with singularities

	use numa__utils
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, tol
	integer, optional, intent(in) :: max_levels
	integer, optional, intent(out) :: iostat
	!********

	character(len = :), allocatable :: msg
	double precision, parameter :: split_ = 1.d0
	integer :: n, neval
	logical :: is_eps_underflow, is_max_level

	if (present(iostat)) iostat = 0
	if (tol <= 0) then
		msg = "tolerance is 0 in gk15i_adaptive_integrator()"
		call PANIC(msg, present(iostat))
		! iostat starts at 2 here for consistency with other integrators
		iostat = 2
		return
	end if

	n = 16
	if (present(max_levels)) n = max_levels

	neval = 0
	is_eps_underflow = .false.
	is_max_level = .false.

	!print *, "starting gk15i_adaptive_integrator()"
	!print *, "xmin = ", xmin

	if (xmin <= 0.d0) then

		! Integrate from xmin to 1 and 1 to infinity.  The choice of where to
		! split the interval is arbitrary, it just can't be 0
		area = &
			gk15_aux (f, xmin, split_, tol/2, 0.d0, n, neval, is_eps_underflow, is_max_level) + &
			gk15i_aux(f, 0.d0, 1.d0/split_, tol/2, 0.d0, n, neval, is_eps_underflow, is_max_level)

	else
		area = gk15i_aux(f, 0.d0, 1.d0/xmin, tol, 0.d0, n, neval, is_eps_underflow, is_max_level)
	end if

	!print *, "neval gk15i = ", neval

	if (is_eps_underflow) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": gk15i_adaptive_integrator() reached numeric tolerance" &
			//" or bound underflow"
	end if
	if (is_max_level) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": gk15i_adaptive_integrator() reached max recursion level"
	end if

end function gk15i_adaptive_integrator

!===============================================================================

recursive double precision function gk15ii_adaptive_integrator &
	( &
		f, tol, max_levels, iostat &
	) &
	result(area)

	! Integrate `f` from -infinity to infinity using the adaptive Gauss-Kronrod method,
	! with 7-point Gauss and 15-point Kronrod

	use numa__utils
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: tol
	integer, optional, intent(in) :: max_levels
	integer, optional, intent(out) :: iostat
	!********

	character(len = :), allocatable :: msg
	double precision, parameter :: left = -1.d0, right = 1.d0
	integer :: n, neval
	logical :: is_eps_underflow, is_max_level

	if (present(iostat)) iostat = 0
	if (tol <= 0) then
		msg = "tolerance is 0 in gk15ii_adaptive_integrator()"
		call PANIC(msg, present(iostat))
		iostat = 2
		return
	end if

	n = 16
	if (present(max_levels)) n = max_levels

	neval = 0
	is_eps_underflow = .false.
	is_max_level = .false.

	!print *, "starting gk15ii_adaptive_integrator()"
	!print *, "xmin = ", xmin

	area = &
		gk15i_aux(f, 1.d0/left, 0.d0, tol/2, 0.d0, n, neval, is_eps_underflow, is_max_level) + &
		gk15_aux (f, left, right, tol/2, 0.d0, n, neval, is_eps_underflow, is_max_level) + &
		gk15i_aux(f, 0.d0, 1.d0/right, tol/2, 0.d0, n, neval, is_eps_underflow, is_max_level)

	!print *, "neval gk15ii = ", neval

	if (is_eps_underflow) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": gk15ii_adaptive_integrator() reached numeric tolerance" &
			//" or bound underflow"
	end if
	if (is_max_level) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": gk15ii_adaptive_integrator() reached max recursion level"
	end if

end function gk15ii_adaptive_integrator

!===============================================================================

recursive double precision function gk15ni_adaptive_integrator &
	( &
		f, xmax, tol, max_levels, iostat &
	) &
	result(area)

	! Integrate `f` from -infinity to xmax using the adaptive Gauss-Kronrod method,
	! with 7-point Gauss and 15-point Kronrod

	use numa__utils
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmax, tol
	integer, optional, intent(in) :: max_levels
	integer, optional, intent(out) :: iostat
	!********

	character(len = :), allocatable :: msg
	double precision, parameter :: split_ = 1.d0
	integer :: n, neval
	logical :: is_eps_underflow, is_max_level

	if (present(iostat)) iostat = 0
	if (tol <= 0) then
		msg = "tolerance is 0 in gk15ni_adaptive_integrator()"
		call PANIC(msg, present(iostat))
		iostat = 2
		return
	end if

	n = 16
	if (present(max_levels)) n = max_levels

	neval = 0
	is_eps_underflow = .false.
	is_max_level = .false.

	!print *, "starting gk15ni_adaptive_integrator()"
	!print *, "xmax = ", xmax

	if (xmax >= 0.d0) then

		area = &
			gk15_aux (f, -split_, xmax, tol/2, 0.d0, n, neval, is_eps_underflow, is_max_level) + &
			gk15i_aux(f, -1.d0/split_, 0.d0, tol/2, 0.d0, n, neval, is_eps_underflow, is_max_level)

	else
	   area = gk15i_aux(f, 1.d0/xmax, 0.d0, tol, 0.d0, n, neval, is_eps_underflow, is_max_level)
	end if

	!print *, "neval gk15ni = ", neval

	if (is_eps_underflow) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": gk15ni_adaptive_integrator() reached numeric tolerance" &
			//" or bound underflow"
	end if
	if (is_max_level) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": gk15ni_adaptive_integrator() reached max recursion level"
	end if

end function gk15ni_adaptive_integrator

!===============================================================================

double precision function kronrod15_single(f, xmin, xmax) result(area)
	! Integrate `f` from xmin to xmax using 15-point Kronrod integration with
	! only a single interval

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax
	!********

	! Tables of abscissas and weights are taken from Jacob Williams' quadpack
	!
	! Note that the first 7 abscissas of Kronrod-15 are the same as Gaussian-7,
	! although the weights are different.  That's the point of Kronrod
	! integration -- it allows improving a Gaussian integration without
	! repeating fn evals
	double precision, parameter :: wx(*,*) = reshape([ &
		6.30920926299785532907006631892042866651d-2,  9.49107912342758524526189684047851262401d-1, &
		6.30920926299785532907006631892042866651d-2, -9.49107912342758524526189684047851262401d-1, &
		1.40653259715525918745189590510237920400d-1,  7.41531185599394439863864773280788407074d-1, &
		1.40653259715525918745189590510237920400d-1, -7.41531185599394439863864773280788407074d-1, &
		1.90350578064785409913256402421013682826d-1,  4.05845151377397166906606412076961463347d-1, &
		1.90350578064785409913256402421013682826d-1, -4.05845151377397166906606412076961463347d-1, &
		2.09482141084727828012999174891714263698d-1,  0.00000000000000000000000000000000000000d0 , &
		2.29353220105292249637320080589695919936d-2,  9.91455371120812639206854697526328516642d-1, &
		2.29353220105292249637320080589695919936d-2, -9.91455371120812639206854697526328516642d-1, &
		1.04790010322250183839876322541518017444d-1,  8.64864423359769072789712788640926201211d-1, &
		1.04790010322250183839876322541518017444d-1, -8.64864423359769072789712788640926201211d-1, &
		1.69004726639267902826583426598550284106d-1,  5.86087235467691130294144838258729598437d-1, &
		1.69004726639267902826583426598550284106d-1, -5.86087235467691130294144838258729598437d-1, &
		2.04432940075298892414161999234649084717d-1,  2.07784955007898467600689403773244913480d-1, &
		2.04432940075298892414161999234649084717d-1, -2.07784955007898467600689403773244913480d-1  &
	], [2, 15])

	area = gauss_general_single(f, xmin, xmax, &
		wx(1,:), &
		wx(2,:)  &
	)

end function kronrod15_single

!===============================================================================

double precision function gauss_general_single(f, xmin, xmax, w, x) result(area)
	! General Gaussian quadrature with given weights `w` and abscissas `x`

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax
	double precision, intent(in) :: w(:), x(:)
	!********

	double precision :: half_h, mid
	integer :: i

	half_h = 0.5d0 * (xmax - xmin)
	mid    = 0.5d0 * (xmin + xmax)

	area = 0.d0
	do i = 1, size(w)
		area = area + w(i) * f(mid + x(i) * half_h)
	end do
	area = area * half_h

end function gauss_general_single

!===============================================================================

recursive double precision function simpson_adapt_aux &
	( &
		f, a, b, eps, whole, fa, fb, fm, rec, &
		neval, is_eps_underflow, is_max_level &
	) &
	result(area)

	! This private fn is not for end users.  See the non-recursive wrapper
	! simpson_adaptive_integrator() below instead
	!
	! Source:  https://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: a, b, eps, whole, fa, fb, fm
	integer, intent(in) :: rec

	integer, intent(inout) :: neval
	logical, intent(inout) :: is_eps_underflow, is_max_level
	!********

	double precision :: h, m, lm, rm, flm, frm, left, right, delta

	m  = 0.5d0 * (a + b)
	h  = 0.5d0 * (b - a)
	lm = 0.5d0 * (a + m)
	rm = 0.5d0 * (m + b)

	if (eps/2 == eps .or. a == lm) then
		is_eps_underflow = .true.
		area = whole
		return
	end if

	flm = f(lm)
	frm = f(rm)
	neval = neval + 2
	left  = h/6 * (fa + 4*flm + fm)
	right = h/6 * (fm + 4*frm + fb)
	delta = left + right - whole

	if (rec <= 0 .and. abs(delta) > 15 * eps) is_max_level = .true.

	if (rec <= 0 .or. abs(delta) <= 15 * eps) then
		area = left + right + delta / 15.d0
		return
	end if

	area = &
		simpson_adapt_aux( &
		f, a, m, eps/2, left , fa, fm, flm, rec-1, neval, is_eps_underflow, is_max_level) + &
		simpson_adapt_aux( &
		f, m, b, eps/2, right, fm, fb, frm, rec-1, neval, is_eps_underflow, is_max_level)

end function simpson_adapt_aux

!===============================================================================

double precision function simpson_adaptive_integrator &
	( &
		f, xmin, xmax, tol, max_levels, iostat &
	) &
	result(area)

	! Integrate `f` from xmin to xmax using the adaptive Simpson method

	use numa__utils
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax, tol
	integer, optional, intent(in) :: max_levels
	integer, optional, intent(out) :: iostat
	!********

	character(len = :), allocatable :: msg
	double precision :: h, s, fa, fb, fm
	integer :: n, neval
	logical :: is_eps_underflow, is_max_level

	if (present(iostat)) iostat = 0

	area = 0.d0
	h = xmax - xmin
	if (h <= 0) then
		msg = "xmax is less than xmin in simpson_adaptive_integrator()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	n = 16
	if (present(max_levels)) n = max_levels

	if (tol <= 0) then
		msg = "tolerance is 0 in simpson_adaptive_integrator()"
		call PANIC(msg, present(iostat))
		iostat = 2
		return
	end if

	! Bootstrap the first level then call the recursive core
	fa = f(xmin)
	fm = f(0.5d0 * (xmin + xmax))
	fb = f(xmax)

	neval = 3  ! number of fn evaluations
	is_eps_underflow = .false.
	is_max_level     = .false.

	s = h / 6 * (fa * 4*fm + fb)

	area = simpson_adapt_aux(f, xmin, xmax, tol, s, fa, fb, fm, n, &
		neval, is_eps_underflow, is_max_level)

	!! Make an optional out arg?
	!print *, "neval simp = ", neval

	if (is_eps_underflow) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": simpson_adaptive_integrator() reached numeric tolerance" &
			//" or bound underflow"
	end if
	if (is_max_level) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": simpson_adaptive_integrator() reached max recursion level"
	end if

end function simpson_adaptive_integrator

!===============================================================================

end module numa__integrate

