
#include "panic.F90"

!> Module for interpolation with polynomials, rational polynomials, fast Fourier
!> transforms, and splines
module numa__interp

	use numa__core
	implicit none

	! Enum for spline_general() switch/case
	integer, parameter :: &
		SPLINE_CASE_NO_CURVE   = 1, &
		SPLINE_CASE_PERIODIC   = 2, &
		SPLINE_CASE_PRESCRIBED = 3

	private :: spline_general

contains

!===============================================================================

!> @brief  Interpolates a function at given points using Lagrange interpolation
!>
!> @details  "Lagrange's formula is, in general, not as suitable for actual
!> calculations as some other methods to be described below ... ."
!>
!> @param[in]  xi  Support point x values where the function is evaluated directly
!> @param[in]  f   Function to be evaluated and interpolated
!> @param[in]  x   Point x values where the function will be interpolated
!> @return         Interpolated y value results
double precision function lagrange_interpolator(xi, f, x) result(fx)
	double precision, intent(in) :: xi(:)
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: x

	!********

	double precision, allocatable :: fi(:)

	integer :: i, n1

	! Degree of polynomial + 1
	n1 = size(xi)

	! Function values `fi` at those support points
	allocate(fi(n1))
	do i = 1, n1
		fi(i) = f(xi(i))
	end do
	!print *, "xi = ", xi
	!print *, "fi = ", fi

	fx = lagrange_interpolator_vals(xi, fi, x)

end function lagrange_interpolator

!===============================================================================
double precision function neville_interpolator(xi, f, x) result(fx)

	double precision, intent(in) :: xi(:)
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: x

	!********

	double precision, allocatable :: fi(:)

	integer :: i, n1

	! Degree of polynomial + 1
	n1 = size(xi)

	! Function values `fi` at those support points
	allocate(fi(n1))
	do i = 1, n1
		fi(i) = f(xi(i))
	end do

	fx = neville_interpolator_vals(xi, fi, x)

end function neville_interpolator

!===============================================================================
double precision function lagrange_interpolator_vals(xi, fi, x) result(fx)

	! Given support points `xi` and function values `fi` at those points,
	! interpolate f at `x` using Lagrange interpolation
	!
	! This version takes the function's values instead of the function itself,
	! which may be useful if you only have the values from a black box

	double precision, intent(in) :: xi(:), fi(:)
	double precision, intent(in) :: x

	!********

	double precision :: prod

	integer :: i, k, n1

	! Degree of polynomial + 1
	n1 = size(xi)

	! Interpolate at f(x) at `x`

	! Lagrange interpolation, eq 2.1.1.4, page 39
	fx = 0.d0
	do i = 1, n1
		prod = 1.d0
		do k = 1, n1
			if (i == k) cycle
			prod = prod * (x - xi(k)) / (xi(i) - xi(k))
		end do
		fx = fx + fi(i) * prod
	end do

end function lagrange_interpolator_vals

!===============================================================================
double precision function neville_interpolator_vals(xi, fi, x) result(fx)

	! Lagrange and Neville are both O(n1**2) time, while Neville has half the
	! constant factor at the cost of an extra O(n1) space for t(:)

	double precision, intent(in) :: xi(:), fi(:)
	double precision, intent(in) :: x

	!********

	double precision, allocatable :: t(:)

	integer :: i, j, n1

	! Degree of polynomial + 1
	n1 = size(xi)

	! Interpolate at f(x) at `x`

	! Neville's algorithm, eq 2.1.2.5, page 42
	allocate(t(n1))
	do i = 1, n1
		t(i) = fi(i)
		do j = i-1, 1, -1
			t(j) = t(j+1) + (t(j+1) - t(j)) * (x - xi(i)) / (xi(i) - xi(j))
		end do
	end do
	fx = t(1)

end function neville_interpolator_vals

!===============================================================================
function newton_interpolator(xi, f, xs) result(fxs)

	double precision, intent(in) :: xi(:)
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xs(:)

	double precision, allocatable :: fxs(:)

	!********

	double precision, allocatable :: fi(:)

	integer :: i, n1

	! Degree of polynomial + 1
	n1 = size(xi)

	! Function values `fi` at those support points
	allocate(fi(n1))
	do i = 1, n1
		fi(i) = f(xi(i))
	end do

	fxs = newton_interpolator_vals(xi, fi, xs)

end function newton_interpolator

!===============================================================================
function newton_interpolator_vals(xi, fi, xs) result(fxs)

	! Newton's interpolation formula:  divided differences
	!
	! Unlike Lagrange and Neville, this takes an array of points xs(:) to be
	! interpolated
	!
	! Alternatively, this could be split into init and eval stages. The whole
	! point of Newton's algorithm is that it's more efficient than Lagrange and
	! Neville for many interpolation points, but only if you don't repeat the
	! O(n1**2) part for every point

	double precision, intent(in) :: xi(:), fi(:)
	double precision, intent(in) :: xs(:)

	double precision, allocatable :: fxs(:)

	!********

	double precision, allocatable :: a(:), t(:)

	integer :: i, j, n1

	! Degree of polynomial + 1
	n1 = size(xi)

	! Find the polynomial coefficients a(:), O(n1**2).  Page 46
	allocate(a(n1))
	allocate(t(n1))
	do i = 1, n1
		t(i) = fi(i)
		do j = i-1, 1, -1
			t(j) = (t(j+1) - t(j)) / (xi(i) - xi(j))
		end do
		a(i) = t(1)
	end do

	! Interpolate f, O(n1) per point
	allocate(fxs( size(xs) ))
	fxs = a(n1)
	do i = n1-1, 1, -1
		fxs = fxs * (xs - xi(i)) + a(i)
	end do
	!print *, "fxs = ", fxs

end function newton_interpolator_vals

!===============================================================================
double precision function neville_rational_interpolator(xi, f, x) result(fx)

	double precision, intent(in) :: xi(:)
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: x

	!********

	double precision, allocatable :: fi(:)

	integer :: i, n1

	! Degree of polynomial + 1
	n1 = size(xi)

	! Function values `fi` at those support points
	allocate(fi(n1))
	do i = 1, n1
		fi(i) = f(xi(i))
	end do

	fx = neville_rational_interpolator_vals(xi, fi, x)

end function neville_rational_interpolator

!===============================================================================
double precision function neville_rational_interpolator_vals(xi, fi, x) result(fx)

	! Lagrange and Neville are both O(n1**2) time, while Neville has half the
	! constant factor at the cost of an extra O(n1) space for t(:)

	double precision, intent(in) :: xi(:), fi(:)
	double precision, intent(in) :: x

	!********

	double precision :: tk2
	double precision, allocatable :: t(:), t0(:)

	integer :: i, k, n1

	! Degree of polynomial + 1
	n1 = size(xi)

	! Interpolate at f(x) at `x`

	! Neville's algorithm with rational polynomials, eq 2.2.3.8, page 72

	allocate(t(n1), t0(n1))
	t0 = 0.d0
	do i = 1, n1
		t(1) = fi(i)

		tk2 = 0.d0
		do k = 2, i
			if (k > 2) tk2 = t0(k-2)

			! "Note that this recursion formula differs from the corresponding
			! polynomial formula (2.1.2.5) only by the expression in brackets*"
			t(k) = t(k-1) + (t(k-1) - t0(k-1)) &
				/ ((x - xi(i-k+1)) / (x - xi(i)) * &
				( &
					1.d0 - (t(k-1) - t0(k-1)) / (t(k-1) - tk2) &  ! "expression in brackets"
				) - 1.d0)

		end do
		t0 = t
	end do
	fx = t(n1)

end function neville_rational_interpolator_vals

!===============================================================================
integer function bit_reverse(j, n) result(k)
	! Reverse the lowest n bits of number j

	integer, intent(in) :: j, n
	!********
	integer :: i

	k = 0
	do i = 0, n - 1
		!if (iand(ishft(j, i), 1) == 1) k = ior(k, ishft(1, -i))
		if (btest(j, i)) k = ibset(k, n - i - 1)
	end do

end function bit_reverse

!===============================================================================

function fft(x) result(xx)
	! Forward fast Fourier transform
	double complex, intent(in) :: x(:)
	double complex, allocatable :: xx(:)
	xx = fft_core(x, invert_ = .false.)
end function fft

!********

function ifft(x) result(xx)
	! Inverse fast Fourier transform
	double complex, intent(in) :: x(:)
	double complex, allocatable :: xx(:)
	xx = fft_core(x, invert_ = .true.)
end function ifft

!===============================================================================

function fft_core(x, invert_) result(xx)
	! Fast Fourier transform, iterative (non-recursive) implementation, using
	! the Cooley-Tukey method

	! Page 86

	double complex, intent(in) :: x(:)
	logical, intent(in) :: invert_
	double complex, allocatable :: xx(:)

	!********

	double precision :: sgn
	double complex :: e, u, v

	integer :: m, j, r, nj, n, n2

	! Pad with zeros to size n2
	n = 0
	do while (2 ** n < size(x))
		n = n + 1
	end do
	n2 = 2 ** n
	allocate(xx(n2))
	xx = 0.d0
	xx(1: size(x)) = x
	!print *, "n, n2 = ", n, n2

	if (invert_) then
		! De-normalize
		xx = xx * n2
	end if

	! Permute array by bit reversal of indices
	do j = 0, n2 - 1
		nj = bit_reverse(j, n)
		!print *, "j, nj = ", j, nj
		if (j < nj) xx([j+1, nj+1]) = xx([nj+1, j+1])
	end do

	sgn = -2.d0
	if (invert_) sgn = 2.d0

	do m = 1, n
		do j = 0, 2 ** (m-1) - 1

			! The main difference between fft and ifft is the sign of this
			! exponent
			e = exp(sgn * IMAG_ * PI * j / 2 ** m)

			do r = 1, n2, 2 ** m
				u = xx(r + j)
				v = xx(r + j + 2 ** (m-1)) * e

				xx(r + j) = u + v
				xx(r + j + 2 ** (m-1)) = u - v

			end do
		end do
	end do

	! Normalize amplitudes
	xx = xx / n2
	!xx = xx / size(x)

end function fft_core

!===============================================================================

function spline_no_curve(xi, yi, x, iostat) result(fx)
	! This function finds the cubic spline with the "no curvature" boundary
	! condition, i.e. the 2nd derivative is 0 at both endpoints, or case (a) in
	! the text
	!
	! This is also known as the natural spline

	use numa__utils
	double precision, intent(in)   :: xi(:)
	double precision, intent(in)   :: yi(:)
	double precision, intent(in)   :: x(:)
	integer, optional, intent(out) :: iostat
	double precision, allocatable  :: fx(:)
	!********
	character(len = :), allocatable :: msg
	integer :: io

	if (present(iostat)) iostat = 0
	fx = spline_general(xi, yi, x, SPLINE_CASE_NO_CURVE, 0.d0, 0.d0, io)
	if (io /= 0) then
		msg = "spline_general() failed in spline_no_curve()"
		call PANIC(msg, present(iostat))
		iostat = io
		return
	end if

end function spline_no_curve

!********

function spline_periodic(xi, yi, x, iostat) result(fx)
	! This function finds the cubic spline with periodic boundary conditions,
	! i.e. the 1st and 2nd derivates match at the start and end points, case (b)
	! in the text

	use numa__utils
	double precision, intent(in)   :: xi(:)
	double precision, intent(in)   :: yi(:)
	double precision, intent(in)   :: x(:)
	integer, optional, intent(out) :: iostat
	double precision, allocatable  :: fx(:)
	!********
	character(len = :), allocatable :: msg
	integer :: io

	if (present(iostat)) iostat = 0
	fx = spline_general(xi, yi, x, SPLINE_CASE_PERIODIC, 0.d0, 0.d0, io)
	if (io /= 0) then
		msg = "spline_general() failed in spline_periodic()"
		call PANIC(msg, present(iostat))
		iostat = io
		return
	end if

end function spline_periodic

!********

function spline_prescribed(xi, yi, x, dy_start, dy_end, iostat) result(fx)
	! This function finds the cubic spline with prescribed first derivatives at
	! the end points, i.e. case (c) in the text

	use numa__utils
	double precision, intent(in)   :: xi(:)
	double precision, intent(in)   :: yi(:)
	double precision, intent(in)   :: x(:)
	double precision, intent(in)   :: dy_start, dy_end
	integer, optional, intent(out) :: iostat
	double precision, allocatable  :: fx(:)
	!********
	character(len = :), allocatable :: msg
	integer :: io

	if (present(iostat)) iostat = 0
	fx = spline_general(xi, yi, x, SPLINE_CASE_PRESCRIBED, dy_start, dy_end, io)
	if (io /= 0) then
		msg = "spline_general() failed in spline_prescribed()"
		call PANIC(msg, present(iostat))
		iostat = io
		return
	end if

end function spline_prescribed

!===============================================================================

function spline_general(xi, yi, x, case_, dy_start, dy_end, iostat) result(fx)

	use numa__utils
	use numa__linalg

	double precision, intent(in)  :: xi(:)
	double precision, intent(in)  :: yi(:)
	double precision, intent(in)  :: x(:)

	double precision, intent(in)  :: dy_start, dy_end
	integer, intent(in) :: case_
	integer, intent(out) :: iostat

	double precision, allocatable :: fx(:)
	!********

	character(len = :), allocatable :: msg
	double precision, allocatable :: h(:), lambda(:), mu(:), &
		d(:), a(:,:), aj, bj
	integer :: i, j, n, io

	iostat = 0
	if (.not. any(case_ == &
		[SPLINE_CASE_NO_CURVE, SPLINE_CASE_PERIODIC, SPLINE_CASE_PRESCRIBED])) then
		msg = "bad spline case in spline_general()"
		call PANIC(msg, .true.)
		iostat = 1
		return
	end if

	n = size(xi)

	h = xi(2:) - xi(1: size(xi) - 1)

	!print *, "n = ", n
	!print *, "h = ", h

	allocate(lambda ( n ))
	allocate(mu     ( n ))
	allocate(d      ( n ))

	do j = 2, n-1
		lambda(j) = h(j) / (h(j-1) + h(j))
		mu(j) = 1.d0 - lambda(j)
		d(j) = 6.d0 / (h(j-1) + h(j)) * &
			((yi(j+1) - yi(j)) / h(j) - (yi(j) - yi(j-1)) / h(j-1))
	end do

	if (case_ == SPLINE_CASE_NO_CURVE) then
		lambda(1) = 0.d0
		d(1) = 0.d0
		mu(n) = 0.d0
		d(n) = 0.d0

	else if (case_ == SPLINE_CASE_PERIODIC) then
		lambda(n) = h(1) / (h(n-1) + h(1))
		mu(n) = 1.d0 - lambda(n)
		d(n) = 6.d0 / (h(n-1) + h(1)) * &
			((yi(2) - yi(n)) / h(1) - (yi(n) - yi(n-1)) / h(n-1))

	else if (case_ == SPLINE_CASE_PRESCRIBED) then
		lambda(1) = 1.d0
		d(1) = 6.d0 / h(1) * ((yi(2) - yi(1)) / h(1) - dy_start)
		mu(n) = 1.d0
		d(n) = 6.d0 / h(n-1) * ((dy_end - (yi(n) - yi(n-1)) / h(n-1)))

	end if

	if (case_ == SPLINE_CASE_PERIODIC) then
		! Trim
		d = d(2:)

		! Compose the tridiagonal system matrix `a`
		allocate(a( 3, n-1 ))
		a(1, 1) = mu(2)
		a(2, 1) = 2.d0
		a(3, 1) = lambda(2)
		do j = 2, n-2
			a(1, j) = mu(j+1)
			a(2, j) = 2.d0
			a(3, j) = lambda(j+1)
		end do
		a(1, n-1) = mu(n)
		a(2, n-1) = 2.d0
		a(3, n-1) = lambda(n)

		!print *, "a = "
		!print "(3es18.6)", a

		call tridiag_corner_invmul(a, d, io)
		if (io /= 0) then
			msg = "tridiag_corner_invmul() failed in spline_general()"
			call PANIC(msg, .true.)
			iostat = 2
			return
		end if

		! Untrim
		d = [d(n-1), d]

	else
		! Prescribed or no-curve
		!
		! Compose the tridiagonal system matrix `a`
		allocate(a( 3, n ))
		a(1, 1) = 0.d0
		a(2, 1) = 2.d0
		a(3, 1) = lambda(1)
		do j = 2, n-1
			a(1, j) = mu(j)
			a(2, j) = 2.d0
			a(3, j) = lambda(j)
		end do
		a(1, n) = mu(n)
		a(2, n) = 2.d0
		a(3, n) = 0.d0

		!print *, "a = "
		!print "(3es18.6)", a
		!print *, "lambda = ", lambda
		!print *, "mu     = ", mu
		!print *, "d       = ", d

		call tridiag_invmul(a, d, io)
		if (io /= 0) then
			msg = "tridiag_corner_invmul() failed in spline_general()"
			call PANIC(msg, .true.)
			iostat = 3
			return
		end if

		!print *, "moments = ", d

	end if

	allocate(fx( size(x) ))

	if (x(1) < xi(1)) then
		msg = "`x` underflows first spline control point in spline_general()"
		call PANIC(msg, .true.)
		iostat = 4
		return
	end if
	if (x(size(x)) > xi(size(xi))) then
		msg = "`x` overflows last spline control point in spline_general()"
		call PANIC(msg, .true.)
		iostat = 5
		return
	end if

	! Evaluate the spline
	j = 1
	aj = (yi(j+1) - yi(j)) / h(j) - h(j) / 6.d0 * (d(j+1) - d(j))
	bj = yi(j) - d(j) * (h(j) ** 2) / 6.d0
	do i = 1, size(x)

		if (i > 1) then
			if (x(i) <= x(i-1)) then
				msg = "`x` is not sorted ascending in spline_general()"
				call PANIC(msg, .true.)
				iostat = 6
				return
			end if
		end if

		! Locate the spline segment
		do while (xi(j+1) <= x(i) .and. j < n)
			!print *, "inc j at x = ", x(i)
			j = j + 1
			aj = (yi(j+1) - yi(j)) / h(j) - h(j) / 6.d0 * (d(j+1) - d(j))
			bj = yi(j) - d(j) * (h(j) ** 2) / 6.d0
		end do

		fx(i) = &
			d(j)   * ((xi(j+1) -  x(i)) ** 3) / (6.d0 * h(j)) + &
			d(j+1) * ((x(i)    - xi(j)) ** 3) / (6.d0 * h(j)) + &
			aj * (x(i) - xi(j)) + bj

	end do

end function spline_general

!===============================================================================

function bezier_curve(xc, t) result(x)
	! Interpolate an n-dimensional bezier curve of degree nc-1 for control
	! points `xc` at interpolation parameters `t` in range [0.0, 1.0]

	double precision, intent(in) :: xc(:,:)
	double precision, intent(in) :: t(:)
	double precision, allocatable :: x(:,:)
	!********

	double precision :: t_
	double precision, allocatable :: v(:,:), v0(:,:)

	integer :: i, j, it, nd, nc, nt

	nd = size(xc, 1)  ! number of dimensions
	nc = size(xc, 2)  ! number of control points
	nt = size(t)

	!print *, "nd, nc, nt = ", nd, nc, nt
	!print *, "xc = "
	!print "(2es18.6)", xc

	allocate(x(nd, nt))
	do it = 1, nt
		t_ = t(it)

		! De Casteljau's algorithm
		v = xc
		do i = nc-1, 1, -1
			v0 = v
			do j = 1, i
				! lerp
				v(:, j) = (1.d0 - t_) * v0(:, j) + t_ * v0(:, j+1)
			end do
		end do
		x(:, it) = v(:, 1)

	end do

	!print *, "x = "
	!print "(2es18.6)", x

end function bezier_curve

!===============================================================================

function cardinal_spline(xc, t, tension, iostat) result(x)

	! Interpolate a cubic cardinal spline through control points `xc` with
	! interpolation parameters `t` in range [0, nc-1], where `nc` is the number
	! of control points, with a given tension.  The tension parameter is
	! consistent with microsoft's documentation:
	!
	!     https://learn.microsoft.com/en-us/windows/win32/gdiplus/-gdiplus-cardinal-splines-about
	!
	! That is, as tension approaches 0, you get straight lines.  Of course you
	! can't use exactly 0, but some small number like 0.001 approximates
	! straight lines.  A tension of 1 is the Catmull-Rom spline
	!
	!     | Note that a tension of 0 corresponds to infinite physical tension,
	!     | forcing the curve to take the shortest way (straight lines) between
	!     | points. A tension of 1 corresponds to no physical tension, allowing
	!     | the spline to take the path of least total bend. With tension values
	!     | greater than 1, the curve behaves like a compressed spring, pushed to
	!     | take a longer path.
	!
	! This seems to be the opposite of the tension mentioned on wikipedia:
	!
	!     https://en.wikipedia.org/wiki/Cubic_Hermite_spline

	use numa__utils
	double precision, intent(in) :: xc(:,:)
	double precision, intent(in) :: t(:)
	double precision, intent(in) :: tension
	integer, optional, intent(out) :: iostat

	double precision, allocatable :: x(:,:)
	!********

	character(len = :), allocatable :: msg
	double precision :: t_, tmod
	double precision, allocatable :: xcl(:,:)
	double precision, allocatable :: v(:,:), v0(:,:)
	integer :: i, j, it, nd, nc, nt, ncl, segment

	if (present(iostat)) iostat = 0

	if (tension <= 0) then
		msg = "tension is negative in cardinal_spline()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	nd = size(xc, 1)  ! number of dimensions
	nc = size(xc, 2)  ! number of control points
	nt = size(t)

	!print *, "nd, nc, nt = ", nd, nc, nt
	!print *, "xc = "
	!print "(2es18.6)", xc

	! Construct the full set of control points xcl from the given interpolation
	! points xc
	ncl = 3 * nc
	allocate(xcl(nd, ncl))

	! First and last points are edge cases

	! First point
	xcl(:,1) = xc(:,1) + 2 * tension / 3 * (xc(:,1) - xc(:,2))  ! not used
	xcl(:,2) = xc(:,1)
	xcl(:,3) = xc(:,1) - 2 * tension / 3 * (xc(:,1) - xc(:,2))

	! Interior control points
	do j = 2, nc-1
		!print *, "j = ", j
		xcl(:, 3*j-2) = xc(:,j) + tension / 3 * (xc(:,j-1) - xc(:,j+1))
		xcl(:, 3*j-1) = xc(:,j)
		xcl(:, 3*j-0) = xc(:,j) - tension / 3 * (xc(:,j-1) - xc(:,j+1))
	end do

	! Last point
	xcl(:,ncl-2) = xc(:,nc) + 2 * tension / 3 * (xc(:,nc-1) - xc(:,nc))
	xcl(:,ncl-1) = xc(:,nc)
	xcl(:,ncl-0) = xc(:,nc) - 2 * tension / 3 * (xc(:,nc-1) - xc(:,nc))  ! not used

	!print *, "xcl = "
	!print "(2es18.6)", xcl

	! Using groups of 4 control points, make cubic bezier splines

	allocate(x(nd, nt))
	do it = 1, nt
		t_ = t(it)

		segment = max(1, min(nc-1, floor(t_) + 1))
		tmod = t_ - (segment-1)  ! can be 1 on final segment

		!print *, "t_, tmod, segment = ", t_, tmod, segment

		! De Casteljau's algorithm
		v = xcl(:, segment * 3 - 1: segment * 3 + 2)
		do i = 3, 1, -1
			v0 = v
			do j = 1, i
				! lerp
				v(:, j) = (1.d0 - tmod) * v0(:, j) + tmod * v0(:, j+1)
			end do
		end do
		x(:, it) = v(:, 1)

	end do

end function cardinal_spline

!===============================================================================

end module numa__interp

