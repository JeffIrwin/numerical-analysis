
module numerical_analysis_m

	implicit none

	double precision, parameter :: PI = 4 * atan(1.d0)
	double complex, parameter :: IMAG = (0.d0, 1.d0)  ! sqrt(-1)

	abstract interface
		! This is a function interface for passing callbacks
		function fn_f64_to_f64(x) result(fx)
			double precision, intent(in) :: x
			double precision :: fx
		end function
	end interface

contains

!===============================================================================
double precision function lagrange_interpolater(xi, f, x) result(fx)

	! Given support points `xi` and function `f`, interpolate f at `x` using
	! Lagrange interpolation
	!
	! "Lagrange's formula is, in general, not as suitable for actual
	! calculations as some other methods to be described below ... ."

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

	fx = lagrange_interpolater_vals(xi, fi, x)

end function lagrange_interpolater

!===============================================================================
double precision function neville_interpolater(xi, f, x) result(fx)

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

	fx = neville_interpolater_vals(xi, fi, x)

end function neville_interpolater

!===============================================================================
double precision function lagrange_interpolater_vals(xi, fi, x) result(fx)

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

end function lagrange_interpolater_vals

!===============================================================================
double precision function neville_interpolater_vals(xi, fi, x) result(fx)

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

end function neville_interpolater_vals

!===============================================================================
function newton_interpolater(xi, f, xs) result(fxs)

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

	fxs = newton_interpolater_vals(xi, fi, xs)

end function newton_interpolater

!===============================================================================
function newton_interpolater_vals(xi, fi, xs) result(fxs)

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

	! Interpolate f, O(n1)
	allocate(fxs( size(xs) ))
	fxs = a(n1)
	do i = n1-1, 1, -1
		fxs = fxs * (xs - xi(i)) + a(i)
	end do
	!print *, "fxs = ", fxs

end function newton_interpolater_vals

!===============================================================================
double precision function neville_rational_interpolater(xi, f, x) result(fx)

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

	fx = neville_rational_interpolater_vals(xi, fi, x)

end function neville_rational_interpolater

!===============================================================================
double precision function neville_rational_interpolater_vals(xi, fi, x) result(fx)

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

end function neville_rational_interpolater_vals

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

recursive function fft(x) result(xx)

	! Fast Fourier transform, iterative (non-recursive) implementation

	! TODO: pad the input to a size of a power of 2

	! Page 86

	double complex, intent(in) :: x(:)
	double complex, allocatable :: xx(:)

	!********

	double complex :: e, u, v

	integer :: m, j, r, nj, n, n2

	n2 = size(x)
	n = 0
	do while (2 ** n < n2)
		n = n + 1
	end do
	!print *, "n, n2 = ", n, n2

	xx = x
	!xx = x / n2

	! Bit reversal of indices
	do j = 0, n2 - 1
		nj = bit_reverse(j, n)
		!print *, "j, nj = ", j, nj
		if (j < nj) xx([j+1, nj+1]) = xx([nj+1, j+1])
	end do

	do m = 1, n
		do j = 0, 2 ** (m-1) - 1
			!e = e_m^j
			!w = exp(-2 * IMAG * PI * [(i, i = 0, n-1)] / n)
			e = exp((-2.d0 * IMAG * PI * j) / (2.d0 ** m))
			do r = 0, n2 - 1, 2 ** m  ! TODO: shift for 1-indexing
				u = xx(r + j + 1)
				v = xx(r + j + 2 ** (m-1) + 1) * e
				!print *, "u, v = ", u, v

				xx(r + j + 1) = u + v
				xx(r + j + 2 ** (m-1) + 1) = u - v

			end do
		end do
	end do

end function fft

!===============================================================================

recursive function ifft(x) result(xx)

	! Inverse fast Fourier transform, iterative (non-recursive) implementation

	double complex, intent(in) :: x(:)
	double complex, allocatable :: xx(:)

	!********

	double complex :: e, u, v

	integer :: m, j, r, nj, n, n2

	n2 = size(x)
	n = 0
	do while (2 ** n < n2)
		n = n + 1
	end do
	!print *, "n, n2 = ", n, n2

	! TODO: I think division should happen at the end, not here.  Test harness
	! will need to change
	xx = x
	!xx = x / n2

	! Bit reversal of indices
	do j = 0, n2 - 1
		nj = bit_reverse(j, n)
		!print *, "j, nj = ", j, nj
		if (j < nj) xx([j+1, nj+1]) = xx([nj+1, j+1])
	end do

	do m = 1, n
		do j = 0, 2 ** (m-1) - 1

			! The only difference between fft and ifft is the sign of this
			! exponent
			e = exp((2.d0 * IMAG * PI * j) / (2.d0 ** m))

			do r = 0, n2 - 1, 2 ** m  ! TODO: shift for 1-indexing
				u = xx(r + j + 1)
				v = xx(r + j + 2 ** (m-1) + 1) * e
				!print *, "u, v = ", u, v

				xx(r + j + 1) = u + v
				xx(r + j + 2 ** (m-1) + 1) = u - v

			end do
		end do
	end do

end function ifft

!===============================================================================

recursive function fft_rec(x) result(xx)

	! Fast Fourier transform, recursive implementation

	! TODO: add a wrapper that pads the input to a size of a power of 2

	double complex, intent(in) :: x(:)
	double complex, allocatable :: xx(:)

	!********

	double complex, allocatable :: xx_even(:), xx_odd(:), w(:)

	integer :: i, n

	n = size(x)
	!print *, "n = ", n

	if (n == 1) then
		xx = x
		return
	end if

	xx_even = fft_rec(x(1: n: 2))
	xx_odd  = fft_rec(x(2: n: 2))
	w = exp(-2 * IMAG * PI * [(i, i = 0, n-1)] / n)

	xx = &
	[    &
		xx_even + w(1    : n/2) * xx_odd, &
		xx_even + w(n/2+1: n  ) * xx_odd  &
	]

end function fft_rec

!===============================================================================

recursive function ifft_rec(x) result(xx)

	! Inverse fast Fourier transform

	! TODO: add a wrapper that pads the input to a size of a power of 2

	double complex, intent(in) :: x(:)
	double complex, allocatable :: xx(:)

	!********

	double complex, allocatable :: xx_even(:), xx_odd(:), w(:)

	integer :: i, n

	n = size(x)
	!print *, "n = ", n

	if (n == 1) then
		xx = x
		return
	end if

	xx_even = ifft_rec(x(1: n: 2))
	xx_odd  = ifft_rec(x(2: n: 2))
	w = exp(2 * IMAG * PI * [(i, i = 0, n-1)] / n)

	xx = &
	[    &
		xx_even + w(1    : n/2) * xx_odd, &
		xx_even + w(n/2+1: n  ) * xx_odd  &
	]

end function ifft_rec

!===============================================================================

! TODO: maybe these example fns (used as callbacks) should be moved to another
! file

double precision function log_fn(x) result(log_x)

	! The built-in Fortran fn `log` cannot be passed directly as a function
	! callback:
	!
	!      145 |         fx = lagrange_interpolater(xi, log, 11.1d0)
	!          |                                          1
	!    Error: Symbol ‘log’ at (1) has no IMPLICIT type
	!
	! So we define this wrapper

	double precision, intent(in) :: x

	log_x = log(x)

end function log_fn

!********

double precision function exp_fn(x) result(exp_x)
	double precision, intent(in) :: x
	exp_x = exp(x)
end function exp_fn

!********

double precision function cotd_fn(x)
	double precision, intent(in) :: x
	cotd_fn = 1 / tand(x)
end function cotd_fn

!===============================================================================

end module numerical_analysis_m

