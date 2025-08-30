
#include "panic.F90"

!> Module for interpolation with polynomials, rational polynomials, fast Fourier
!> transforms, and splines
module numa__interp

	use numa__core
	implicit none

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

end module numa__interp

