
module numa

	implicit none

	double precision, parameter :: PI = 4 * atan(1.d0)
	double complex, parameter :: IMAG = (0.d0, 1.d0)  ! sqrt(-1)

	integer, parameter :: &
		SPLINE_CASE_NO_CURVE   = 1, &
		SPLINE_CASE_PERIODIC   = 2, &
		SPLINE_CASE_PRESCRIBED = 3

	abstract interface
		! This is a function interface for passing callbacks
		function fn_f64_to_f64(x) result(fx)
			double precision, intent(in) :: x
			double precision :: fx
		end function
	end interface

	! If this file gets too long, it might be good to split it up, e.g. into
	! interpolate.f90, fft.f90, (and in the future integrate.f90), etc.

contains

!===============================================================================
double precision function lagrange_interpolator(xi, f, x) result(fx)

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

	! Interpolate f, O(n1)
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
double precision function norm2c(v)

	double complex, intent(in) :: v(:)

	! It's weird that Fortran's intrinsic norm2() can't handle complex args.
	! Intel MKL has dznrm2() which is equivalent
	!
	! Note that dot_product() already conjugates one of the arguments, so there
	! is no need for additional conjugation

	norm2c = dble(sqrt(dot_product(v, v)))

end function norm2c

!===============================================================================
function fft(x) result(xx)
	! Forward fast Fourier transform
	double complex, intent(in) :: x(:)
	double complex, allocatable :: xx(:)
	xx = fft_core(x, invert = .false.)
end function fft

!********

function ifft(x) result(xx)
	! Inverse fast Fourier transform
	double complex, intent(in) :: x(:)
	double complex, allocatable :: xx(:)
	xx = fft_core(x, invert = .true.)
end function ifft

!===============================================================================

function fft_core(x, invert) result(xx)
	! Fast Fourier transform, iterative (non-recursive) implementation, using
	! the Cooley-Tukey method

	! Page 86

	double complex, intent(in) :: x(:)
	logical, intent(in) :: invert
	double complex, allocatable :: xx(:)

	!********

	double precision :: sign_
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

	if (invert) then
		! De-normalize
		xx = xx * n2
	end if

	! Permute array by bit reversal of indices
	do j = 0, n2 - 1
		nj = bit_reverse(j, n)
		!print *, "j, nj = ", j, nj
		if (j < nj) xx([j+1, nj+1]) = xx([nj+1, j+1])
	end do

	sign_ = -2.d0
	if (invert) sign_ = 2.d0

	do m = 1, n
		do j = 0, 2 ** (m-1) - 1

			! The main difference between fft and ifft is the sign of this
			! exponent
			e = exp(sign_ * IMAG * PI * j / 2 ** m)

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

subroutine tridiag_factor(a)
	! This is the factorization phase of tridiag_invmul(), without pivoting
	!
	! If you are only solving once, the interface for tridiag_invmul() is
	! easier.  However, separate factor and solve phases are better if you need
	! to factor once and then solve multiple times without repeating the
	! factorization
	!
	! Compare the interface of lapack routines dgetrf() (factor) and dgetrs()
	! (solve) vs dgesv() (combined factor-solve).  For example:
	!
	!     https://github.com/JeffIrwin/ribbit/blob/ec868aa3db96258b95909f3101434128fc42428f/src/ribbit.f90#L1722

	double precision, intent(inout) :: a(:,:)
	!********
	integer :: i, n

	n = size(a, 2)  ! *not* same as size 1, which is always 3 for this tridiag storage scheme!

	! TODO: check that size 1 is 3?

	! This implementation is based on wikipedia's "Tridiagonal
	! matrix algorithm" page.  Note that variable names and
	! subscripts differ.

	! l(1) = a(1,1) is 0, u(n) = a(3,n) is 0

	! TODO: check a(1,1) == 0 and a(3,n) == 0?

	a(3,1) = a(3,1) / a(2,1)
	do i = 2, n - 1
		a(2, i) = a(2, i) - a(1, i) * a(3, i-1)
		a(3, i) = a(3, i) / a(2, i)
	end do

end subroutine tridiag_factor

!********

subroutine tridiag_solve(a, bx)
	! This is the solving phase of tridiag_invmul()

	double precision, intent(in) :: a(:,:)
	double precision, intent(inout) :: bx(:)  ! could be extended to rank-2 for multiple RHS's

	!********

	integer :: i, n

	n = size(a, 2)  ! *not* same as size 1, which is always 3 for this tridiag storage scheme!

	! Forward sweep
	bx(1)  = bx(1)  / a(2,1)
	do i = 2, n - 1
		bx(i) = (bx(i) - a(1, i) * bx(i-1)) / a(2, i)
	end do
	bx(n) = (bx(n) - a(1,n) * bx(n-1)) / (a(2,n) - a(1,n) * a(3,n-1))

	! Back substitution
	do i = n - 1, 1, -1
		bx(i) = bx(i) - a(3,i) * bx(i+1)
	end do

end subroutine tridiag_solve

!********

subroutine tridiag_invmul(a, bx)
	! Solve the linear algebra problem for `x`:
	!
	!     a * x = b
	!
	! Tridiagonal matrix `a` is packed into a 3xn array.  On input, `bx` is the
	! RHS `b`.  On output, `bx` is the solution `x`
	!
	! See exercises.f90 for an example of how to pack the `a` matrix

	double precision, intent(inout) :: a(:,:)
	double precision, intent(inout) :: bx(:)  ! could be extended to rank-2 for multiple RHS's

	! I'm not sure if it's better for performance to use the current
	! representation of `a` or its transpose instead.  Most expressions in
	! factor/solve work on a 3x2 or 2x1 submatrix of `a` currently, so it might
	! be better as-is, but not much of a difference

	!********

	call tridiag_factor(a)
	call tridiag_solve(a, bx)

end subroutine tridiag_invmul

!===============================================================================
subroutine tridiag_corner_invmul(aa, bx)

	! Solve an augmented tridiagonal system with off-diagonal opposite corner
	! elements, AKA the cyclic Thomas algorithm
	!
	! The cyclic code on the Tridiagonal_matrix_algorithm wikipedia page is
	! straight up wrong.  Here are better references:
	!
	!   - https://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
	!   - https://sachinashanbhag.blogspot.com/2014/03/cyclic-tridiagonal-matrices.html

	double precision, intent(inout) :: aa(:,:)
	double precision, intent(inout) :: bx(:)

	!********
	double precision :: a1, b1, cn, v1, vn
	double precision, allocatable :: u(:)

	integer :: i, n

	! Make a modified tridiagonal system with no corners and modified first and
	! last main diagonal elements
	n = size(aa, 2)

	a1 = aa(1, 1)  ! top right elem
	b1 = aa(2, 1)  ! first main diag elem
	cn = aa(3, n)  ! lower left elem

	! Make vectors `u` and `v` for Sherman-Morrison formula
	!
	! We need a scratch vector for u.  Although it is initially sparse, we solve
	! `aa \ u` later which fills it in.  For `v` vector we only need v(1) and
	! v(n)
	u = [(0.d0, i = 1, n)]
	u(1) = -b1;  u(n) = cn
	v1 = 1.d0 ;  vn = -a1 / b1

	aa(1,1) = 0.d0
	aa(3,n) = 0.d0
	aa(2,1) = aa(2,1) - u(1) * v1
	aa(2,n) = aa(2,n) - u(n) * vn

	!print *, "aa = "
	!print "(3es16.6)", aa

	! Solve two pure tridiag problems
	call tridiag_factor(aa)
	call tridiag_solve(aa, bx)
	call tridiag_solve(aa, u)

	!print *, "bx = ", bx
	!print *, "u = ", u

	! Solution to cylcic problem
	bx = bx - ((v1*bx(1) + vn*bx(n)) / (1.d0 + (v1*u(1) + vn*u(n)))) * u

	!print *, "solution = ", bx

end subroutine tridiag_corner_invmul

!===============================================================================

subroutine banded_factor(a, al, indx, nl, nu)
	! This is the factorization phase of banded_invmul(), with pivoting
	!
	! Compare the interface of lapack routines dgbtrf() (factor) and dgbtrs()
	! (solve) vs dgbsv() (combined factor-solve)

	use numa__utils

	! Reconsider the order of arguments -- maybe in(out) first and out last.  Or
	! perhaps it's better to leave as-is for consistency with banded_solve()
	! args
	double precision, intent(inout) :: a(:,:)
	double precision, intent(out), allocatable :: al(:,:)
	integer, intent(out), allocatable :: indx(:)
	integer, intent(in) :: nl, nu

	!********
	double precision :: d, dum
	integer :: i, j, k, l, n, mm

	!print *, "nl, nu = ", nl, nu

	n = size(a, 2)  ! *not* same as size 1

	if (nl + nu + 1 /= size(a, 1)) then
		! TODO: return code?  This should probably be a fatal error
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": bad size of matrix `a` in banded_factor()"
	end if

	! l(1) = a(1,1) is 0, u(n) = a(3,n) is 0

	! TODO: check a(1,1) == 0 and a(3,n) == 0?
	mm = nl + 1 + nu

	! Left-shift the top nl rows by the number of zeroes
	do i = 1, nl
		a(1: mm-nl+i-1, i) = a(1+nl-i+1: mm, i)
		a(mm - nl + i: mm, i) = 0
	end do
	d = 1  ! not used for anything, but it shows even vs odd row pivot permutations

	!print *, "a first loop = "
	!print "(5es18.6)", a

	allocate(indx(n))
	allocate(al(nl, n))
	l = nl
	do k = 1, n
		dum = a(1, k)
		i = k
		if (l < n) l = l + 1

		do j = k+1, l
			if (abs(a(1, j)) > abs(dum)) then
				dum = a(1, j)
				i = j
			end if
		end do

		indx(k) = i
		if (i /= k) then
			d = -d
			a(:, [i, k]) = a(:, [k, i])
		end if

		do i = k+1, l
			dum = a(1,i) / a(1,k)
			al(i-k, k) = dum
			a(1:mm-1, i) = a(2:mm, i) - dum * a(2:mm, k)
			a(mm, i) = 0
		end do
	end do
	!print *, "a = "
	!print "(5es18.6)", a
	!print *, "al = "
	!print "(2es18.6)", al

end subroutine banded_factor

!********

subroutine banded_solve(a, al, indx, nl, nu, bx)

	! This is the solving phase of banded_invmul()

	double precision, intent(in) :: a(:,:)
	double precision, intent(in) :: al(:,:)
	integer, intent(in) :: indx(:)

	double precision, intent(inout) :: bx(:)  ! could be extended to rank-2 for multiple RHS's

	integer, intent(in) :: nl, nu

	!********
	double precision :: dum
	integer :: i, j, k, l, n, mm

	!********
	! Solve stage

	n = size(a, 2)  ! *not* same as size 1
	mm = nl + 1 + nu

	l = nl
	do k = 1, n
		j = indx(k)
		if (j /= k) bx([j, k]) = bx([k, j])
		if (l < n) l = l + 1
		bx(k+1: l) = bx(k+1: l) - al(1: l-k, k) * bx(k)
	end do

	l = 1
	do i = n, 1, -1
		dum = bx(i) - dot_product(a(2:l, i), bx(i+1: l+i-1))
		bx(i) = dum / a(1, i)
		if (l < mm) l = l + 1
	end do

end subroutine banded_solve

!********

subroutine banded_invmul(a, bx, nl, nu)
	! Solve the linear algebra problem for `x`:
	!
	!     a * x = b
	!
	! Banded matrix `a` is packed.  TODO: notes on packed representation
	!
	! TODO: should {a, nl, nu} be a struct?  One of the n*'s could be redundant

	double precision, intent(inout) :: a(:,:)
	double precision, intent(inout) :: bx(:)  ! could be extended to rank-2 for multiple RHS's
	integer, intent(in) :: nl, nu

	!********

	double precision, allocatable :: al(:,:)
	integer, allocatable :: indx(:)

	call banded_factor(a, al, indx, nl, nu)
	call banded_solve (a, al, indx, nl, nu, bx)

end subroutine banded_invmul

!===============================================================================

subroutine lu_invmul(a, bx)
	! Solve the linear algebra problem for `x`:
	!
	!     a * x = b

	double precision, intent(inout) :: a(:,:)
	double precision, intent(inout) :: bx(:)  ! could be extended to rank-2 for multiple RHS's

	!********

	integer :: i
	integer, allocatable :: pivot(:)

	! Initialize pivot to identity
	pivot = [(i, i = 1, size(a,1))]

	call lu_factor(a, pivot)

	!print *, "pivot = ", pivot
	!print *, "lu_factor(a) = "
	!!print "(5es18.6)", a
	!print "(6es15.5)", transpose(a)

	call lu_solve(a, bx, pivot)

end subroutine lu_invmul

!********

subroutine lu_factor(a, pivot)

	!! Make pivoting optional (for comparison to other solvers)
	!logical, parameter :: DO_PIVOT = .false.

	double precision, intent(inout) :: a(:,:)
	integer, intent(inout) :: pivot(:)

	!********

	integer :: i, j, k, n, max_index

	!print *, "pivot = ", pivot

	n = size(a, 1)

	do i = 1, n
		! Find max value in column i
		max_index = i
		!if (DO_PIVOT) then
			do j = i+1, n
				if (abs(a(j,i)) > abs(a(i,i))) max_index = j
			end do
		!end if

		! Swap rows
		pivot([i, max_index]) = pivot([max_index, i])

		do j = i+1, n
			a    (pivot(j), i) = &
				a(pivot(j), i) / &
				a(pivot(i), i)
			do k = i+1, n
				a    (pivot(j), k) = &
					a(pivot(j), k) - &
					a(pivot(j), i) * &
					a(pivot(i), k)
			end do
		end do

	end do

end subroutine lu_factor

!********

subroutine lu_solve(a, bx, pivot)

	double precision, intent(in) :: a(:,:)
	double precision, intent(inout) :: bx(:)
	integer, intent(in) :: pivot(:)

	integer :: i, j, n

	!print *, "pivot lu_solve = ", pivot

	n = size(a, 1)
	!print *, "bx init = ", bx

	! Forward substitution
	do i = 1, n
		do j = 1, i-1
			bx    (pivot(i)) = &
				bx(pivot(i)) - &
				a (pivot(i), j) * &
				bx(pivot(j))
		end do
	end do
	!print *, "bx forward = ", bx

	! Back substitution
	do i = n, 1, -1
		do j = i+1, n
			bx    (pivot(i)) = &
				bx(pivot(i)) - &
				a (pivot(i), j) * &
				bx(pivot(j))
		end do
		bx    (pivot(i)) = &
			bx(pivot(i)) / &
			a (pivot(i), i)
	end do
	!print *, "bx back = ", bx

	! Unpivot
	bx = bx(pivot)
	!print *, "bx unpivot = ", bx

end subroutine lu_solve

!===============================================================================

function spline_no_curve(xi, yi, x) result(fx)
	! This function finds the cubic spline with the "no curvature" boundary
	! condition, i.e. the 2nd derivative is 0 at both endpoints, or case (a) in
	! the text
	!
	! This is also known as the natural spline

	double precision, intent(in)  :: xi(:)
	double precision, intent(in)  :: yi(:)
	double precision, intent(in)  :: x(:)
	double precision, allocatable :: fx(:)

	fx = spline_general(xi, yi, x, SPLINE_CASE_NO_CURVE, 0.d0, 0.d0)

end function spline_no_curve

!********

function spline_periodic(xi, yi, x) result(fx)
	! This function finds the cubic spline with periodic boundary conditions,
	! i.e. the 1st and 2nd derivates match at the start and end points, case (b)
	! in the text

	double precision, intent(in)  :: xi(:)
	double precision, intent(in)  :: yi(:)
	double precision, intent(in)  :: x(:)
	double precision, allocatable :: fx(:)

	fx = spline_general(xi, yi, x, SPLINE_CASE_PERIODIC, 0.d0, 0.d0)

end function spline_periodic

!********

function spline_prescribed(xi, yi, x, dy_start, dy_end) result(fx)
	! This function finds the cubic spline with prescribed first derivatives at
	! the end points, i.e. case (c) in the text

	double precision, intent(in)  :: xi(:)
	double precision, intent(in)  :: yi(:)
	double precision, intent(in)  :: x(:)
	double precision, intent(in)  :: dy_start, dy_end
	double precision, allocatable :: fx(:)

	fx = spline_general(xi, yi, x, SPLINE_CASE_PRESCRIBED, dy_start, dy_end)

end function spline_prescribed

!===============================================================================

function spline_general(xi, yi, x, case_, dy_start, dy_end) result(fx)

	use numa__utils

	double precision, intent(in)  :: xi(:)
	double precision, intent(in)  :: yi(:)
	double precision, intent(in)  :: x(:)

	double precision, intent(in)  :: dy_start, dy_end
	integer, intent(in) :: case_

	double precision, allocatable :: fx(:)
	!********

	double precision, allocatable :: h(:), lambda(:), mu(:), &
		d(:), a(:,:), aj, bj
	integer :: i, j, n

	if (.not. any(case_ == &
		[SPLINE_CASE_NO_CURVE, SPLINE_CASE_PERIODIC, SPLINE_CASE_PRESCRIBED])) then
		write(*,*) RED // "Error" // COLOR_RESET // &
			": bad spline case in spline_general()"
		call exit(1)
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

		call tridiag_corner_invmul(a, d)

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

		call tridiag_invmul(a, d)
		!print *, "moments = ", d

	end if

	allocate(fx( size(x) ))

	if (x(1) < xi(1)) then
		! TODO: all these warnings should probably be fatal
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": `x` underflows first spline control point"
	end if
	if (x(size(x)) > xi(size(xi))) then
		write(*,*) YELLOW // "Warning" // COLOR_RESET // &
			": `x` overflows last spline control point"
	end if

	! Evaluate the spline
	j = 1
	aj = (yi(j+1) - yi(j)) / h(j) - h(j) / 6.d0 * (d(j+1) - d(j))
	bj = yi(j) - d(j) * (h(j) ** 2) / 6.d0
	do i = 1, size(x)

		if (i > 1) then
			if (x(i) <= x(i-1)) then
				write(*,*) YELLOW // "Warning" // COLOR_RESET // &
					": `x` is not sorted ascending in spline_general()"
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

function cardinal_spline(xc, t, tension) result(x)

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

	double precision, intent(in) :: xc(:,:)
	double precision, intent(in) :: t(:)
	double precision, intent(in) :: tension

	double precision, allocatable :: x(:,:)
	!********

	double precision :: t_, tmod
	double precision, allocatable :: xcl(:,:)
	double precision, allocatable :: v(:,:), v0(:,:)

	integer :: i, j, it, nd, nc, nt, ncl, segment

	! TODO: panic if tension <= 0

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

double precision function romberg_integrator(f, xmin, xmax, tol, nmax) result(area)
	! Integrate `f` from xmin to xmax using Romberg integration with tolerance
	! `tol`

	use numa__utils
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax
	double precision, intent(in) :: tol
	integer, optional, intent(in) :: nmax
	!********

	double precision :: dx, s, q, area0
	double precision, allocatable :: t(:)
	integer :: i, k, nt, n
	logical :: converged

	n = 10
	if (present(nmax)) n = nmax

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

	double precision :: half_h, mid, offset

	half_h = 0.5d0 * (xmax - xmin)
	mid = xmin + half_h
	offset = sqrt(3.d0) / 3.d0 * half_h

	area = half_h * (f(mid - offset) + f(mid + offset))

end function gauss2_single

!===============================================================================

double precision function gauss3_single(f, xmin, xmax) result(area)
	! Integrate `f` from xmin to xmax using 3-point Gaussian integration with
	! only a single interval

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax
	!********

	double precision :: half_h, mid, offset

	half_h = 0.5d0 * (xmax - xmin)
	mid = xmin + half_h
	offset = 0.7745966692414834d0 * half_h

	area = half_h * ( &
		5.d0/9 * f(mid - offset) + &
		8.d0/9 * f(mid) + &
		5.d0/9 * f(mid + offset) &
	)

end function gauss3_single

!===============================================================================

double precision function gauss4_single(f, xmin, xmax) result(area)
	! Integrate `f` from xmin to xmax using 4-point Gaussian integration with
	! only a single interval

	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax
	!********

	double precision :: half_h, mid

	double precision, parameter :: &
		x1 = 0.8611363115940526d0, &
		x2 = 0.3399810435848563d0
	double precision, parameter :: &
		w1 = 0.3478548451374538d0, &
		w2 = 0.6521451548625461d0

	half_h = 0.5d0 * (xmax - xmin)
	mid = xmin + half_h

	area = half_h * ( &
		w1 * f(mid - x1 * half_h) + &
		w1 * f(mid + x1 * half_h) + &
		w2 * f(mid - x2 * half_h) + &
		w2 * f(mid + x2 * half_h)   &
	)

end function gauss4_single

!===============================================================================

end module numa

