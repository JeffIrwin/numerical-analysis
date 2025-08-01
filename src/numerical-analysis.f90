
module numa

	implicit none

	! TODO:
	! - There are several places where I note a routine should "panic" or issue
	!   a "fatal" error, only logging a warning for now and continuing.  These
	!   should actually stop (or wrap a stop), unless an optional `iostat` arg
	!   is given to the function, in which case we can continue and let the
	!   caller decide how to recover.  I guess this is the Fortranic way

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

	interface invmul
		! The built-in `matmul()` works on matrices and vectors similarly
		!
		! TODO: make a fn version.  Name?
		procedure :: lu_invmul_vec
		procedure :: lu_invmul_mat
	end interface invmul

	interface zeros
		procedure :: zeros_vec
		procedure :: zeros_mat
	end interface zeros

	interface diag
		procedure :: diag_set
		procedure :: diag_get
	end interface diag

	! If this file gets too long, it might be good to split it up roughly
	! per-chapter, e.g. into interpolate.f90, (fft.f90,) integrate.f90, etc.

	!********

	private :: simpson_adapt_aux, gk15_aux, gk15i_aux

	!********

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
			e = exp(sgn * IMAG * PI * j / 2 ** m)

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

subroutine lu_invmul_vec(a, bx)
	! Solve the linear algebra problem for `x`:
	!
	!     a * x = b

	double precision, intent(inout) :: a(:,:)
	double precision, intent(inout) :: bx(:)

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

end subroutine lu_invmul_vec

!********

subroutine lu_invmul_mat(a, bx)
	! Solve the linear algebra problem for `x`:
	!
	!     a * x = b

	double precision, intent(inout) :: a(:,:)
	double precision, intent(inout) :: bx(:,:)

	!********

	integer :: i, nrhs
	integer, allocatable :: pivot(:)

	! Initialize pivot to identity
	pivot = [(i, i = 1, size(a,1))]

	call lu_factor(a, pivot)

	!print *, "pivot = ", pivot
	!print *, "lu_factor(a) = "
	!!print "(5es18.6)", a
	!print "(6es15.5)", transpose(a)

	nrhs = size(bx, 2)
	do i = 1, nrhs
		call lu_solve(a, bx(:,i), pivot)
	end do

end subroutine lu_invmul_mat

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

subroutine cholesky_factor(a)
	! Cholesky factorization is only possible for positive definite matrices
	! `a`!  There is no enforcement here

	double precision, intent(inout) :: a(:,:)

	!********

	double precision :: x, p
	integer :: i, j, k, n

	n = size(a, 1)

	do i = 1, n
	do k = i, n
		x = a(i,k)
		do j = i-1, 1, -1
			x = x - a(k,j) * a(i,j)
		end do
		if (i == k) then
			if (x <= 0.d0) then
				! TODO: panic
				print *, "Error: singular matrix in cholesky_factor()"
			end if
			p = 1.d0 / sqrt(x)
		end if
		a(k,i) = x * p
	end do
	end do

end subroutine cholesky_factor

!********
double precision function sign_(x)
	double precision, intent(in) :: x
	sign_ = sign(1.d0, x)
end function sign_

!********

subroutine hess(a)
	! Apply the Householder Hessenberg reduction to `a`
	!
	! Source:  https://dspace.mit.edu/bitstream/handle/1721.1/75282/18-335j-fall-2006/contents/lecture-notes/lec14.pdf
	double precision, intent(inout) :: a(:,:)
	!********

	integer :: k, m
	double precision, allocatable :: x(:), v(:)

	m = size(a,1)

	do k = 1, m-2

		x = a(k+1:, k)
		v = x
		v(1) = v(1) + sign_(x(1)) * norm2(x)
		v = v / norm2(v)

		! TODO: avoid outer_product() at least, and maybe `x` temp array.
		! Avoiding `v` temp array might be a little more work than worthwhile

		a(k+1:, k:) = a(k+1:, k:) - 2 * outer_product(v, matmul(v, a(k+1:, k:)))
		a(:, k+1:)  = a(:, k+1:)  - 2 * outer_product(matmul(a(:, k+1:), v), v)

	end do

end subroutine hess

!********

subroutine qr_factor(a, diag_)
	! Replace `a` with its QR factorization using Householder transformations
	double precision, intent(inout) :: a(:,:)
	double precision, allocatable, intent(out) :: diag_(:)

	!********

	double precision :: s, normx, u1, wa
	integer :: j, k, n

	n = size(a, 1)
	allocate(diag_(n))
	diag_ = 0.d0

	! Ref:  https://www.cs.cornell.edu/~bindel/class/cs6210-f09/lec18.pdf
	do j = 1, n

		normx = norm2(a(j:, j))
		! TODO: panic if normx == 0

		s = -sign_(a(j,j))
		u1 = a(j,j) - s * normx
		a(j+1:, j) = a(j+1:, j) / u1
		a(j,j) = s * normx
		diag_(j) = -s * u1 / normx

		do k = j+1, n

			! TODO: for Hessenberg `a`, dot_product() can be applied to a
			! smaller slice here
			wa = diag_(j) * (dot_product(a(j+1:, j), a(j+1:, k)) + a(j,k))

			a(j,k) = a(j,k) - wa
			a(j+1:, k) = a(j+1:, k) - a(j+1:, j) * wa
		end do

	end do
	!print *, "diag_ = "
	!print "(es15.5)", diag_

	! In the end, `R` is formed by zeroing appropriate elements of `a`
	!
	! And then `Q` can be expensively calculated as inv(r) * a (or some similar
	! transpose), or more efficiently like in the loop below
	!
	! Results can be checked by verifying that q' * q == eye() or that a = q * r

end subroutine qr_factor

!********

function qr_get_q_expl(qr, diag_) result(q)
	! Explicitly get the unitary Q matrix from a previously QR-decomposed
	! matrix.  In most cases you will want to avoid using this fn and just
	! implicitly multiply something by Q instead using qr_mul()
	double precision, intent(in) :: qr(:,:), diag_(:)
	double precision, allocatable :: q(:,:)

	q = qr_mul(qr, diag_, eye(size(qr,1)))

end function qr_get_q_expl

function qr_get_r_expl(qr) result(r)
	! Explicitly get R by zeroing lower triangle
	double precision, intent(in) :: qr(:,:)
	double precision, allocatable :: r(:,:)
	!********
	integer :: i

	r = qr
	do i = 1, size(qr, 1)
		r(i+1:, i) = 0
	end do

end function qr_get_r_expl

!********

function qr_mul(qr, diag_, x) result(qx)
	! Implicitly multiply Q * x with a previously computed QR factorization `qr`
	double precision, intent(in) :: qr(:,:), diag_(:), x(:,:)
	double precision, allocatable :: qx(:,:)
	!********
	double precision :: wq
	integer :: j, k, n

	!print *, "qr = "
	!print "(5es15.5)", qr

	n = size(qr, 1)
	qx = x  ! could make a subroutine version which replaces x instead

	! To multiply by transpose(Q) instead, just loop from 1 up to n instead
	do j = n, 1, -1
		do k = 1, n
			wq = diag_(j) * (dot_product(qr(j+1:, j), qx(j+1:, k)) + qx(j,k))
			qx(j, k) = qx(j, k) - wq
			qx(j+1:, k) = qx(j+1:, k) - qr(j+1:, j) * wq
		end do
	end do

	!print *, "qx = "
	!print "(5es15.5)", qx

end function qr_mul

!********

function outer_product(a, b) result(c)
	double precision, intent(in) :: a(:), b(:)
	double precision, allocatable :: c(:,:)

	integer :: i, j, na, nb
	na = size(a)
	nb = size(b)

	allocate(c(na, nb))
	do j = 1, nb
	do i = 1, na
		c(i,j) = a(i) * b(j)
	end do
	end do

end function outer_product

!********

subroutine cholesky_solve(a, bx)

	double precision, intent(in) :: a(:,:)
	double precision, intent(inout) :: bx(:)

	integer :: i, j, n

	n = size(a, 1)
	!print *, "bx init = ", bx

	! This is the same as lu_solve(), except:
	! - The pivot is the identity
	! - One of the `a` terms is transposed, because U == transpose(L)
	! - Forward substitution has an extra division step because U has non-zero diagonals

	! Forward substitution
	do i = 1, n
		do j = 1, i-1
			bx    (i) = &
				bx(i) - &
				a (i, j) * &
				bx(j)
		end do
		bx    (i) = &
			bx(i) / &
			a (i, i)
	end do
	!print *, "bx forward = ", bx

	! Back substitution
	do i = n, 1, -1
		do j = i+1, n
			bx    (i) = &
				bx(i) - &
				a (j, i) * &
				bx(j)
		end do
		bx    (i) = &
			bx(i) / &
			a (i, i)
	end do
	!print *, "bx back = ", bx

	! No pivoting for Cholesky

end subroutine cholesky_solve

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

double precision function romberg_integrator(f, xmin, xmax, tol, max_levels) &
		result(area)
	! Integrate `f` from xmin to xmax using Romberg integration with tolerance
	! `tol`
	!
	! TODO: make tol optional?

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
		f, xmin, xmax, tol, max_levels &
	) &
	result(area)

	! Integrate `f` from xmin to xmax using the adaptive Gauss-Kronrod method,
	! with 7-point Gauss and 15-point Kronrod

	use numa__utils
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax, tol
	integer, optional, intent(in) :: max_levels
	!********

	integer :: n, neval

	logical :: is_eps_underflow, is_max_level

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
		f, xmin, tol, max_levels &
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
	!********

	double precision, parameter :: split_ = 1.d0

	integer :: n, neval

	logical :: is_eps_underflow, is_max_level

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
		f, tol, max_levels &
	) &
	result(area)

	! Integrate `f` from -infinity to infinity using the adaptive Gauss-Kronrod method,
	! with 7-point Gauss and 15-point Kronrod

	use numa__utils
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: tol
	integer, optional, intent(in) :: max_levels
	!********

	double precision, parameter :: left = -1.d0, right = 1.d0

	integer :: n, neval

	logical :: is_eps_underflow, is_max_level

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
		f, xmax, tol, max_levels &
	) &
	result(area)

	! Integrate `f` from -infinity to xmax using the adaptive Gauss-Kronrod method,
	! with 7-point Gauss and 15-point Kronrod

	use numa__utils
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmax, tol
	integer, optional, intent(in) :: max_levels
	!********

	double precision, parameter :: split_ = 1.d0

	integer :: n, neval

	logical :: is_eps_underflow, is_max_level

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
		f, xmin, xmax, tol, max_levels &
	) &
	result(area)

	! Integrate `f` from xmin to xmax using the adaptive Simpson method

	use numa__utils
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: xmin, xmax, tol
	integer, optional, intent(in) :: max_levels
	!********

	double precision :: h, s, fa, fb, fm
	integer :: n, neval
	logical :: is_eps_underflow, is_max_level

	area = 0.d0
	h = xmax - xmin
	if (h == 0) return

	n = 16
	if (present(max_levels)) n = max_levels

	! TODO: panic if tol == 0.  Same for gk15*_adaptive_integrator.  Recursive core will barf anyway

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

function eye(n)
	! n x n identity matrix

	integer, intent(in) :: n
	double precision, allocatable :: eye(:,:)

	integer :: i, j
	allocate(eye(n, n))
	do i = 1, n
		do j = 1, n
			if (i == j) then
				eye(i,j) = 1.d0
			else
				eye(i,j) = 0.d0
			end if
		end do
	end do

end function eye

!********

function zeros_mat(m, n) result(a)
	! m x n matrix of 0
	integer, intent(in) :: m, n
	double precision, allocatable :: a(:,:)
	allocate(a(m, n))
	a = 0
end function zeros_mat

function zeros_vec(n) result(a)
	! Size n vector of 0
	integer, intent(in) :: n
	double precision, allocatable :: a(:)
	allocate(a(n))
	a = 0
end function zeros_vec

!********

function diag_set(v) result(d)
	! Spread a diagonal vector `v` into a matrix
	double precision, intent(in) :: v(:)
	double precision, allocatable :: d(:,:)
	!********
	integer :: i, n

	n = size(v)
	d = zeros(n, n)
	do i = 1, n
		d(i,i) = v(i)
	end do

end function diag_set

function diag_get(a) result(v)
	! Get the diagonal vector `v` from a matrix `a`
	double precision, intent(in) :: a(:,:)
	double precision, allocatable :: v(:)
	!********
	integer :: i, n

	n = min(size(a,1), size(a,2))
	allocate(v(n))
	do i = 1, n
		v(i) = a(i,i)
	end do

end function diag_get

!===============================================================================

function inv(a)
	! Return the inverse of matrix `a` without modifying `a`
	double precision, allocatable, intent(in) :: a(:,:)
	double precision, allocatable :: inv(:,:)

	!double precision, allocatable :: copy_(:,:)

	!! Use LU decomposition
	!copy_ = a
	!inv = eye(size(a,1))
	!call invmul(copy_, inv)

	! Use Gauss-Jordan.  I suspect it does less copying than my LU inverse
	! implementation.  It's around 40% faster, although both methods are O(n**3)
	inv = a
	call gauss_jordan(inv)

end function inv

!===============================================================================

subroutine invert(a)
	! Replace `a` with its inverse
	double precision, allocatable, intent(inout) :: a(:,:)
	!double precision, allocatable :: copy_(:,:)

	!********
	! Use Gauss-Jordan
	call gauss_jordan(a)

	!********
	! Use LU

	!! One move and one copy
	!call move_alloc(a, copy_)
	!a = eye(size(copy_, 1))
	!call invmul(copy_, a)

	!!! Two copies
	!!copy_ = a
	!!a = eye(size(a,1))
	!!call invmul(copy_, a)

	!!! Also two copies
	!!copy_ = eye(size(a,1))
	!!call invmul(a, copy_)
	!!a = copy_

	!!! Break tests
	!!a(1,1) = 1.01d0 * a(1,1)

end subroutine invert

!===============================================================================

subroutine gauss_jordan(a)
	! Use the Gauss-Jordan method to replace matrix `a` with its inverse
	double precision, allocatable, intent(inout) :: a(:,:)

	!********

	double precision :: max_, hr

	integer, allocatable :: p(:)  ! pivot
	integer :: i, j, k, n, r, r1(1)

	n = size(a, 1)
	p = [(i, i = 1, n)]

	do j = 1, n

		r1 = maxloc(abs(a(j+1: n, j))) + j
		r = r1(1)
		max_ = abs(a(r, j))
		if (max_ == 0) then
			! TODO: panic
			print *, "Error: singular matrix in gauss_jordan()"
		end if

		if (r > j) then
			! Swap rows
			a([j, r], :) = a([r, j], :)
			p([j, r]) = p([r, j])
		end if

		! Transform
		hr = 1.d0 / a(j,j)
		a(:,j) = a(:,j) * hr
		a(j,j) = hr
		do k = 1, n
			if (k == j) cycle
			do i = 1, n
				if (i == j) cycle

				a(i,k) = a(i,k) - a(i,j) * a(j,k)
			end do
			a(j,k) = -hr * a(j,k)
		end do

	end do

	! Swap columns
	a(:,p) = a

	!print *, "In Gauss-Jordan:"
	!print *, "a = "
	!print "(5es18.6)", a

end subroutine gauss_jordan

!===============================================================================

function eig_basic_qr(a, iters) result(eigvals)
	! Get the real eigenvalues of `a` using `iters` iterations of the basic QR
	! algorithm.  Many iterations may be required for good convergence with
	! large sized `a`
	!
	! The matrix `a` is modified in the process by QR decomposition
	!
	! This only supports real eigenvalues.  The same algorithm could be adapted
	! for complex eigenvalues:
	!
	!     https://math.stackexchange.com/a/3072781/232771
	!
	! The basic idea is that `a` converges to a "quasi-triangular" real matrix
	! in this algorithm.  Where the off-triangle approaches 0, we have real
	! eigenvalues.  Where there is a non-zero in the off-diagonal, we have a
	! complex conjugate pair of eigenvalues, which could be computed as the same
	! eigenvalue of its 2x2 sub-matrix block

	use numa__utils, only:  sorted
	double precision, intent(inout) :: a(:,:)
	double precision, allocatable :: eigvals(:)
	integer, intent(in) :: iters
	!********

	double precision, allocatable :: diag_(:), r(:,:), q(:,:)
	integer :: i

	do i = 1, iters
		call qr_factor(a, diag_)
		r = qr_get_r_expl(a)

		! Could also do a = transpose(qr_mul_transpose(), transpose(r)),
		! but two transposes seems expensive
		q = qr_get_q_expl(a, diag_)
		a = matmul(r, q)
		! TODO: naive matmul does ~2x as many operations as needed here, given
		! that r is half zeros.  Should roll my own triangular matmul.  Then
		! calling qr_get_r_expl() can also be eliminated as we would just
		! implicitly ignore zeros, using `a` directly instead of `r`

	end do
	!print *, "a = "
	!print "(4es15.5)", a

	eigvals = sorted(diag(a))

end function eig_basic_qr

!===============================================================================

end module numa

