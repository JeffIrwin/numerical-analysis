
module numa__exercises

	use numa
	use numa__functions
	use numa__utils

	implicit none

contains

!===============================================================================

integer function assert(condition)
	! Return 0 if true or 1 if false
	logical, intent(in) :: condition
	if (condition) then
		assert = 0
	else
		assert = 1
	end if
end function assert

!===============================================================================

subroutine test(val, expect, tol, nfail, msg)
	double precision, intent(in) :: val, expect, tol
	integer, intent(inout) :: nfail
	character(len = *), intent(in) :: msg

	!********
	integer :: io

	io = assert(abs(val - expect) < tol)
	nfail = nfail + io
	if (io /= 0) write(*,*) ERROR // "assert failed for value " // &
		to_str(val) // " in "// msg

end subroutine test

!********

subroutine tests(vals, expects, tol, nfail, msg)
	double precision, intent(in) :: vals(:), expects(:), tol
	integer, intent(inout) :: nfail
	character(len = *), intent(in) :: msg
	!********
	integer :: i
	do i = 1, size(vals, 1)
		call test(vals(i), expects(i), tol, nfail, msg)
	end do
end subroutine tests

!===============================================================================

integer function chapter_2_exercise_2() result(nfail)

	! 2. Interpolate the function ln x by a quadratic polynomial at x = 10, 11,
	!    12.
	!    (a) Estimate the error commited for x = 11.1 when approximating ln x by
	!        the interpolating polynomial

	double precision, parameter :: tol = 1.d-3
	!double precision, parameter :: tol = 4.d-5
	double precision :: x, fx
	double precision, allocatable :: xi(:), xs(:), fxs(:)

	write(*,*) CYAN // "Starting chapter_2_exercise_2()" // COLOR_RESET

	nfail = 0

	! Support points `xi`
	xi = [10.d0, 11.d0, 12.d0]

	x = 11.1d0

	! The acceptable tolerance depends on xi, x, and the fn being interpolated

	fx = lagrange_interpolator(xi, log_fn, x)
	!nfail = nfail + assert(abs(fx - log_fn(x)) < tol)
	call test(fx, log_fn(x), tol, nfail, "lagrange_interpolator()")
	print *, "fx   = ", fx, " (Lagrange)"
	print *, "f(x) = ", log_fn(x)
	print *, ""

	!fx = lagrange_interpolator(xi, exp_fn, x)  ! not from the text book, but let's see what happens
	!print *, "fx   = ", fx
	!print *, "f(x) = ", exp_fn(x)
	!print *, ""

	fx = neville_interpolator(xi, log_fn, x)
	!fx = neville_interpolator_vals(xi, log(xi), x)
	call test(fx, log_fn(x), tol, nfail, "neville_interpolator()")
	print *, "fx   = ", fx, " (Neville)"
	print *, "f(x) = ", log_fn(x)
	print *, ""

	xs = [11.1d0, 11.2d0]
	!xs = [11.1d0, 11.2d0, 11.3d0, 11.4d0]
	fxs = newton_interpolator(xi, log_fn, xs)
	!fxs = newton_interpolator_vals(xi, log(xi), xs)
	call tests(fxs, log(xs), tol, nfail, "newton_interpolator()")
	print *, "fxs   = ", fxs, " (Newton)"
	print *, "f(xs) = ", log(xs)
	print *, ""

end function chapter_2_exercise_2

!===============================================================================

integer function chapter_2_example_2p2p4() result(nfail)

	! Compare rational and polynomial interpolation, from the example in section
	! 2.2.4 on page 73
	!
	! For the function f(x) = cot(x) th evalues at [1: 5] degrees have been
	! tabulated.  Determine an approximate value for cot(2.5 degrees) using both
	! polynomial and rational interpolation

	double precision, parameter :: tol = 1.d-4
	double precision :: x, fx
	double precision, allocatable :: xi(:)

	write(*,*) CYAN // "Starting chapter_2_example_2p2p4()" // COLOR_RESET

	nfail = 0

	! Support points `xi`
	xi = [1, 2, 3, 4, 5]

	x = 2.5d0

	! The acceptable tolerance depends on xi, x, and the fn being interpolated

	! The expected value for polynomial interpolation, 22.6352, is quite far
	! from the actual cotangent value because cot has a singularity at x == 0.
	! Rational interpolation works better near singularities

	fx = neville_interpolator(xi, cotd_fn, x)
	call test(fx, 22.6352d0, tol, nfail, "neville_interpolator()")
	print *, "fx   = ", fx, " (Neville polynomial)"
	print *, "f(x) = ", cotd_fn(x)
	print *, ""

	fx = neville_rational_interpolator(xi, cotd_fn, x)
	call test(fx, cotd_fn(x), tol, nfail, "neville_rational_interpolator()")
	print *, "fx   = ", fx, " (Neville rational)"
	print *, "f(x) = ", cotd_fn(x)
	print *, ""

	xi = [1, 2, 3, 4, 5, 6, 7, 8] * 1.d-3
	x = 0.5d0 * sum(xi([1, 2]))

	fx = neville_interpolator(xi, log_fn, x)
	call test(fx, log(x), 100 * tol, nfail, "neville_interpolator()")  ! big tol
	print *, "fx   = ", fx, " (Neville polynomial)"
	print *, "f(x) = ", log(x)
	print *, ""

	fx = neville_rational_interpolator(xi, log_fn, x)
	call test(fx, log(x), tol, nfail, "neville_rational_interpolator()")
	print *, "fx   = ", fx, " (Neville rational)"
	print *, "f(x) = ", log(x)
	print *, ""

end function chapter_2_example_2p2p4

!===============================================================================

integer function chapter_2_fft_1() result(nfail)

	double precision, parameter :: tol = 1.d-14
	double precision :: ts, freq, amp1, amp4, amp7, norm_diff
	double precision, allocatable :: t(:)

	double complex, allocatable :: x(:), xr(:), xx(:)

	integer :: i, sr, n

	write(*,*) CYAN // "Starting chapter_2_fft_1()" // COLOR_RESET

	nfail = 0

	! Reference for this fft example:
	!
	!     https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter24.03-Fast-Fourier-Transform.html
	!
	! The idea is to start with a signal composed of a few simple sine waves.
	! Then the expected fft of that signal is just a few spikes in the frequency
	! domain with zeroes (or near machine zero) for the other frequency
	! components

	sr = 256  ! sampling rate
	n = sr

	ts = 1.0d0 / sr  ! sampling interval
	t = [(i, i = 0, sr-1)] * ts
	!print *, "t = ", t

	freq = 1.d0 ; amp1 = 3.d0
	x = amp1 * sin(2 * PI * freq * t)
	!x = x + amp1 * cos(2 * PI * freq * t)
	!print *, "x = ", x

	freq = 4.d0 ; amp4 = 1.d0
	x = x + amp4 * sin(2 * PI * freq * t)
	!x = x + amp4 * cos(2 * PI * freq * t)

	freq = 7.d0 ; amp7 = 0.5d0
	x = x + amp7 * sin(2 * PI * freq * t)
	!x = x + amp7 * cos(2 * PI * freq * t)

	! Fortran has no complex edit descriptor, so print as pairs of reals
	print *, "x = "
	print "(2es18.6)", x(1: 5)
	!print "(2es18.6)", x

	xx = fft(x)

	!print *, "abs(xx) = ", abs(xx(1: 5))
	print *, "xx = "
	print "(2es18.6)", xx(1: 5)
	!print "(2es18.6)", xx

	print *, "abs(xx) = "
	print "(es18.6)", abs(xx(1: 5))

	! These indices of xx assume that the sample period is 1.  See
	! chapter_2_fft_2() for a better, more general example of matching fft
	! output with meaningful frequencies
	call test(abs(2 * xx(1+1)), amp1, tol, nfail, "fft() amp1")
	call test(abs(2 * xx(4+1)), amp4, tol, nfail, "fft() amp4")
	call test(abs(2 * xx(7+1)), amp7, tol, nfail, "fft() amp7")

	print *, "sum abs = ", sum(abs(xx))
	call test(sum(abs(xx)), amp1+amp4+amp7, tol, nfail, "fft() sum")
	print *, ""

	xr = ifft(xx)
	!xr = ifft(xx) / sr

	print *, "xr = "
	print "(2es18.6)", xr(1: 5)

	! Get the norm of the difference between the original signal x and the
	! round-trip fft() + ifft() signal xr
	norm_diff = norm2c(xr - x)
	print *, "norm_diff = ", norm_diff
	call test(norm_diff, 0.d0, 1.d-12, nfail, "fft() round trip")
	print *, ""

end function chapter_2_fft_1

!===============================================================================

integer function chapter_2_fft_2() result(nfail)

	character(len = *), parameter :: label = "chapter_2_fft_2"
	character(len = :), allocatable :: fout

	!double precision, parameter :: tol = 1.d-14
	double precision :: fs, dt, norm_diff, df, a1, a2, a3, a4, f1, f2, f3, f4
	double precision, allocatable :: t(:), r(:), freqs(:), ymag(:)

	double complex, allocatable :: x(:), y(:), xr(:)

	integer :: i, l, ny, fid, nrng
	integer, allocatable :: iy_max(:)

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	! Reference for this fft example:
	!
	!     https://www.mathworks.com/help/matlab/ref/fft.html

	!Fs = 1000;            % Sampling frequency                    
	!T = 1/Fs;             % Sampling period       
	!L = 1500;             % Length of signal
	!t = (0:L-1)*T;        % Time vector
	fs = 1000       ! sampling frequency
	dt = 1.d0 / fs  ! sampling period

	!! Number of signal sample points
	!l = 2048
	!l = 1500
	!!l = 256
	!l = 1024 * 16
	l = 15000

	! Time vector
	t = [(i, i = 0, l-1)] * dt
	!print *, "t = ", t(1: 5), " ... ", t(l-5: l)

	allocate(r(l))
	call random_seed(size = nrng)
	call random_seed(put = [(0, i = 1, nrng)])
	call random_number(r)

	! Tests expect these amplitudes to be in descending order
	a1 = 0.8  ! DC offset
	a2 = 0.5
	a3 = 0.35
	a4 = 0.3

	f1 = 0
	f2 = 5
	f3 = 3
	f4 = 8

	!! Create a signal
	!x = 0.8 + 0.7*sin(2*PI*3*t) + sin(2*PI*5*t)
	x = a1 &
		+ 2*a2 * sin(2*PI*f2*t) &
		+ 2*a3 * sin(2*PI*f3*t) &
		+ 2*a4 * sin(2*PI*f4*t)

	!! Add some random noise
	x = x + 1.5 * (2 * (r - 0.5d0))
	!x = x + 0.5 * (2 * (r - 0.5d0))

	! Write to file for plotting, e.g. with gnuplot
	fout = label // "_tx.out"
	open(newunit = fid, file = fout)
	write(fid, "(a)") "# t, x"
	do i = 1, l
		write(fid, "(3es18.6)") t(i), x(i)
	end do
	close(fid)

	y = fft(x)
	ny = size(y)

	! Frequency resolution
	df = fs / ny

	freqs = df * [(i, i = 0, ny-1)]

	! Normalize amplitudes to account for zero-padding
	ymag = abs(y) * size(y) / size(x)

	fout = label // "_fy.out"
	open(newunit = fid, file = fout)
	write(fid, "(a)") "# f, abs(y)"
	!do i = 1, ny
	do i = 1, min(200, ny)
		write(fid, "(2es18.6)") freqs(i), ymag(i)
	end do
	close(fid)

	! Find largest frequencies from ymag.  They should match the input signal
	! (including the dc offset)

	iy_max = select_max_n(ymag(1: ny/2), 5)

	!print *, "iy_max = ", iy_max
	print *, "most prominent amplitudes = "
	print "(es18.6)", [(ymag(iy_max(i)), i = 1, size(iy_max))]

	print *, "most prominent frequencies = "
	print "(es18.6)", [(freqs(iy_max(i)), i = 1, size(iy_max))]
	!print *, "df = ", df
	!print *, ""

	! These test assertions can be very sensitive to things like the sample size
	! `l`, due to artifacts like FFT leakage.  This doesn't mean that anything
	! is wrong, it's just that the tests have to be written carefully, or you
	! need to do more work to equivalence nearby frequencies which look like
	! multiple spikes when they're really just one leaky spike

	call test(ymag(iy_max(1)), a1, 0.1d0, nfail, "fft() amplitude 1")
	call test(ymag(iy_max(2)), a2, 0.1d0, nfail, "fft() amplitude 2")
	call test(ymag(iy_max(3)), a3, 0.1d0, nfail, "fft() amplitude 3")
	call test(ymag(iy_max(4)), a4, 0.1d0, nfail, "fft() amplitude 4")

	call test(freqs(iy_max(1)), f1, df, nfail, "fft() frequency 1")
	call test(freqs(iy_max(2)), f2, df, nfail, "fft() frequency 2")
	call test(freqs(iy_max(3)), f3, df, nfail, "fft() frequency 3")
	call test(freqs(iy_max(4)), f4, df, nfail, "fft() frequency 4")

	xr = ifft(y)
	fout = label // "_txr.out"
	open(newunit = fid, file = fout)
	write(fid, "(a)") "# t, xr"
	do i = 1, size(xr)
		write(fid, "(3es18.6)") dt * (i-1), xr(i)
	end do
	close(fid)

	! Truncate xr to undo the initial padding and verify round-trip correctness
	norm_diff = norm2c(x - xr(1: l))
	print *, "norm_diff = ", norm_diff
	call test(norm_diff, 0.d0, 1.d-12, nfail, "fft() 2 round trip")
	print *, ""

end function chapter_2_fft_2

!===============================================================================

integer function chapter_2_tridiag() result(nfail)
	! The general tridiagonal matrix algorithm is not part of chapter 2 per se,
	! but it's a nice way to solve for moments in cubic spline problems

	character(len = *), parameter :: label = "chapter_2_tridiag"

	double precision, allocatable :: a(:,:), bx(:), x(:)

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	! Small test
	allocate(a(3, 3))
	a(:,1) = [0, 2, 1]
	a(:,2) = [3, 5, 4]
	a(:,3) = [6, 7, 0]

	bx = [8, 9, 10]

	call tridiag_invmul(a, bx)
	print "(a,3es18.6)", "bx = ", bx
	call test(norm2(bx - [65, -122, 106]), 0.d0, 1.d-11, nfail, "tridiag_invmul() 3x3")

	deallocate(a)

	! Slightly bigger test
	!
	! The lower band is a(1,:) == [3, 6, 8, ...], the main diagonal is a(2,:) ==
	! [2, 5, 7, ...] (i.e. the middle column below), and the upper band is
	! a(3,:) == [1, 4, 5, ...].  In a traditional layout, this is the matrix:
	!
	!     [2, 1, 0, 0, 0, 0]
	!     [3, 5, 4, 0, 0, 0]
	!     [0, 6, 7, 5, 0, 0]
	!     [0, 0, 8, 9, 4, 0]
	!     [0, 0, 0, 7, 9, 5]
	!     [0, 0, 0, 0, 7, 8]
	!
	allocate(a(3, 6))
	a(:,1) = [0, 2, 1]
	a(:,2) = [3, 5, 4]
	a(:,3) = [6, 7, 5]
	a(:,4) = [8, 9, 4]
	a(:,5) = [7, 9, 5]
	a(:,6) = [7, 8, 0]

	bx = [8, 9, 10, 11, 12, 13]

	call tridiag_invmul(a, bx)
	print "(a,6es18.6)", "bx = ", bx
	x = &
	[   &
		 4.3366500829187418d+00, &
		-6.7330016583748276d-01, &
		-1.6086235489220257d-01, &
		 3.0331674958540629d+00, &
		-3.7529021558872304d+00, &
		 4.9087893864013266d+00  &
	]
	call test(norm2(bx - x), 0.d0, 1.d-11, nfail, "tridiag_invmul() 6x6")

	print *, ""

end function chapter_2_tridiag

!===============================================================================

integer function chapter_4_lu() result(nfail)

	character(len = *), parameter :: label = "chapter_4_lu"

	double precision, allocatable :: a(:,:), bx(:), x(:)
	integer :: i, n, nrng, irep

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Test a fixed case
	allocate(a(5, 5))

	a(:,1) = [-4.130000d+00, 6.000000d+00, 1.100000d+01, 1.600000d+01, 2.100000d+01]
	a(:,2) = [ 2.000000d+00, 7.960000d+00, 1.200000d+01, 1.700000d+01, 2.200000d+01]
	a(:,3) = [ 3.000000d+00, 8.000000d+00, 1.280000d+00, 1.800000d+01, 2.300000d+01]
	a(:,4) = [ 4.000000d+00, 9.000000d+00, 1.400000d+01, 3.400000d+00, 2.400000d+01]
	a(:,5) = [ 5.000000d+00, 1.000000d+01, 1.500000d+01, 2.000000d+01, 5.600000d+00]

	! TODO: fill the above matrix in with some zeros to compare tridiag/banded solvers

	x = [+4.0000d+00, +2.0000d+00, -1.0000d+00, +7.0000d+00, +1.8000d+01]

	! Print this to check that the solution is correct in the transpose sense
	print *, "a * x = ", matmul(a, x)

	bx = [102.48d0, 274.92d0, 434.72d0, 463.8d0, 373.8d0]

	print *, "a = "
	print "(5es18.6)", a

	call lu_invmul(a, bx)
	print "(a,5es18.6)", "bx = ", bx
	call test(norm2(bx - x), 0.d0, 1.d-11, nfail, "lu_invmul() 5x5")

	deallocate(a, x, bx)

	!********
	! Do some fuzz testing with random data.  Chances are almost 0 for randomly
	! creating a singular matrix
	!
	! TODO: add fuzz testing for some other problems, e.g. tridiagonal and
	! banded solvers.  We will need special matmul() implementations to verify
	! results

	call random_seed(size = nrng)
	call random_seed(put = [(0, i = 1, nrng)])

	do n = 1, 90

		! LU decomposition is still fast for n >> 90, but maybe not fast enough
		! for multiple reps in a unit test

		allocate(a(n, n), x(n))

		if (mod(n, 10) == 0) then
			print *, "Testing lu_invmul() with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 10
			call random_number(a)  ! random matrix
			call random_number(x)  ! random expected answer

			bx = matmul(a, x)      ! calculate rhs

			!print *, "a = "
			!print "(6es15.5)", a
			!print "(a,6es15.5)", "x  = ", x
			!print "(a,6es15.5)", "b  = ", bx

			call lu_invmul(a, bx)
			!print "(a,6es15.5)", "bx = ", bx

			! Rounding error grows with `n`
			call test(norm2(bx - x), 0.d0, 1.d-9 * n, nfail, "lu_invmul() fuzz n x n")

		end do

		deallocate(a, x, bx)
	end do

	!********
	print *, ""

end function chapter_4_lu

!===============================================================================

integer function chapter_2_tridiag_corner() result(nfail)
	! The general tridiagonal matrix algorithm is not part of chapter 2 per se,
	! but it's a nice way to solve for moments in cubic spline problems

	character(len = *), parameter :: label = "chapter_2_tridiag_corner"

	double precision, allocatable :: a(:,:), bx(:), x(:), ad(:,:)

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET
	nfail = 0

	! Slightly bigger test
	!
	! The lower band is a(1,:) == [3, 6, 8, ...], the main diagonal is a(2,:) ==
	! [2, 5, 7, ...] (i.e. the middle column below), and the upper band is
	! a(3,:) == [1, 4, 5, ...].  In a traditional layout, this is the matrix:
	!
	!     [ 2,  1,  0,  0,  0, 11]
	!     [ 3,  5,  4,  0,  0,  0]
	!     [ 0,  6,  7,  5,  0,  0]
	!     [ 0,  0,  8,  9,  4,  0]
	!     [ 0,  0,  0,  7,  9,  5]
	!     [12,  0,  0,  0,  7,  8]
	!

	! Dense LU version for comparison
	allocate(ad(6, 6))
	ad(:,1) = [ 2,  1,  0,  0,  0, 11]
	ad(:,2) = [ 3,  5,  4,  0,  0,  0]
	ad(:,3) = [ 0,  6,  7,  5,  0,  0]
	ad(:,4) = [ 0,  0,  8,  9,  4,  0]
	ad(:,5) = [ 0,  0,  0,  7,  9,  5]
	ad(:,6) = [12,  0,  0,  0,  7,  8]
	ad = transpose(ad)  ! TODO: i think all my tridiag algos are transposed
	print *, "ad = "
	print "(6es15.5)", transpose(ad)

	! Tridiagonal packing
	allocate(a(3, 6))
	a(:,1) = [11.d0,  2.d0,  1.d0]
	a(:,2) = [ 3.d0,  5.d0,  4.d0]
	a(:,3) = [ 6.d0,  7.d0,  5.d0]
	a(:,4) = [ 8.d0,  9.d0,  4.d0]
	a(:,5) = [ 7.d0,  9.d0,  5.d0]
	a(:,6) = [ 7.d0,  8.d0, 12.d0]

	bx = [8.d0, 9.d0, 10.d0, 11.d0, 12.d0, 13.d0]
	call lu_invmul(ad, bx)
	print "(a,6es18.6)", "bx (dense) = ", bx

	bx = [8.d0, 9.d0, 10.d0, 11.d0, 12.d0, 13.d0]
	call tridiag_corner_invmul(a, bx)
	print "(a,6es18.6)", "bx         = ", bx
	x = &
	[   &
		-0.67986230636832845d0, &
		 0.54056519116760138d0, &
		 2.0841902408167448d0, &
		-1.5665445665445641d0, &
		 2.1063447930917798d0, &
		 0.80174176559718680d0 &
	]
	call test(norm2(bx - x), 0.d0, 1.d-11, nfail, "tridiag_corner_invmul() 6x6")
	!print *, "expected = ", x

	print *, ""

end function chapter_2_tridiag_corner

!===============================================================================

integer function chapter_2_banded() result(nfail)
	! Again, not part of chapter 2, but I've gone off on a slight tangent after
	! the tridiagonal matrix algorithm

	character(len = *), parameter :: label = "chapter_2_banded"

	double precision, allocatable :: a(:,:), bx(:), x(:)
	integer :: nl, nu, n

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	nl = 2  ! number of lower bands
	nu = 2  ! number of upper bands
	n = 11  ! size of x
	allocate(a(nl+nu+1, n))
	a(:,  1) = [0, 0, 5, 2, 1]
	a(:,  2) = [0, 3, 6, 2, 1]
	a(:,  3) = [2, 3, 7, 2, 1]
	a(:,  4) = [2, 3, 8, 2, 1]
	a(:,  5) = [2, 3, 9, 2, 1]
	a(:,  6) = [2, 3, 8, 2, 1]
	a(:,  7) = [2, 3, 7, 2, 1]
	a(:,  8) = [2, 3, 6, 2, 1]
	a(:,  9) = [2, 3, 5, 2, 1]
	a(:, 10) = [2, 3, 4, 2, 0]
	a(:, 11) = [2, 3, 3, 0, 0]

	bx = [1, 2, 3, 4, 2, 3, 4, 5, 6, 7, 8]

	x = &
	[   &
		 9.3517985708705731d-02, &
		 1.4687004667042000d-01, &
		 2.3866997811563129d-01, &
		 3.6088580662010022d-01, &
		-2.0107571478290644d-02, &
		 1.4337866230804550d-01, &
		 3.3421344259696145d-01, &
		 5.2309491753638426d-01, &
		 2.2439522278094537d-01, &
		 1.2324239681282793d-01, &
		 2.3938274546665417d+00  &
	]

	call banded_invmul(a, bx, nl, nu)
	print "(a)", "bx = "
	print "(es28.16)", bx
	!print "(a)", "bx, x = "
	!do i = 1, n
	!	print "(2es18.6)", bx(i), x(i)
	!end do
	call test(norm2(bx - x), 0.d0, 1.d-11, nfail, "banded_invmul() (2+2)x11")

	deallocate(a)

	!********

	! Test an asymmetric band pattern

	nl = 3  ! number of lower bands
	nu = 2  ! number of upper bands
	n = 11  ! size of x
	allocate(a(nl+nu+1, n))
	a(:,  1) = [0, 0, 0, 5, 2, 1]
	a(:,  2) = [0, 0, 3, 6, 2, 1]
	a(:,  3) = [0, 2, 3, 7, 2, 1]
	a(:,  4) = [1, 2, 3, 8, 2, 1]
	a(:,  5) = [2, 2, 3, 9, 2, 1]
	a(:,  6) = [1, 2, 3, 8, 2, 1]
	a(:,  7) = [1, 2, 3, 7, 2, 1]
	a(:,  8) = [2, 2, 3, 6, 2, 1]
	a(:,  9) = [1, 2, 3, 5, 2, 1]
	a(:, 10) = [2, 2, 3, 4, 2, 0]
	a(:, 11) = [1, 2, 3, 3, 0, 0]

	bx = [1, 2, 3, 4, 2, 3, 4, 5, 6, 7, 8]

	x = &
	[   &
		 9.2539838515047296d-02, &
		 1.4666377015277052d-01, &
		 2.4397326711922243d-01, &
		 3.5445132929979012d-01, &
		-4.1786515922543226d-02, &
		 1.3017521726850928d-01, &
		 2.7110014632251389d-01, &
		 5.8888159125572470d-01, &
		 2.6313184397072215d-01, &
		-4.9631077135265822d-02, &
		 2.3445826507362098d+00  &
	]

	call banded_invmul(a, bx, nl, nu)
	print "(a)", "bx2 = "
	print "(es28.16)", bx
	call test(norm2(bx - x), 0.d0, 1.d-11, nfail, "banded_invmul() (3+2)x11")

	deallocate(a)

	!********

	! Test an asymmetric band pattern with more upper bands than lower
	! (transpose of the last test)

	nl = 2  ! number of lower bands
	nu = 3  ! number of upper bands
	n = 11  ! size of x
	allocate(a(nl+nu+1, n))
	a(:,  1) = [0, 0, 5, 3, 2, 1]
	a(:,  2) = [0, 2, 6, 3, 2, 2]
	a(:,  3) = [1, 2, 7, 3, 2, 1]
	a(:,  4) = [1, 2, 8, 3, 2, 1]
	a(:,  5) = [1, 2, 9, 3, 2, 2]
	a(:,  6) = [1, 2, 8, 3, 2, 1]
	a(:,  7) = [1, 2, 7, 3, 2, 2]
	a(:,  8) = [1, 2, 6, 3, 2, 1]
	a(:,  9) = [1, 2, 5, 3, 2, 0]
	a(:, 10) = [1, 2, 4, 3, 0, 0]
	a(:, 11) = [1, 2, 3, 0, 0, 0]

	bx = [1, 2, 3, 4, 2, 3, 4, 5, 6, 7, 8]

	x = &
	[   &
		-5.7197635714385166d-02, &
		 1.2962566645714274d-01, &
		 2.6474979163633455d-01, &
		 3.6761159592782849d-01, &
		-9.6415647039373428d-02, &
		 3.4694267641018961d-02, &
		 6.1984038868364255d-01, &
		 2.6200212978601356d-01, &
		 1.6414013139981204d-01, &
		-7.1307113059291349d-01, &
		 3.0873340432620049d+00  &
	]

	call banded_invmul(a, bx, nl, nu)
	print "(a)", "bx3 = "
	print "(es28.16)", bx
	call test(norm2(bx - x), 0.d0, 1.d-11, nfail, "banded_invmul() (2+3)x11")

	deallocate(a)

	print *, ""
end function chapter_2_banded

!===============================================================================

integer function chapter_2_cubic_splines() result(nfail)

	character(len = *), parameter :: label = "chapter_2_cubic_splines"

	double precision :: diff, dy_start, dy_end, phase
	double precision, allocatable :: xi(:), yi(:), x(:), fx(:)

	integer :: i, fid

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	! Set support points `xi`.  They don't need to be evenly spaced
	!
	! For the no-curve BC, carefully choose a fn to interpolate which happens to
	! have 0 2nd derivative at end-points of interpolated interval, i.e. sin
	! from 0 to PI

	!xi = [0.d0, 1.2d0, 2.4d0, PI]
	!xi = [0.d0, PI/3, 2*PI/3, PI]
	xi = [0.d0, PI/4, PI/2, 3*PI/4, PI]
	!xi = [0.d0, PI/2, PI]

	yi = sin(xi)
	!print *, "xi = ", xi
	!print *, "yi = ", yi

	x = 3.1415d0 * [(i, i = 0, 100)] / 100.d0
	fx = spline_no_curve(xi, yi, x)

	diff = sum(abs(fx - sin(x)))
	print *, "diff = ", diff
	call test(diff, 4.1936105024329519d-002, 1.d-11, nfail, "spline_no_curve 1")

	open(file = "plot-spline-1.txt", newunit = fid)
	write(fid, *) "# x, f(x), sin(x)"
	!write(fid, "(2es18.6)") x(i), fx(i)
	write(fid, "(3es18.6)") [(x(i), fx(i), sin(x(i)), i = 1, size(x))]
	close(fid)

	!********
	! Different number of support/control points

	xi = [0.d0, PI/5, 2*PI/5, 3*PI/5, 4*PI/5, PI]
	yi = sin(xi)
	x = 3.1415d0 * [(i, i = 0, 100)] / 100.d0
	fx = spline_no_curve(xi, yi, x)

	diff = sum(abs(fx - sin(x)))
	print *, "diff = ", diff
	call test(diff, 1.5928426416190206d-002, 1.d-11, nfail, "spline_no_curve 2")

	open(file = "plot-spline-2.txt", newunit = fid)
	write(fid, *) "# x, f(x), sin(x)"
	write(fid, "(3es18.6)") [(x(i), fx(i), sin(x(i)), i = 1, size(x))]
	close(fid)

	!********
	! Unevenly spaced control points

	xi = [0.d0, PI/5, 2*PI/5, 4*PI/5, PI]
	yi = sin(xi)
	x = 3.1415d0 * [(i, i = 0, 100)] / 100.d0
	fx = spline_no_curve(xi, yi, x)

	diff = sum(abs(fx - sin(x)))
	print *, "diff = ", diff
	call test(diff, 0.29268525551874491d0, 1.d-11, nfail, "spline_no_curve 3")

	open(file = "plot-spline-3.txt", newunit = fid)
	write(fid, *) "# x, f(x), sin(x)"
	write(fid, "(3es18.6)") [(x(i), fx(i), sin(x(i)), i = 1, size(x))]
	close(fid)

	!********
	! Prescribed 1st derivatives
	xi = [0.d0, PI/4, PI/2]
	yi = sin(xi)
	x = 0.5d0 * 3.1415d0 * [(i, i = 0, 100)] / 100.d0
	dy_start = 1.d0
	dy_end   = 0.d0
	fx = spline_prescribed(xi, yi, x, dy_start, dy_end)

	diff = sum(abs(fx - sin(x)))
	print *, "diff = ", diff
	call test(diff, 3.5556039482586774d-002, 1.d-11, nfail, "spline_prescribed 4")

	open(file = "plot-spline-4.txt", newunit = fid)
	write(fid, *) "# x, f(x), sin(x)"
	write(fid, "(3es18.6)") [(x(i), fx(i), sin(x(i)), i = 1, size(x))]
	close(fid)

	!********
	! Periodic spline with matching derivatives
	xi = [0.d0, PI/2, PI, 3*PI/2, 2*PI]
	yi = cos(xi)
	x = 2 * 3.1415d0 * [(i, i = 0, 100)] / 100.d0
	fx = spline_periodic(xi, yi, x)

	diff = sum(abs(fx - cos(x)))
	print *, "diff = ", diff
	call test(diff, 1.161d0, 1d-3, nfail, "spline_periodic 5")

	open(file = "plot-spline-5.txt", newunit = fid)
	write(fid, *) "# x, f(x), cos(x)"
	write(fid, "(3es18.6)") [(x(i), fx(i), cos(x(i)), i = 1, size(x))]
	close(fid)

	!********
	! Periodic with more support points
	xi = [0.d0, PI/4, PI/2, 3*PI/4, PI, 5*PI/4, 3*PI/2, 7*PI/4, 2*PI]
	yi = cos(xi)
	x = 2 * 3.1415d0 * [(i, i = 0, 100)] / 100.d0
	fx = spline_periodic(xi, yi, x)

	diff = sum(abs(fx - cos(x)))
	print *, "diff = ", diff
	call test(diff, 4.196d-002, 1d-4, nfail, "spline_periodic 6")

	open(file = "plot-spline-6.txt", newunit = fid)
	write(fid, *) "# x, f(x), cos(x)"
	write(fid, "(3es18.6)") [(x(i), fx(i), cos(x(i)), i = 1, size(x))]
	close(fid)

	!********
	! Phase-shifted cosine
	xi = [0.d0, PI/4, PI/2, 3*PI/4, PI, 5*PI/4, 3*PI/2, 7*PI/4, 2*PI]
	phase = -PI/3
	yi = cos(xi + phase)
	x = 2 * 3.1415d0 * [(i, i = 0, 100)] / 100.d0
	fx = spline_periodic(xi, yi, x)

	diff = sum(abs(fx - cos(x + phase)))
	print *, "diff = ", diff
	call test(diff, 4.105d-2, 1d-4, nfail, "spline_periodic 7")

	open(file = "plot-spline-7.txt", newunit = fid)
	write(fid, *) "# x, f(x), cos(x + phase)"
	write(fid, "(3es18.6)") [(x(i), fx(i), cos(x(i) + phase), i = 1, size(x))]
	close(fid)

	!********
	! Phase-shifted cosine and uneven support points
	xi = [0.d0, PI/4, PI/2, PI, 5*PI/4, 3*PI/2, 7*PI/4, 2*PI]
	phase = -PI/3
	yi = cos(xi + phase)
	x = 2 * 3.1415d0 * [(i, i = 0, 100)] / 100.d0
	fx = spline_periodic(xi, yi, x)

	diff = sum(abs(fx - cos(x + phase)))
	print *, "diff = ", diff
	call test(diff, 0.1374d0, 1.d-4, nfail, "spline_periodic 8")

	open(file = "plot-spline-8.txt", newunit = fid)
	write(fid, *) "# x, f(x), cos(x + phase)"
	write(fid, "(3es18.6)") [(x(i), fx(i), cos(x(i) + phase), i = 1, size(x))]
	close(fid)

	!********
	print *, ""

end function chapter_2_cubic_splines

!===============================================================================

integer function chapter_2_bezier_splines() result(nfail)

	character(len = *), parameter :: label = "chapter_2_bezier_splines"

	integer :: i, nc, fid
	double precision, allocatable :: xyc(:,:), t(:), xy(:,:)

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Make a plot

	nc = 4  ! number of control points
	allocate(xyc(2, nc))
	xyc(1,:) = [0, 0, 1, 1]  ! x control coordinates
	xyc(2,:) = [0, 1, 1, 0]  ! y control coordinates (xyc can be n-dimensional)

	t = [(i, i = 0, 100)] / 100.d0  ! interpolation parameter in range [0, 1]

	xy = bezier_curve(xyc, t)

	open(file = "plot-bezier-1.txt", newunit = fid)
	write(fid, *) "# x, y"
	write(fid, "(2es18.6)") [(xy(1,i), xy(2,i), i = 1, size(xy, 2))]
	close(fid)

	!********
	! Test a value from the curve above
	t = [0.5d0]
	xy = bezier_curve(xyc, t)
	!print *, "xy = ", xy

	call test(xy(1,1), 0.50d0, 1.d-12, nfail, "bezier_curve 1")
	call test(xy(2,1), 0.75d0, 1.d-12, nfail, "bezier_curve 2")

	deallocate(xyc)
	!********
	! Plot a cardinal spline

	nc = 4  ! number of control points
	allocate(xyc(2, nc))
	xyc(1,:) = [0, 0, 1, 1]  ! x control coordinates
	xyc(2,:) = [0, 1, 1, 0]  ! y control coordinates (xyc can be n-dimensional)

	t = (nc-1) * [(i, i = 0, 100)] / 100.d0  ! interpolation parameter in range [0, 1]

	xy = cardinal_spline(xyc, t, 1.d0)
	!xy = cardinal_spline(xyc, t, 1.5d0)

	open(file = "plot-bezier-2.txt", newunit = fid)
	write(fid, *) "# x, y"
	write(fid, "(2es18.6)") [(xy(1,i), xy(2,i), i = 1, size(xy, 2))]
	close(fid)

	!********
	! Test some values from the curve above
	t = [0.5d0, 1.5d0, 2.5d0]
	xy = cardinal_spline(xyc, t, 1.d0)
	!print *, "xy = ", xy

	call test(xy(1,1), -0.125d0, 1.d-12, nfail, "cardinal_spline x1")
	call test(xy(2,1),  0.625d0, 1.d-12, nfail, "cardinal_spline y1")

	call test(xy(1,2),  0.500d0, 1.d-12, nfail, "cardinal_spline x2")
	call test(xy(2,2),  1.250d0, 1.d-12, nfail, "cardinal_spline y2")

	call test(xy(1,3),  1.125d0, 1.d-12, nfail, "cardinal_spline x3")
	call test(xy(2,3),  0.625d0, 1.d-12, nfail, "cardinal_spline y3")

	deallocate(xyc)

	!********
	! Cardinal spline with more points

	nc = 7  ! number of control points
	allocate(xyc(2, nc))
	xyc(1,:) = [0, 0, 1, 1,  1,  2, 2]  ! x control coordinates
	xyc(2,:) = [0, 1, 1, 0, -1, -1, 0]  ! y control coordinates (xyc can be n-dimensional)

	t = (nc-1) * [(i, i = 0, 100)] / 100.d0  ! interpolation parameter in range [0, 1]
	xy = cardinal_spline(xyc, t, 1.d0)
	!xy = cardinal_spline(xyc, t, 0.5d0)
	!print *, "xy = ", xy

	open(file = "plot-bezier-3.txt", newunit = fid)
	write(fid, *) "# x, y"
	write(fid, "(2es18.6)") [(xy(1,i), xy(2,i), i = 1, size(xy, 2))]
	close(fid)

	!********
	! Test some values from the curve above

	! Tension = 1
	t = [2.5d0]
	xy = cardinal_spline(xyc, t, 1.d0)
	!print *, "xy = ", xy

	call test(xy(1,1),  1.125d0, 1.d-12, nfail, "cardinal_spline x3 t 1")
	call test(xy(2,1),  0.625d0, 1.d-12, nfail, "cardinal_spline y3 t 1")

	! Tension = 0.001
	t = [2.5d0]
	xy = cardinal_spline(xyc, t, 0.001d0)
	!print *, "xy = ", xy

	call test(xy(1,1),  1.000125d0, 1.d-12, nfail, "cardinal_spline x3 t 0.001")
	call test(xy(2,1),  0.500125d0, 1.d-12, nfail, "cardinal_spline y3 t 0.001")

	! Tension = 0.5
	t = [2.5d0]
	xy = cardinal_spline(xyc, t, 0.5d0)
	!print *, "xy = ", xy

	call test(xy(1,1),  1.0625d0, 1.d-12, nfail, "cardinal_spline x3 t 0.5")
	call test(xy(2,1),  0.5625d0, 1.d-12, nfail, "cardinal_spline y3 t 0.5")

	! Tension = 1.5
	t = [2.5d0]
	xy = cardinal_spline(xyc, t, 1.5d0)
	!print *, "xy = ", xy

	call test(xy(1,1),  1.1875d0, 1.d-12, nfail, "cardinal_spline x3 t 1.5")
	call test(xy(2,1),  0.6875d0, 1.d-12, nfail, "cardinal_spline y3 t 1.5")

	deallocate(xyc)

	!********
	print *, ""

end function chapter_2_bezier_splines

!===============================================================================

integer function chapter_3_newton_cotes() result(nfail)

	! The Newton-Cotes formulas are generalizations of integration rules like
	! the trapezoid rule, and higher-order methods like Simpson's (1/3) rule,
	! Simpson's 3/8 rule, etc.

	character(len = *), parameter :: label = "chapter_3_newton_cotes"

	double precision :: area
	double precision, allocatable :: xi(:), yi(:)

	integer :: i

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0
	!********

	area = simpson_13_integrator(sin_fn, 0.d0, PI, PI)
	print *, "area = ", area
	call test(area, 2.094395d0, 1.d-6, nfail, "simpson_13_integrator 1")

	area = simpson_13_integrator(sin_fn, 0.d0, PI, 0.1d0)
	print *, "area = ", area
	call test(area, 2.0000000732694589d0, 1.d-12, nfail, "simpson_13_integrator 2")

	area = simpson_38_integrator(sin_fn, 0.d0, PI, PI)
	print *, "area = ", area
	call test(area, 2.04052428d0, 1.d-6, nfail, "simpson_38_integrator 3")

	area = simpson_38_integrator(sin_fn, 0.d0, PI, 0.1d0)
	print *, "area = ", area
	call test(area, 2.0000000325630984d0, 1.d-12, nfail, "simpson_38_integrator 4")

	!********

	area = milne_integrator(sin_fn, 0.d0, PI, PI)
	print *, "area = ", area
	call test(area, 1.99857d0, 1.d-6, nfail, "milne_integrator 5")

	area = milne_integrator(sin_fn, 0.d0, PI, 0.1d0)
	print *, "area = ", area
	call test(area, 1.9999999999988802d0, 1.d-12, nfail, "milne_integrator 6")

	!********

	area = weddle_integrator(sin_fn, 0.d0, PI, PI)
	print *, "area = ", area
	call test(area, 2.0000178136366551d0, 1.d-12, nfail, "weddle_integrator 5")

	area = weddle_integrator(sin_fn, 0.d0, PI, 0.1d0)
	print *, "area = ", area
	call test(area, 1.9999999999999998d0, 1.d-12, nfail, "weddle_integrator 6")

	area = weddle_integrator(sinx_x, 0.d0, 1.d0, 0.1d0)
	print *, "area = ", area
	call test(area, 0.94608307036718320d0, 1.d-12, nfail, "weddle_integrator 6.1")

	!********
	! Value array integrator(s)

	!print *, "value integrators:"

	xi = PI * [(i, i = 0, 2)] / 2.d0
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0943951023931953d0, 1.d-12, nfail, "simpson_integrator_vals 7")

	xi = PI * [(i, i = 0, 3)] / 3.d0
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0405242847634950d0, 1.d-12, nfail, "simpson_integrator_vals 8")

	xi = PI * [(i, i = 0, 4)] / 4.d0
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0045597549844207d0, 1.d-12, nfail, "simpson_integrator_vals 9")

	xi = PI * [(i, i = 0, 5)] / 5.d0
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0034411937297119d0, 1.d-12, nfail, "simpson_integrator_vals 10")

	xi = PI * [(i, i = 0, 6)] / 6.d0
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0008631896735367d0, 1.d-12, nfail, "simpson_integrator_vals 11")

	xi = PI * [(i, i = 0, 7)] / 7.d0
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0006963918546892d0, 1.d-12, nfail, "simpson_integrator_vals 12")

	xi = PI * [(i, i = 0, 30)] / 30.d0
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0000013379479498d0, 1.d-12, nfail, "simpson_integrator_vals 13")

	xi = PI * [(i, i = 0, 31)] / 31.d0
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0000012070944040d0, 1.d-12, nfail, "simpson_integrator_vals 14")

	xi = PI * [(i, i = 0, 32)] / 32.d0
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0000010333694127d0, 1.d-12, nfail, "simpson_integrator_vals 15")

	xi = PI * [(i, i = 0, 33)] / 33.d0
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0000009368045015d0, 1.d-12, nfail, "simpson_integrator_vals 16")

	!********
	print *, ""

end function chapter_3_newton_cotes

!===============================================================================

integer function chapter_3_romberg() result(nfail)

	character(len = *), parameter :: label = "chapter_3_romberg"

	double precision :: area
	double precision, allocatable :: areas(:), expect(:)

	integer :: i, n

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0
	!********

	area = romberg_integrator_fixed(sin_fn, 0.d0, PI, 4)
	print *, "area = ", area
	call test(area, 1.9999999945872906d0, 1.d-12, nfail, "romberg_integrator_fixed 1")

	!********
	print *, "integrating sin(x)/x:"

	n = 6
	allocate(areas(n))
	do i = 1, n
		areas(i) = romberg_integrator_fixed(sinx_x, 0.d0, 1.d0, i)
	end do
	print *, "areas = "
	print "(es24.14)", areas

	expect = [ &
		9.46145882273587d-01, &
		9.46083004063674d-01, &
		9.46083070387223d-01, &
		9.46083070367181d-01, &
		9.46083070367183d-01, &
		9.46083070367183d-01  &
	]

	call tests(areas, expect, 1.d-12, nfail, "romberg_integrator_fixed 2")

	!********

	area = romberg_integrator(sin_fn, 0.d0, PI, 1.d-3)
	print *, "area = ", area
	call test(area, 1.9999999945872902d0, 1.d-12, nfail, "romberg_integrator 3")

	area = romberg_integrator(sin_fn, 0.d0, PI, 1.d-6)
	print *, "area = ", area
	call test(area, 2.0000000000013212d0, 1.d-12, nfail, "romberg_integrator 4")

	area = romberg_integrator(sin_fn, 0.d0, PI, 1.d-9)
	print *, "area = ", area
	call test(area, 1.9999999999999991d0, 1.d-12, nfail, "romberg_integrator 5")

	area = romberg_integrator(sinx_x, 0.d0, 1.d0, 1.d-9)
	print *, "area = ", area
	call test(area, 9.46083070367183d-01, 1.d-12, nfail, "romberg_integrator 5.1")

	area = romberg_integrator(fn_example_pg171, 0.d0, PI/2, 1.d-12)
	print *, "area = ", area
	call test(area, 1.d0, 1.d-11, nfail, "romberg_integrator 5.2")

	!********
	! These tests will throw convergence warnings
	print *, GREEN // "Note" // COLOR_RESET // ": the following two test warnings are expected"

	area = romberg_integrator(sin_fn, 0.d0, PI, 0.d0)
	print *, "area = ", area
	call test(area, 2.d0, 1.d-14, nfail, "romberg_integrator 6")

	area = romberg_integrator(sin_fn, 0.d0, PI, tol = 0.d0, max_levels = 4)
	print *, "area = ", area
	call test(area, 1.9999999945872902d0, 1.d-14, nfail, "romberg_integrator 7")

	!********
	print *, ""

end function chapter_3_romberg

!===============================================================================

integer function chapter_3_gauss() result(nfail)

	character(len = *), parameter :: label = "chapter_3_gauss"

	double precision :: area

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0
	!********

	area = gauss2_single(sin_fn, 0.d0, PI)
	print *, "area = ", area
	call test(area, 1.9358195746511373d0, 1.d-12, nfail, "gauss2_single 1")

	area = gauss3_single(sin_fn, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0013889136077436d0, 1.d-12, nfail, "gauss3_single 2")

	area = gauss4_single(sin_fn, 0.d0, PI)
	print *, "area = ", area
	call test(area, 1.9999842284577216d0, 1.d-12, nfail, "gauss4_single 3")

	area = gauss5_single(sin_fn, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0000001102844722d0, 1.d-12, nfail, "gauss5_single 4")

	area = gauss7_single(sin_fn, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0000000000017906d0, 1.d-12, nfail, "gauss7_single 5")

	area = kronrod15_single(sin_fn, 0.d0, PI)
	print *, "area = ", area
	call test(area, 1.9999999999999942d0, 1.d-12, nfail, "kronrod15_single 6")

	!********
	print *, ""

end function chapter_3_gauss

!===============================================================================

integer function chapter_3_adaptive() result(nfail)

	character(len = *), parameter :: label = "chapter_3_adaptive"

	double precision :: area

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0
	!********

	area = simpson_adaptive_integrator(sin_fn, 0.d0, PI, 1.d-10)
	print *, "area = ", area
	call test(area, 1.9999999999999993d0, 1.d-12, nfail, "simpson_adaptive_integrator 1")

	!********
	print *, ""

end function chapter_3_adaptive

!===============================================================================

end module numa__exercises

