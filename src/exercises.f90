
module exercises_m

	use numerical_analysis_m
	use utils_m

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

	fx = lagrange_interpolater(xi, log_fn, x)
	!nfail = nfail + assert(abs(fx - log_fn(x)) < tol)
	call test(fx, log_fn(x), tol, nfail, "lagrange_interpolater()")
	print *, "fx   = ", fx, " (Lagrange)"
	print *, "f(x) = ", log_fn(x)
	print *, ""

	!fx = lagrange_interpolater(xi, exp_fn, x)  ! not from the text book, but let's see what happens
	!print *, "fx   = ", fx
	!print *, "f(x) = ", exp_fn(x)
	!print *, ""

	fx = neville_interpolater(xi, log_fn, x)
	!fx = neville_interpolater_vals(xi, log(xi), x)
	call test(fx, log_fn(x), tol, nfail, "neville_interpolater()")
	print *, "fx   = ", fx, " (Neville)"
	print *, "f(x) = ", log_fn(x)
	print *, ""

	xs = [11.1d0, 11.2d0]
	!xs = [11.1d0, 11.2d0, 11.3d0, 11.4d0]
	fxs = newton_interpolater(xi, log_fn, xs)
	!fxs = newton_interpolater_vals(xi, log(xi), xs)
	call tests(fxs, log(xs), tol, nfail, "newton_interpolater()")
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

	fx = neville_interpolater(xi, cotd_fn, x)
	call test(fx, 22.6352d0, tol, nfail, "neville_interpolater()")
	print *, "fx   = ", fx, " (Neville polynomial)"
	print *, "f(x) = ", cotd_fn(x)
	print *, ""

	fx = neville_rational_interpolater(xi, cotd_fn, x)
	call test(fx, cotd_fn(x), tol, nfail, "neville_rational_interpolater()")
	print *, "fx   = ", fx, " (Neville rational)"
	print *, "f(x) = ", cotd_fn(x)
	print *, ""

	xi = [1, 2, 3, 4, 5, 6, 7, 8] * 1.d-3
	x = 0.5d0 * sum(xi([1, 2]))

	fx = neville_interpolater(xi, log_fn, x)
	call test(fx, log(x), 100 * tol, nfail, "neville_interpolater()")  ! big tol
	print *, "fx   = ", fx, " (Neville polynomial)"
	print *, "f(x) = ", log(x)
	print *, ""

	fx = neville_rational_interpolater(xi, log_fn, x)
	call test(fx, log(x), tol, nfail, "neville_rational_interpolater()")
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
	n = sr    ! TODO: this needs be be disentangled for better testing

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
	print "(2es18.6)", x(1: 10)
	!print "(2es18.6)", x

	xx = fft(x)

	!print *, "abs(xx) = ", abs(xx(1: 10))
	print *, "xx = "
	print "(2es18.6)", xx(1: 10)
	!print "(2es18.6)", xx

	print *, "abs(xx) = "
	print "(es18.6)", abs(xx(1: 10))

	! TODO: generalize.  These indices of xx assume that the sample period is 1
	call test(abs(2 * xx(1+1)), amp1, tol, nfail, "fft() amp1")
	call test(abs(2 * xx(4+1)), amp4, tol, nfail, "fft() amp4")
	call test(abs(2 * xx(7+1)), amp7, tol, nfail, "fft() amp7")

	print *, "sum abs = ", sum(abs(xx))
	call test(sum(abs(xx)), amp1+amp4+amp7, tol, nfail, "fft() sum")
	print *, ""

	xr = ifft(xx)
	!xr = ifft(xx) / sr

	print *, "xr = "
	print "(2es18.6)", xr(1: 10)
	!print *, "xr = ", xr
	!print "(a,20es18.6)", "xr = ", xr(1: 10)
	print *, ""

	! Get the norm of the difference between the original signal x and the
	! round-trip fft() + ifft() signal xr

	!print *, "norm_diff = ", norm2(xr - x)
	!print *, "norm_diff = ", nrm2(xr - x)
	!print *, "norm_diff = ", nrm2(n, xr - x, 1)  ! undefined

	norm_diff = norm2c(xr - x)
	print *, "norm_diff = ", norm_diff
	call test(norm_diff, 0.d0, 1.d-12, nfail, "fft() round trip")

end function chapter_2_fft_1

!===============================================================================

integer function chapter_2_fft_2() result(nfail)

	character(len = *), parameter :: label = "chapter_2_fft_2"
	character(len = :), allocatable :: fout

	!double precision, parameter :: tol = 1.d-14
	double precision :: fs, dt, norm_diff
	double precision, allocatable :: t(:), r(:), freqs(:), ymag(:)

	double complex, allocatable :: x(:), y(:), xr(:)

	integer :: i, l, ny, fid

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	! Reference for this fft example:
	!
	!     https://www.mathworks.com/help/matlab/ref/fft.html

	!Fs = 1000;            % Sampling frequency                    
	!T = 1/Fs;             % Sampling period       
	!L = 1500;             % Length of signal
	!t = (0:L-1)*T;        % Time vector
	fs = 1000
	dt = 1.d0 / fs

	!l = 2048
	!l = 1500
	!!l = 256
	!l = 1024 * 16
	l = 15000

	t = [(i, i = 0, l-1)] * dt
	!print *, "t = ", t(1: 5), " ... ", t(l-5: l)

	allocate(r(l))
	call random_number(r)
	!x = x + 2*randn(size(t))

	!! Create a signal
	!x = 0.8 + 0.7*sin(2*PI*50*t) + sin(2*PI*120*t)
	x = 0.8 + 0.7*sin(2*PI*3*t) + sin(2*PI*5*t)
	!x = x + 0.8 + 0.7*sin(2*PI*50*t) + sin(2*PI*120*t)

	!! Add some random noise
	!x = x + 0.5 * (2 * (r - 0.5d0))
	x = x + 0.01 * (2 * (r - 0.5d0))

	fout = label // "_tx.out"
	open(newunit = fid, file = fout)
	write(fid, "(a)") "# t, x"
	do i = 1, l
		write(fid, "(3es18.6)") t(i), x(i)
	end do
	close(fid)

	y = fft(x)
	ny = size(y)

	!plot(Fs/L*(0:L-1),abs(Y),"LineWidth",3)
	!title("Complex Magnitude of fft Spectrum")
	!xlabel("f (Hz)")
	!ylabel("|fft(X)|")

	!freqs = fs / l * [(i, i = 0, l-1)]
	freqs = fs / ny * [(i, i = 0, ny-1)]

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

	! TODO: find 3 largest frequencies from ymag.  They should match the input
	! signal (including the dc offset)

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

end function chapter_2_fft_2

!===============================================================================

end module exercises_m

