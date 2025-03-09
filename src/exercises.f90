
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

integer function chapter_2_fft() result(nfail)

	double precision, parameter :: tol = 1.d-14
	double precision :: ts, freq, amp1, amp4, amp7, norm_diff
	double precision, allocatable :: t(:)

	double complex, allocatable :: x(:), x0(:), xx(:), xx_oneside(:)

#if defined(__INTEL_COMPILER)
	double precision, external :: nrm2, dznrm2
#endif

	integer :: i, sr, n, n_oneside

	write(*,*) CYAN // "Starting chapter_2_fft()" // COLOR_RESET

	nfail = 0

	sr = 256  ! sampling rate
	n = sr    ! TODO: this needs be be disentangled for better testing

	ts = 1.0d0 / sr  ! sampling interval
	!t = np.arange(0,1,ts)
	t = [(i, i = 0, sr-1)] * ts
	!print *, "t = ", t

	freq = 1.d0 ; amp1 = 3.d0
	x = amp1 * sin(2 * PI * freq * t)
	!print *, "x = ", x

	freq = 4.d0 ; amp4 = 1.d0
	x = x + amp4 * sin(2 * PI * freq * t)

	freq = 7.d0 ; amp7 = 0.5d0
	x = x + amp7 * sin(2 * PI * freq * t)

	! Backup original signal to compare with round trip fft() then ifft()
	x0 = x

	! Fortran has no complex edit descriptor, so print as pairs of reals
	print *, "x = "
	!print "(2es18.6)", x(1: 10)
	print "(2es18.6)", x
	print *, ""

	!xx = fft_rec(x)
	xx = fft(x)

	!print *, "abs(xx) = ", abs(xx(1: 10))
	print *, "xx = "
	print "(2es18.6)", xx

	! Normalize the amplitude and truncate to one side
	n_oneside = size(x) / 2
	xx_oneside = xx(1: n_oneside) / n_oneside
	!print *, "abs(xx_oneside) = ", abs(xx_oneside)
	!print *, "abs(xx_oneside) = ", abs(xx_oneside(1: 10))
	print *, "abs(xx_oneside) = "
	print "(es18.6)", abs(xx_oneside(1: 10))

	! TODO: generalize.  These indices of xx_oneside assume that the sample period is 1
	call test(abs(xx_oneside(1+1)), amp1, tol, nfail, "fft()")
	call test(abs(xx_oneside(4+1)), amp4, tol, nfail, "fft()")
	call test(abs(xx_oneside(7+1)), amp7, tol, nfail, "fft()")

	print *, "sum abs = ", sum(abs(xx_oneside))
	call test(sum(abs(xx_oneside)), amp1+amp4+amp7, tol, nfail, "fft()")
	print *, ""

	!return  ! TODO

	x = ifft(xx) / sr
	!x = fft(xx) / sr

	print *, "x (round trip) = "
	print "(2es18.6)", x(1: 10)
	!print *, "x = ", x
	!print "(a,20es18.6)", "x = ", x(1: 10)
	print *, ""

	! Get the norm of the difference between the original signal x0 and the
	! round-trip fft() + ifft() signal x

	!print *, "norm diff = ", norm2(x - x0)
	!print *, "norm diff = ", nrm2(x - x0)
	!print *, "norm diff = ", nrm2(n, x - x0, 1)  ! undefined

#if defined(__INTEL_COMPILER)
	print *, "norm diff = ", dznrm2(n, x - x0, 1), " (mkl)"
#endif
	norm_diff = real(sqrt(dot_product(x - x0, x - x0)))
	print *, "norm_diff = ", norm_diff
	call test(norm_diff, 0.d0, 1.d-12, nfail, "fft()")

	!! dot_product() already performs a conjg.  D'oh!
	!print *, "norm diff = ", abs(sqrt(dot_product(conjg(x - x0), x - x0)))
	!print *, "norm diff = ", sqrt(dot_product(conjg(x - x0), x - x0))
	!print *, "norm diff = ", sqrt(dot_product(x - x0, conjg(x - x0)))

	!# calculate the frequency
	!N = len(xx)
	!n = np.arange(N)
	!T = N/sr
	!freq = n/T 
	!
	!# Get the one-sided specturm
	!n_oneside = N//2
	!
	!# Get the one side frequency
	!f_oneside = freq[:n_oneside]
	!
	!# normalize the amplitude
	!xx_oneside =xx[:n_oneside]/n_oneside
	!
	!#print("freq = ", freq)
	!#print("abs(xx) = ", abs(xx))
	!
	!#print("abs(xx_oneside) = ", abs(xx_oneside))
	!print("abs(xx_oneside) = ", abs(xx_oneside[0:10]))

end function chapter_2_fft

!===============================================================================

end module exercises_m

