
!> Main program module for the unit tests
module numa__test

	use numa
	use numa__exercises
	use numa__utils

end module numa__test

!> Main program for the unit tests
program main

	use numa__test

	implicit none

	integer :: nfail

	real :: t0, t

	call cpu_time(t0)

	write(*,*) MAGENTA // "Starting numerical-analysis test" // COLOR_RESET

	nfail = 0

	!-------------------------------------------------------
	! Test some pre-requisite basic functionality before getting into the actual
	! numerical analysis
	nfail = nfail + test_bad_numa_usage()

	!-------------------------------------------------------
	! Chapter 1: Error Analysis (no exercises here)

	!-------------------------------------------------------
	! Chapter 2: Interpolation
	nfail = nfail + chapter_2_exercise_2()
	nfail = nfail + chapter_2_example_2p2p4()
	nfail = nfail + chapter_2_fft_1()
	nfail = nfail + chapter_2_fft_2()

	! These might belong in chapter 4, but the tridiagonal matrix algorithm is a
	! prerequisite for cubic splines, and banded matrices are a small extension
	! of tridiagonal matrices
	nfail = nfail + chapter_2_tridiag()
	nfail = nfail + chapter_2_tridiag_corner()
	nfail = nfail + chapter_2_banded()
	nfail = nfail + chapter_2_cubic_splines()
	nfail = nfail + chapter_2_bezier_splines()

	!-------------------------------------------------------
	! Chapter 3: Topics in Integration
	nfail = nfail + chapter_3_newton_cotes()
	nfail = nfail + chapter_3_romberg()
	nfail = nfail + chapter_3_gauss()
	nfail = nfail + chapter_3_adaptive()

	!-------------------------------------------------------
	! Chapter 6: Eigenvalue Problems
	! TODO: move in chapter order

	! TODO: add unit tests for sorting routines

	nfail = nfail + chapter_6_basic_qr()
	nfail = nfail + chapter_6_hessenberg_qr()
	nfail = nfail + chapter_6_francis_qr()

	!-------------------------------------------------------
	! Chapter 4: Systems of Linear Equations
	nfail = nfail + chapter_4_lu()
	nfail = nfail + chapter_4_lu_c64()
	nfail = nfail + chapter_4_inv()
	nfail = nfail + chapter_4_cholesky()
	nfail = nfail + chapter_4_qr()
	nfail = nfail + chapter_4_qr_c64()
	nfail = nfail + chapter_4_gram_schmidt()
	nfail = nfail + chapter_4_lls()

	!-------------------------------------------------------
	! Chapter 5: Finding Zeros and Minimum Points by Iterative Methods

	!-------------------------------------------------------
	! Chapter 7: Ordinary Differential Equations

	!-------------------------------------------------------
	! Chapter 8: Iterative Methods for the Solution of Large Systems of Linear
	! Equations

	!-------------------------------------------------------
	if (nfail == 0) then
		write(*,*) GREEN // "Success!" // COLOR_RESET
	else
		write(*,*) ERROR // to_str(nfail) // " test(s) failed!" // COLOR_RESET
	end if

	write(*,*) MAGENTA // "Ending numerical-analysis test" // COLOR_RESET
	call cpu_time(t)
	write(*,*) "CPU time = ", to_str(t - t0), " s"

	call exit(nfail)

end program main

