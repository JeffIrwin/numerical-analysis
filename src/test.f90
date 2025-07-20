
module numa__test

	use numa
	use numa__exercises
	use numa__utils

end module numa__test

program main

	use numa__test

	implicit none

	integer :: nfail

	write(*,*) MAGENTA // "Starting numerical-analysis test" // COLOR_RESET

	nfail = 0

	! Chapter 1: Error Analysis (no exercises here)

	! Chapter 2: Interpolation
	nfail = nfail + chapter_2_exercise_2()
	nfail = nfail + chapter_2_example_2p2p4()
	nfail = nfail + chapter_2_fft_1()
	nfail = nfail + chapter_2_fft_2()

	! TODO: put in chapter order.  But for development logging, it's easier to
	! have the thing I'm working get tested last
	nfail = nfail + chapter_4_lu()

	nfail = nfail + chapter_2_tridiag()
	nfail = nfail + chapter_2_tridiag_corner()
	nfail = nfail + chapter_2_banded()
	nfail = nfail + chapter_2_cubic_splines()
	nfail = nfail + chapter_2_bezier_splines()

	! Chapter 3: Topics in Integration
	nfail = nfail + chapter_3_newton_cotes()
	nfail = nfail + chapter_3_romberg()
	nfail = nfail + chapter_3_gauss()

	! Chapter 4: Systems of Linear Equations

	! Chapter 5: Finding Zeros and Minimum Points by Iterative Methods

	! Chapter 6: Eigenvalue Problems

	! Chapter 7: Ordinary Differential Equations

	! Chapter 8: Iterative Methods for the Solution of Large Systems of Linear
	! Equations

	if (nfail == 0) then
		write(*,*) GREEN // "Success!" // COLOR_RESET
	else
		write(*,*) ERROR // to_str(nfail) // " test(s) failed!" // COLOR_RESET
	end if

	write(*,*) MAGENTA // "Ending numerical-analysis test" // COLOR_RESET
	call exit(nfail)

end program main

