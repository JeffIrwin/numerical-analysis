
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
	nfail = nfail + chapter_2_b_splines()

	if (nfail == 0) then
		write(*,*) GREEN // "Success!" // COLOR_RESET
	else
		write(*,*) ERROR // to_str(nfail) // " test(s) failed!" // COLOR_RESET
	end if

	write(*,*) MAGENTA // "Ending numerical-analysis test" // COLOR_RESET
	call exit(nfail)

end program main

