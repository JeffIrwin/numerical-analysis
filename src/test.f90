
module test_m

	use exercises_m
	use numerical_analysis_m
	use utils_m

end module test_m

program main

	use test_m

	implicit none

	integer :: nfail

	write(*,*) MAGENTA // "Starting numerical-analysis test" // COLOR_RESET

	!print *, "PI = ", PI

	nfail = 0
	nfail = nfail + chapter_2_exercise_2()
	nfail = nfail + chapter_2_example_2p2p4()

	if (nfail == 0) then
		write(*,*) GREEN // "Success!" // COLOR_RESET
	else
		write(*,*) ERROR // to_str(nfail) // " test(s) failed!" // COLOR_RESET
	end if

	write(*,*) MAGENTA // "Ending numerical-analysis test" // COLOR_RESET
	call exit(nfail)

end program main

