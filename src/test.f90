
module test_m

	use exercises_m
	use numerical_analysis_m
	use utils_m

end module test_m

program test

	use test_m

	implicit none

	integer :: io

	write(*,*) MAGENTA // "Starting numerical-analysis test" // COLOR_RESET

	!print *, "PI = ", PI

	io = 0
	io = io + chapter_2_exercise_2()

	if (io == 0) then
		write(*,*) GREEN // "Success!" // COLOR_RESET
	else
		write(*,*) RED // "Error!" // COLOR_RESET
	end if

	write(*,*) MAGENTA // "Ending numerical-analysis test" // COLOR_RESET
	call exit(io)

end program test

