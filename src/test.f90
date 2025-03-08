
module test_m
end module test_m

program test

	use numerical_analysis_m
	use exercises_m
	use test_m

	implicit none

	integer :: io

	write(*,*) "Starting numerical-analysis test"

	!print *, "PI = ", PI

	io = 0
	io = io + chapter_2_exercise_2()

	write(*,*) "Ending numerical-analysis test"
	call exit(io)

end program test

