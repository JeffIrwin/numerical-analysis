
module test_m
end module test_m

program test

	use numerical_analysis_m
	use test_m
	implicit none

	integer :: io

	write(*,*) "Starting numerical-analysis test"

	!print *, "PI = ", PI

	io = chapter_2_exercise_2()

	write(*,*) "Ending numerical-analysis test"

end program test

