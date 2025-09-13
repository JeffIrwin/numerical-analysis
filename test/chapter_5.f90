
!> Unit tests for chapter 5: finding zeros and minimum points by iterative
!> methods
module numa__chapter_5

	use numa
	use numa__exercises
	use numa__functions

	implicit none

contains

!===============================================================================

integer function chapter_5_nr() result(nfail)

	character(len = *), parameter :: label = "chapter_5_nr"
	double precision :: x

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********

	!x = newton_raphson(f, df)
	x = newton_raphson(f_nr_ex1, df_nr_ex1, x0 = 1.d0)
	print *, "x = ", x

end function chapter_5_nr

!===============================================================================

end module numa__chapter_5

