
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
	double precision, allocatable :: xn(:), xe(:)

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! 1D

	x = newton_raphson(f_nr_ex1, df_nr_ex1, x0 = 1.d0)
	print *, "x = ", x
	call test(x, sqrt(612.d0), 1.d-9, nfail, "newton_raphson 1")

	x = bisect_root(f_nr_ex1, 0.d0, 1000.d0, tol = 1.d-9, maxiters = 50)
	print *, "x = ", x
	call test(x, sqrt(612.d0), 1.d-8, nfail, "bisect_root 1.1")

	x = fzero(f_nr_ex1, 0.d0, 1000.d0, 1.d-9)
	print *, "x = ", x
	call test(x, sqrt(612.d0), 1.d-8, nfail, "fzero 1")
	!stop

	!********

	x = newton_raphson(sin_fn, cos_fn, x0 = 2.0d0)
	print *, "x = ", x
	call test(x, PI, 1.d-9, nfail, "newton_raphson 1.1")

	!********

	x = newton_raphson(f_nr_ex2, df_nr_ex2, x0 = 1.d0)
	print *, "x = ", x
	call test(x, 0.865474033d0, 1.d-7, nfail, "newton_raphson 2")

	x = bisect_root(f_nr_ex2, 0.d0, 1000.d0, tol = 1.d-9, maxiters = 50)
	print *, "x = ", x
	call test(x, 0.865474033d0, 1.d-7, nfail, "bisect_root 2")

	x = fzero(f_nr_ex2, 0.d0, 1000.d0, 1.d-9)
	print *, "x = ", x
	call test(x, 0.865474033d0, 1.d-7, nfail, "fzero 2")
	!stop

	!********

	x = bisect_root(sin_fn, 2.d0, 4.d0, tol = 1.d-7)
	print *, "x = ", x
	call test(x, PI, 1.d-6, nfail, "bisect_root 1.1")

	x = bisect_root(sin_fn, 2.d0, 4.d0, maxiters = 35)
	print *, "x = ", x
	call test(x, PI, 1.d-8, nfail, "bisect_root 1.2")

	!********

	x = bisect_root(f_nr_ex2, 0.d0, 2.d0, tol = 1.d-7)
	print *, "x = ", x
	call test(x, 0.865474033d0, 1.d-6, nfail, "bisect_root 2.1")

	!********
	! Multi-dimensional

	xe = [0.8332816, 3.533462e-002, -0.4985493]  ! expected root

	xn = newton_raphson_nd(f_nr_ex3, df_nr_ex3, x0 = [1.d0, 1.d0, 1.d0])
	!xn = newton_raphson_nd(f_nr_ex3, df_nr_ex3, nx = 3)
	!xn = newton_raphson_nd(f_nr_ex3, df_nr_ex3)
	print *, "xn = ", xn
	call test(norm2(xn - xe), 0.d0, 1.d-7, nfail, "newton_raphson 3")

	!********

	xe = [1.0989425808889501d0, 0.36761667884567795d0, 0.14493165687848802d0]
	xn = newton_raphson_nd(f_nr_ex4, df_nr_ex4, nx = 3)
	print *, "xn = ", xn
	call test(norm2(xn - xe), 0.d0, 1.d-7, nfail, "newton_raphson 4.1")

	!********
	! Same as above but with different initial guess, and thus a different root
	! solution

	xe = ones(3)
	xn = newton_raphson_nd(f_nr_ex4, df_nr_ex4, [1.d0, 2.d0, 3.d0])
	print *, "xn = ", xn
	call test(norm2(xn - xe), 0.d0, 1.d-7, nfail, "newton_raphson 4.1")

	!********
	print *, ""

end function chapter_5_nr

!===============================================================================

end module numa__chapter_5

