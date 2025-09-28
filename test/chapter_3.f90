
!> Unit tests for chapter 3: integration
module numa__chapter_3

	use numa
	use numa__exercises
	use numa__functions

	implicit none

contains

!===============================================================================

integer function chapter_3_newton_cotes() result(nfail)

	! The Newton-Cotes formulas are generalizations of integration rules like
	! the trapezoid rule, and higher-order methods like Simpson's (1/3) rule,
	! Simpson's 3/8 rule, etc.

	character(len = *), parameter :: label = "chapter_3_newton_cotes"

	double precision :: area
	double precision, allocatable :: xi(:), yi(:)

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0
	!********

	area = simpson_13_integrator(sin_fn, 0.d0, PI, PI)
	print *, "area = ", area
	call test(area, 2.094395d0, 1.d-6, nfail, "simpson_13_integrator 1")

	area = simpson_13_integrator(sin_fn, 0.d0, PI, 0.1d0)
	print *, "area = ", area
	call test(area, 2.0000000732694589d0, 1.d-12, nfail, "simpson_13_integrator 2")

	area = simpson_38_integrator(sin_fn, 0.d0, PI, PI)
	print *, "area = ", area
	call test(area, 2.04052428d0, 1.d-6, nfail, "simpson_38_integrator 3")

	area = simpson_38_integrator(sin_fn, 0.d0, PI, 0.1d0)
	print *, "area = ", area
	call test(area, 2.0000000325630984d0, 1.d-12, nfail, "simpson_38_integrator 4")

	!********

	area = milne_integrator(sin_fn, 0.d0, PI, PI)
	print *, "area = ", area
	call test(area, 1.99857d0, 1.d-6, nfail, "milne_integrator 5")

	area = milne_integrator(sin_fn, 0.d0, PI, 0.1d0)
	print *, "area = ", area
	call test(area, 1.9999999999988802d0, 1.d-12, nfail, "milne_integrator 6")

	!********

	area = weddle_integrator(sin_fn, 0.d0, PI, PI)
	print *, "area = ", area
	call test(area, 2.0000178136366551d0, 1.d-12, nfail, "weddle_integrator 5")

	area = weddle_integrator(sin_fn, 0.d0, PI, 0.1d0)
	print *, "area = ", area
	call test(area, 1.9999999999999998d0, 1.d-12, nfail, "weddle_integrator 6")

	area = weddle_integrator(sinx_x, 0.d0, 1.d0, 0.1d0)
	print *, "area = ", area
	call test(area, 0.94608307036718320d0, 1.d-12, nfail, "weddle_integrator 6.1")

	!********
	! Value array integrator(s)

	!print *, "value integrators:"

	xi = range_f64(0.d0, PI, 3)
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0943951023931953d0, 1.d-12, nfail, "simpson_integrator_vals 7")

	xi = range_f64(0.d0, PI, 4)
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0405242847634950d0, 1.d-12, nfail, "simpson_integrator_vals 8")

	xi = range_f64(0.d0, PI, 5)
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0045597549844207d0, 1.d-12, nfail, "simpson_integrator_vals 9")

	xi = range_f64(0.d0, PI, 6)
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0034411937297119d0, 1.d-12, nfail, "simpson_integrator_vals 10")

	xi = range_f64(0.d0, PI, 7)
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0008631896735367d0, 1.d-12, nfail, "simpson_integrator_vals 11")

	xi = range_f64(0.d0, PI, 8)
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0006963918546892d0, 1.d-12, nfail, "simpson_integrator_vals 12")

	xi = range_f64(0.d0, PI, 31)
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0000013379479498d0, 1.d-12, nfail, "simpson_integrator_vals 13")

	xi = range_f64(0.d0, PI, 32)
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0000012070944040d0, 1.d-12, nfail, "simpson_integrator_vals 14")

	xi = range_f64(0.d0, PI, 33)
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0000010333694127d0, 1.d-12, nfail, "simpson_integrator_vals 15")

	xi = range_f64(0.d0, PI, 34)
	yi = sin(xi)
	area = simpson_integrator_vals(yi, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0000009368045015d0, 1.d-12, nfail, "simpson_integrator_vals 16")

	!********
	print *, ""

end function chapter_3_newton_cotes

!===============================================================================

integer function chapter_3_romberg() result(nfail)

	character(len = *), parameter :: label = "chapter_3_romberg"

	double precision :: area
	double precision, allocatable :: areas(:), expect(:)

	integer :: i, n

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0
	!********

	area = romberg_integrator_fixed(sin_fn, 0.d0, PI, 4)
	print *, "area = ", area
	call test(area, 1.9999999945872906d0, 1.d-12, nfail, "romberg_integrator_fixed 1")

	!********
	print *, "integrating sin(x)/x:"

	n = 6
	allocate(areas(n))
	do i = 1, n
		areas(i) = romberg_integrator_fixed(sinx_x, 0.d0, 1.d0, i)
	end do
	print *, "areas = "
	print "(es24.14)", areas

	expect = [ &
		9.46145882273587d-01, &
		9.46083004063674d-01, &
		9.46083070387223d-01, &
		9.46083070367181d-01, &
		9.46083070367183d-01, &
		9.46083070367183d-01  &
	]

	call tests(areas, expect, 1.d-12, nfail, "romberg_integrator_fixed 2")

	!********

	area = romberg_integrator(sin_fn, 0.d0, PI, 1.d-3)
	print *, "area = ", area
	call test(area, 1.9999999945872902d0, 1.d-12, nfail, "romberg_integrator 3")

	area = romberg_integrator(sin_fn, 0.d0, PI, 1.d-6)
	print *, "area = ", area
	call test(area, 2.0000000000013212d0, 1.d-12, nfail, "romberg_integrator 4")

	area = romberg_integrator(sin_fn, 0.d0, PI, 1.d-9)
	print *, "area = ", area
	call test(area, 1.9999999999999991d0, 1.d-12, nfail, "romberg_integrator 5")

	area = romberg_integrator(sinx_x, 0.d0, 1.d0, 1.d-9)
	print *, "area = ", area
	call test(area, 9.46083070367183d-01, 1.d-12, nfail, "romberg_integrator 5.1")

	area = romberg_integrator(fn_example_pg171, 0.d0, PI/2, 1.d-12)
	print *, "area = ", area
	call test(area, 1.d0, 1.d-11, nfail, "romberg_integrator 5.2")

	!********
	! These tests will throw convergence warnings
	print *, GREEN // "Note" // COLOR_RESET // ": the following two test warnings are expected"

	area = romberg_integrator(sin_fn, 0.d0, PI, 0.d0)
	print *, "area = ", area
	call test(area, 2.d0, 1.d-14, nfail, "romberg_integrator 6")

	area = romberg_integrator(sin_fn, 0.d0, PI, tol = 0.d0, max_levels = 4)
	print *, "area = ", area
	call test(area, 1.9999999945872902d0, 1.d-14, nfail, "romberg_integrator 7")

	!********
	print *, ""

end function chapter_3_romberg

!===============================================================================

integer function chapter_3_gauss() result(nfail)

	character(len = *), parameter :: label = "chapter_3_gauss"

	double precision :: area
	!integer :: i

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0
	!********

	! Single-interval Gaussian quadrature isn't very useful, but it's the basis
	! of adaptive Gauss-Kronrod integration
	!
	! Alternatively you could use composite Gaussian integration, which is not
	! implemented here.  The idea is that you statically subdivide the
	! integration interval into a fixed number of subintervals.  The problem is
	! how do you decide how many subintervals to divide without having an error
	! estimate?  I've done composite Newton-Cotes integrators above, but
	! adaptive integrators seem more interesting than static composite
	! integrators

	area = gauss2_single(sin_fn, 0.d0, PI)
	print *, "area = ", area
	call test(area, 1.9358195746511373d0, 1.d-12, nfail, "gauss2_single 1")

	area = gauss3_single(sin_fn, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0013889136077436d0, 1.d-12, nfail, "gauss3_single 2")

	area = gauss4_single(sin_fn, 0.d0, PI)
	print *, "area = ", area
	call test(area, 1.9999842284577216d0, 1.d-12, nfail, "gauss4_single 3")

	area = gauss5_single(sin_fn, 0.d0, PI)
	print *, "area = ", area
	call test(area, 2.0000001102844722d0, 1.d-12, nfail, "gauss5_single 4")

	area = gauss7_single(sin_fn, 0.d0, PI)
	print *, "area g7  = ", area
	call test(area, 2.0000000000017906d0, 1.d-12, nfail, "gauss7_single 5")

	area = kronrod15_single(sin_fn, 0.d0, PI)
	print *, "area k15 = ", area
	call test(area, 1.9999999999999942d0, 1.d-12, nfail, "kronrod15_single 6")

	!do i = 400, 410, 2
	!	area = gauss7_single(sin_fn, 0.d0, i*PI)
	!	print *, "area g7  "//to_str(i)//" = ", area

	!	area = kronrod15_single(sin_fn, 0.d0, i*PI)
	!	print *, "area k15 "//to_str(i)//" = ", area
	!end do

	! Expect 2/3
	area = gauss7_single(sqrt_fn, 0.d0, 1.d0)
	print *, "area g7  sqrt = ", area
	call test(area, 0.66691308508873925d0, 1.d-12, nfail, "gauss7_single 7")

	area = kronrod15_single(sqrt_fn, 0.d0, 1.d0)
	print *, "area k15 sqrt = ", area
	call test(area, 0.66668012554841749d0, 1.d-12, nfail, "kronrod15_single 7")

	!********
	print *, ""

end function chapter_3_gauss

!===============================================================================

integer function chapter_3_adaptive() result(nfail)

	character(len = *), parameter :: label = "chapter_3_adaptive"

	double precision :: area, expect

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0
	!********

	! Simpson
	area = simpson_adaptive_integrator(sin_fn, 0.d0, PI, 1.d-10)
	print *, "area = ", area
	call test(area, 1.9999999999999993d0, 1.d-12, nfail, "simpson_adaptive_integrator 1")

	! Gauss-Kronrod
	area = gk15_adaptive_integrator(sin_fn, 0.d0, PI, 1.d-10)
	print *, "area = ", area
	call test(area, 2.d0, 1.d-12, nfail, "gk15_adaptive_integrator 1")
	print *, ""

	area = simpson_adaptive_integrator(sqrt_fn, 0.d0, 1.d0, 1.d-10)
	print *, "area = ", area
	call test(area, 0.66666666613549719d0, 1.d-12, nfail, "simpson_adaptive_integrator 2")

	area = gk15_adaptive_integrator(sqrt_fn, 0.d0, 1.d0, 1.d-10)
	print *, "area = ", area
	call test(area, 2.d0/3, 1.d-12, nfail, "gk15_adaptive_integrator 2")
	print *, ""

	!********
	! Log from 0 -- nearly impossible!

	! Simpson can't actually evaluate from singularities because it tries to
	! evaluate the end points.  Gauss-Kronrod has an advantage here because
	! Gauss points are always interior, never on the boundary

	area = simpson_adaptive_integrator(log_fn, 1.d-10, 1.d0, 1.d-2, 32)
	print *, "area simp = ", area
	print *, "diff = ", area + 1.d0
	call test(area, -1.0000064653249559d0, 1.d-14, nfail, "gk15_adaptive_integrator 3")
	print *, ""

	area = gk15_adaptive_integrator(log_fn, 0.d0, 1.d0, 1.d-10)
	print *, "area gk15 = ", area
	print *, "diff = ", area + 1.d0
	call test(area, -0.99999997417185815d0, 1.d-14, nfail, "gk15_adaptive_integrator 4")
	print *, ""

	!********
	! 1 / sqrt(x)

	area = gk15_adaptive_integrator(inv_sqrt_fn, 0.d0, 1.d0, 1.d-10)
	!area = gk15_adaptive_integrator(inv_sqrt_fn, 0.d0, 1.d0, 1.d-10, 96)  ! eps underflow
	print *, "area gk15 = ", area
	print *, "diff = ", area - 2.d0
	print *, ""
	call test(area, 1.9998215687092518d0, 1.d-14, nfail, "gk15_adaptive_integrator 5")

	!********
	! Integrate 1 / x ** 2 from a lower bound to infinity

	area = gk15i_adaptive_integrator(inv_square_fn, 1.d0, 1.d-10)
	print *, "area gk15i = ", area
	print *, "diff = ", area - 1.d0
	call test(area, 1.d0, 1.d-14, nfail, "gk15i_adaptive_integrator 5")
	print *, ""

	area = gk15i_adaptive_integrator(inv_square_fn, 2.d0, 1.d-10)
	print *, "area gk15i = ", area
	print *, "diff = ", area - 1.d0/2
	call test(area, 1.d0/2, 1.d-14, nfail, "gk15i_adaptive_integrator 6")
	print *, ""

	area = gk15i_adaptive_integrator(inv_square_fn, 3.d0, 1.d-10)
	print *, "area gk15i = ", area
	print *, "diff = ", area - 1.d0/3
	call test(area, 1.d0/3, 1.d-14, nfail, "gk15i_adaptive_integrator 7")
	print *, ""

	!********

	area = gk15i_adaptive_integrator(exp_nx2, 1.d0, 1.d-10)
	expect = 0.5d0 * sqrt(PI) * erfc(1.d0)
	print *, "area exp_nx2 = ", area
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15i_adaptive_integrator 8")
	print *, ""

	area = gk15i_adaptive_integrator(exp_nx2, 2.d0, 1.d-10)
	expect = 0.5d0 * sqrt(PI) * erfc(2.d0)
	print *, "area exp_nx2 = ", area
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15i_adaptive_integrator 9")
	print *, ""

	area = gk15i_adaptive_integrator(exp_nx2, 3.d0, 1.d-10)
	expect = 0.5d0 * sqrt(PI) * erfc(3.d0)
	print *, "area exp_nx2 = ", area
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15i_adaptive_integrator 10")
	print *, ""

	area = gk15i_adaptive_integrator(exp_nx2, -1.d0, 1.d-10)
	expect = 0.5d0 * sqrt(PI) * erfc(-1.d0)
	print *, "area exp_nx2 = ", area
	print *, "expect       = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15i_adaptive_integrator 11")
	print *, ""

	area = gk15i_adaptive_integrator(exp_nx2, 0.d0, 1.d-10)
	expect = 0.5d0 * sqrt(PI)
	print *, "area exp_nx2 = ", area
	print *, "expect       = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15i_adaptive_integrator 12")
	print *, ""

	!********
	! Integrate from -infinity to an upper bound

	area = gk15ni_adaptive_integrator(exp_nx2, -2.d0, 1.d-10)
	expect = 0.5d0 * sqrt(PI) * erfc(2.d0)
	print *, "area exp_nx2 = ", area
	print *, "expect       = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15ni_adaptive_integrator 13")
	print *, ""

	area = gk15ni_adaptive_integrator(exp_nx2, -1.d0, 1.d-10)
	expect = 0.5d0 * sqrt(PI) * erfc(1.d0)
	print *, "area exp_nx2 = ", area
	print *, "expect       = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15ni_adaptive_integrator 14")
	print *, ""

	area = gk15ni_adaptive_integrator(exp_nx2, 0.d0, 1.d-10)
	expect = 0.5d0 * sqrt(PI)
	print *, "area exp_nx2 = ", area
	print *, "expect       = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15ni_adaptive_integrator 15")
	print *, ""

	area = gk15ni_adaptive_integrator(exp_nx2, 1.d0, 1.d-10)
	expect = 0.5d0 * sqrt(PI) * erfc(-1.d0)
	print *, "area exp_nx2 = ", area
	print *, "expect       = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15ni_adaptive_integrator 16")
	print *, ""

	area = gk15ni_adaptive_integrator(exp_nx2, 2.d0, 1.d-10)
	expect = 0.5d0 * sqrt(PI) * erfc(-2.d0)
	print *, "area exp_nx2 = ", area
	print *, "expect       = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15ni_adaptive_integrator 17")
	print *, ""

	!********

	area = gk15ii_adaptive_integrator(exp_nx2, 1.d-10)
	expect = sqrt(PI)
	print *, "area exp_nx2 = ", area
	print *, "expect       = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15ii_adaptive_integrator 17")
	print *, ""

	area = gk15i_adaptive_integrator(inv_1px2, 0.d0, 1.d-10)
	expect = 0.5d0 * PI
	print *, "area   = ", area
	print *, "expect = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15i_adaptive_integrator 18")
	print *, ""

	area = gk15ii_adaptive_integrator(inv_1px2, 1.d-10)
	expect = PI
	print *, "area   = ", area
	print *, "expect = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15ii_adaptive_integrator 19")
	print *, ""

	area = gk15i_adaptive_integrator(inv_1px1p5, 0.d0, 1.d-10)
	expect = 2.4182207210215423d0
	print *, "area   = ", area
	print *, "expect = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15i_adaptive_integrator 20")
	print *, ""

	area = gk15ii_adaptive_integrator(inv_1px1p5, 1.d-10)
	expect = 2 * 2.4182207210215423d0
	print *, "area   = ", area
	print *, "expect = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15ii_adaptive_integrator 21")
	print *, ""

	area = gk15i_adaptive_integrator(inv_1px4, 0.d0, 1.d-10)
	expect = 0.5d0 * PI / sqrt(2.d0)
	print *, "area   = ", area
	print *, "expect = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15i_adaptive_integrator 22")
	print *, ""

	area = gk15ii_adaptive_integrator(inv_1px4, 1.d-10)
	expect = PI / sqrt(2.d0)
	print *, "area   = ", area
	print *, "expect = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15ii_adaptive_integrator 23")
	print *, ""

	!********
	print *, "Bessel's integral:"

	! TODO: the fn interface `fn_f64_to_f64` is somewhat limited.  If you wanted
	! to roll your own Bessel fn calculator, you would need a different fn to
	! integrate for every real value of x (and int n)
	!
	! Maybe there should be another gk15 integrator that takes a "user data"
	! array of doubles (and ints).  This would need another fn interface, say
	! fn_f64_to_f64_dble_data, which also takes the data
	!
	! Actually I have created exactly this in fn_f64_params_to_f64.  A new
	! gk15_adaptive_integrator() is still needed

	area = gk15_adaptive_integrator(bessel_3_2p5, 0.d0, PI, 1.d-10)
	expect = bessel_jn(3, 2.5d0)
	print *, "area   = ", area
	print *, "expect = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15_adaptive_integrator 24")
	print *, ""

	area = gk15_integrator_params(bessel_3, [2.5d0], 0.d0, PI, 1.d-10)
	expect = bessel_jn(3, 2.5d0)
	print *, "area   = ", area
	print *, "expect = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15 bessel(3, 2.5)")
	print *, ""

	! Could also add integrators that take int params, e.g. for generalization
	! beyond bessel_3()
	area = gk15_integrator_params(bessel_3, [2.8d0], 0.d0, PI, 1.d-10)
	expect = bessel_jn(3, 2.8d0)
	print *, "area   = ", area
	print *, "expect = ", expect
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15 bessel(3, 2.8)")
	print *, ""

	!********

	area = gk15i_adaptive_integrator(gamma_int_5, 0.d0, 1.d-10)
	expect = 24.d0  ! gamma(5) == factorial(4)
	print *, "area = ", area
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15i gamma 5")
	print *, ""

	area = gk15i_adaptive_integrator(gamma_int_6, 0.d0, 1.d-10)
	expect = 120.d0
	print *, "area = ", area
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15i gamma 6")
	print *, ""

	area = gk15i_integrator_params(gamma_int, [5.d0], 0.d0, 1.d-10)
	expect = 24.d0
	print *, "area = ", area
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15i gamma params 5")
	print *, ""

	area = gk15i_integrator_params(gamma_int, [5.5d0], 0.d0, 1.d-10)
	expect = gamma(5.5d0)
	print *, "area = ", area
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15i gamma params 5.5")
	print *, ""

	area = gk15i_integrator_params(gamma_int, [6.d0], 0.d0, 1.d-10)
	expect = 120.d0
	print *, "area = ", area
	print *, "diff = ", area - expect
	call test(area, expect, 1.d-10, nfail, "gk15i gamma params 5")
	print *, ""

	!********

	! This could always be expanded with more tests.  Quadpack has some
	! interesting tests:
	!
	!     	https://github.com/jacobwilliams/quadpack/blob/master/test/quadpack_test_module.F90

	!********

end function chapter_3_adaptive

!===============================================================================

end module numa__chapter_3

