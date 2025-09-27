
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
	double precision :: x, x2
	double precision, allocatable :: xn(:), xe(:)
	integer :: i, iters, iters2

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

	x  = fzero  (f_nr_ex1, 0.d0, 1000.d0, 1.d-9, iters)
	x2 = fzero77(f_nr_ex1, 0.d0, 1000.d0, 1.d-9, iters2)
	print *, "xf, iters = ", x, iters
	print *, "x2 = ", x2
	call test(x, sqrt(612.d0), 1.d-8, nfail, "fzero 1")
	call test(x, x2, 1.d-200, nfail, "fzero77 1")
	call test(iters, iters2, nfail, "fzero iters 1")
	!stop

	!********

	x = newton_raphson(f_nr_ex2, df_nr_ex2, x0 = 1.d0)
	print *, "x = ", x
	call test(x, 0.865474033d0, 1.d-7, nfail, "newton_raphson 2")

	x = bisect_root(f_nr_ex2, 0.d0, 1000.d0, tol = 1.d-9, maxiters = 50)
	print *, "x = ", x
	call test(x, 0.865474033d0, 1.d-7, nfail, "bisect_root 2")

	x = bisect_root(f_nr_ex2, 0.d0, 2.d0, tol = 1.d-7)
	print *, "x = ", x
	call test(x, 0.865474033d0, 1.d-6, nfail, "bisect_root 2.1")

	x  = fzero  (f_nr_ex2, 0.d0, 1000.d0, 1.d-9, iters)
	x2 = fzero77(f_nr_ex2, 0.d0, 1000.d0, 1.d-9, iters2)
	print *, "xf= ", x
	print *, "x2= ", x2
	call test(x, 0.865474033d0, 1.d-7, nfail, "fzero 2")
	call test(x, x2, 1.d-200, nfail, "fzero77 2")
	call test(iters, iters2, nfail, "fzero iters 2")

	!********

	x = newton_raphson(sin_fn, cos_fn, x0 = 2.0d0)
	print *, "x = ", x
	call test(x, PI, 1.d-9, nfail, "newton_raphson 1.1")

	x = bisect_root(sin_fn, 2.d0, 4.d0, tol = 1.d-7)
	print *, "x = ", x
	call test(x, PI, 1.d-6, nfail, "bisect_root 1.1")

	x = bisect_root(sin_fn, 2.d0, 4.d0, maxiters = 35)
	print *, "x = ", x
	call test(x, PI, 1.d-8, nfail, "bisect_root 1.2")

	x  = fzero  (sin_fn, 2.d0, 4.d0, 1.d-8, iters)
	x2 = fzero77(sin_fn, 2.d0, 4.d0, 1.d-8, iters2)
	print *, "xf= ", x
	print *, "x2= ", x2
	call test(x, PI, 1.d-8, nfail, "fzero 1.2")
	call test(x, x2, 1.d-200, nfail, "fzero77 1.2")
	call test(iters, iters2, nfail, "fzero iters 1.2")

	!********

	x  = fzero  (log_fn, 0.1d0, 10.d0, 1.d-8, iters)
	x2 = fzero77(log_fn, 0.1d0, 10.d0, 1.d-8, iters2)
	print *, "xf= ", x
	print *, "x2= ", x2
	call test(x, 1.d0, 1.d-8, nfail, "fzero 3")
	call test(x, x2, 1.d-200, nfail, "fzero77 3")
	call test(iters, iters2, nfail, "fzero iters 3")

	!********

	do i = 1, 20
		x  = fzero  (wilkinson_fn, 0.1d0 + i - 1, 1.5d0 + i - 1, 1.d-9, iters)
		x2 = fzero77(wilkinson_fn, 0.1d0 + i - 1, 1.5d0 + i - 1, 1.d-9, iters2)
		!print *, "xf, iters = ", x, iters
		!print *, "x2, iters = ", x2, iters2
		call test(x, 0.d0 + i, 1.d-8, nfail, "fzero Wilkinson")
		call test(x, x2, 1.d-200, nfail, "fzero77 Wilkinson")
		call test(iters, iters2, nfail, "fzero iters Wilkinson")
	end do

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
	! Broyden's method, multi-dimensional

	xe = [0.8332816, 3.533462e-002, -0.4985493]  ! expected root

	xn = broyden(f_nr_ex3, df_nr_ex3, x0 = [1.d0, 1.d0, 1.d0])
	print *, "xn = ", xn
	print *, "xe = ", xe
	call test(norm2(xn - xe), 0.d0, 1.d-7, nfail, "broyden 1")

	!********

	xe = [1.0989425808889501d0, 0.36761667884567795d0, 0.14493165687848802d0]
	xn = broyden(f_nr_ex4, df_nr_ex4, x0 = [0.1d0, 0.1d0, 0.1d0])
	print *, "xn = ", xn
	call test(norm2(xn - xe), 0.d0, 1.d-7, nfail, "newton_raphson 4.1")

	!********
	! Test 1D minimization using golden_search

	x = golden_search(sin_fn, 0.d0, 2*PI, 1.d-5)
	print *, "x = ", x
	call test(x, 3*PI/2, 1.d-5, nfail, "golden_search 1")

	x = golden_search(sin_fn, 1*PI, 3*PI, 1.d-5)
	print *, "x = ", x
	call test(x, 3*PI/2, 1.d-5, nfail, "golden_search 1")

	x = golden_search(cos_fn, 0.d0, 2*PI, 1.d-5)
	print *, "x = ", x
	call test(x, PI, 1.d-5, nfail, "golden_search 2")

	x = golden_search(cos_fn, PI, 2.99*PI, 1.d-5)
	print *, "x = ", x
	call test(x, PI, 1.d-5, nfail, "golden_search 2")

	!********
	print *, ""

end function chapter_5_nr

!===============================================================================

double precision function fzero77(f, ax, bx, tol, iters)
	! TODO: test more fzero/fzero77 comparisons
	!
	! This is the original Fortran 77 version.  There is an updated version in
	! the roots module
	!
	! Source:  https://github.com/dr-nikolai/FMM/blob/master/fmm/zeroin.f
	procedure(fn_f64_to_f64) :: f
	double precision ax,bx,tol
	integer, optional, intent(out) :: iters

      double precision  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
	  integer :: iter
!  compute eps, the relative machine precision
      eps = 1.0d0
   10 eps = eps/2.0d0
      tol1 = 1.0d0 + eps
      if (tol1 .gt. 1.0d0) go to 10
!
! initialization
!
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
	  iter = 0
!
! begin step
!
   20 c = a

      fc = fa

      d = b - a
      e = d

   30 continue
		iter = iter + 1
		!print *, "iter 77 = ", iter
      if (dabs(fc) .ge. dabs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
!
! convergence test
!
   40 tol1 = 2.0d0*eps*dabs(b) + 0.5d0*tol
      xm = .5*(c - b)

      if (dabs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0d0) go to 90
!
! is bisection necessary
!
      if (dabs(e) .lt. tol1) go to 70
      if (dabs(fa) .le. dabs(fb)) go to 70
!
! is quadratic interpolation possible
!
      if (a .ne. c) go to 50
!
! linear interpolation
!
      !print *, "linear interpolation 77"
      s = fb/fa
      p = 2.0d0*xm*s
      q = 1.0d0 - s
      go to 60
!
! inverse quadratic interpolation
!
   50 continue
      !print *, "quadratic interpolation 77"
      q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
      q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
!
! adjust signs
!
   60 if (p .gt. 0.0d0) q = -q
      p = dabs(p)
!
! is interpolation acceptable
!
      if ((2.0d0*p) .ge. (3.0d0*xm*q - dabs(tol1*q))) go to 70

      if (p .ge. dabs(0.5d0*e*q)) go to 70
      e = d
      d = p/q
      go to 80
!
! bisection
!

   70 continue
      !print *, "bisection 77"
      d = xm
      e = d
!
! complete step
!
   80 a = b
      fa = fb
      if (dabs(d) .gt. tol1) b = b + d
      if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/dabs(fc))) .gt. 0.0d0) go to 20
      go to 30
! done
   90 fzero77 = b
      if (present(iters)) iters = iter
      return
end function fzero77

!===============================================================================

end module numa__chapter_5

