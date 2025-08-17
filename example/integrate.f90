
module numa__example_integrator
	implicit none

	contains
		double precision function log_fn(x)
			double precision, intent(in) :: x
			log_fn = log(x)
		end function log_fn

end module numa__example_integrator

program example

	use numa, only:  gk15_adaptive_integrator
	use numa__example_integrator
	implicit none

	double precision :: area

	print *, "starting example integrate.f90"

	! Integrate log_fn from 0 to 1 with a tolerance of 0.01
	area = gk15_adaptive_integrator(log_fn, 0.d0, 1.d0, 1.d-2)
	print *, "area = ", area
	print *, "diff = ", area - (-1.d0)

	! Integrate log_fn from 0 to 1 with a tolerance of 0.001 and up to 32 recursion levels
	area = gk15_adaptive_integrator(log_fn, 0.d0, 1.d0, 1.d-3, 32)
	print *, "area = ", area
	print *, "diff = ", area - (-1.d0)

	print *, "ending example integrate.f90"
	print *, ""

end program example

