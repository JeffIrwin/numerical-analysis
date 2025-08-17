
program example

	use numa, only:  invmul, inv
	implicit none

	double precision, allocatable :: a(:,:), x(:), b(:)

	print *, "starting example matrix-solve.f90"

	a = reshape([ &
		-1,  2,  3,  4, &
		 2, -3,  4,  1, &
		 3,  4, -1,  2, &
		 4,  1,  2, -3  &
		], &
		[4, 4] &
	)

	b = [13, 12, 11, 10]

	! Solve for x in the matrix equation a*x == b
	x = invmul(a, b)

	!! Alternatively, use more expensive explicit inverse `inv()`
	!x = matmul(inv(a), b)

	print "(a,*(es16.6))", " x = ", x

	print *, "diff = ", norm2(b - matmul(a, x))

	print *, "ending example matrix-solve.f90"
	print *, ""

end program example

