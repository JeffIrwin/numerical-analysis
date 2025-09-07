
!> Unit tests and exercises
module numa__exercises

	use numa
	use numa__functions

	implicit none

	interface test
		procedure :: test_f64
		procedure :: test_i32
		procedure :: test_bool
	end interface test

contains

!===============================================================================

integer function assert(condition)
	! Return 0 if true or 1 if false
	logical, intent(in) :: condition
	if (condition) then
		assert = 0
	else
		assert = 1
	end if
end function assert

!===============================================================================

! This could be a macro like PANIC, but it has so many arguments that there
! might be line truncation issues
subroutine test_f64(val, expect, tol, nfail, msg)
	double precision, intent(in) :: val, expect, tol
	integer, intent(inout) :: nfail
	character(len = *), intent(in) :: msg

	!********
	integer :: io

	io = assert(abs(val - expect) < tol)
	nfail = nfail + io
	if (io /= 0) then

		write(*,*) ERROR // "assert failed for value " // &
			to_str(val) // " in "// msg

		!! Could be helpful sometimes, but maybe too noisy.  Maybe add a cmd arg
		!! to enable backtrace for unit test failures.  Also this is gfortran
		!! extension, intel has `tracebackqq()`
		!call backtrace()

	end if

end subroutine test_f64

subroutine test_i32(val, expect, nfail, msg)
	integer, intent(in) :: val, expect
	integer, intent(inout) :: nfail
	character(len = *), intent(in) :: msg

	!********
	integer :: io

	io = assert(val == expect)
	nfail = nfail + io
	if (io /= 0) then

		write(*,*) ERROR // "assert failed for value " // &
			to_str(val) // " in "// msg

		!! Could be helpful sometimes, but maybe too noisy.  Maybe add a cmd arg
		!! to enable backtrace for unit test failures.  Also this is gfortran
		!! extension, intel has `tracebackqq()`
		!call backtrace()

	end if

end subroutine test_i32

subroutine test_bool(val, expect, nfail, msg)
	logical, intent(in) :: val, expect
	integer, intent(inout) :: nfail
	character(len = *), intent(in) :: msg

	!********
	integer :: io

	io = assert(val .eqv. expect)
	nfail = nfail + io
	if (io /= 0) then

		write(*,*) ERROR // "assert failed for value " // &
			to_str(val) // " in "// msg

		!! Could be helpful sometimes, but maybe too noisy.  Maybe add a cmd arg
		!! to enable backtrace for unit test failures.  Also this is gfortran
		!! extension, intel has `tracebackqq()`
		!call backtrace()

	end if

end subroutine test_bool

!********

subroutine test_ne(val, unexpect, nfail, msg)
	integer, intent(in) :: val, unexpect
	integer, intent(inout) :: nfail
	character(len = *), intent(in) :: msg

	!********
	integer :: io

	io = assert(val /= unexpect)
	nfail = nfail + io
	if (io /= 0) then
		write(*,*) ERROR // "assert failed for value " // &
			to_str(val) // " in "// msg
	end if

end subroutine test_ne

!********

subroutine tests(vals, expects, tol, nfail, msg)
	double precision, intent(in) :: vals(:), expects(:), tol
	integer, intent(inout) :: nfail
	character(len = *), intent(in) :: msg
	!********
	integer :: i
	do i = 1, size(vals, 1)
		call test(vals(i), expects(i), tol, nfail, msg)
	end do
end subroutine tests

!===============================================================================

integer function test_compiler_sanity() result(nfail)

	character(len = *), parameter :: label = "test_compiler_sanity"

	character(len = :), allocatable :: link
	double precision, allocatable :: a(:)
	integer :: n, i1(1)
	logical, allocatable :: mask(:)

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0
	!-----------------------------------

	n = 5
	allocate(a(n), mask(n))
	a = 1.d0
	mask = .false.

	i1 = minloc(a, mask)
	if (i1(1) /= 0) then

		! Actually exit because other tests need a conforming minloc() fn.  At
		! least, linprog() needs it

		write(*,*) RED // "Catastrophic error" // COLOR_RESET // &
			":  minloc() does not return 0 for all false mask."
		write(*,*) "If using an Intel compiler, use `-standard-semantics`."
		write(*,*) "See the discussion at the following link:"
		write(*,*)
		link = "https://community.intel.com/t5/Intel-Fortran-Compiler/Incorrect-minloc-maxloc-results-with-ifx-2024-0/m-p/1549936"
		write(*,"(a)") "    " // link
		write(*,*)

		call exit(1)

	end if

end function test_compiler_sanity

!===============================================================================

integer function test_basics() result(nfail)

	character(len = *), parameter :: label = "test_basics"

	double precision, allocatable :: v(:)

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0
	!-----------------------------------
	! Test range (linspace) with count and step overloads

	! Counts

	v = range_f64(0.d0, 10.d0, 11)
	call test(v(1), 0.d0, 1.d-6, nfail, "range_f64 1.1")
	call test(v(size(v)), 10.d0, 1.d-6, nfail, "range_f64 1.2")
	call test(size(v), 11, nfail, "range_f64 1.3")

	v = range_f64(4.d0, 10.d0, 7)
	call test(v(1), 4.d0, 1.d-6, nfail, "range_f64 2.1")
	call test(v(size(v)), 10.d0, 1.d-6, nfail, "range_f64 2.2")
	call test(size(v), 7, nfail, "range_f64 2.3")

	v = range_f64(0.d0, 1.d0, 11)
	call test(v(1), 0.d0, 1.d-6, nfail, "range_f64 3.1")
	call test(v(size(v)), 1.d0, 1.d-6, nfail, "range_f64 3.2")
	call test(size(v), 11, nfail, "range_f64 3.3")

	v = range_f64(0.4d0, 1.d0, 7)
	call test(v(1), 0.4d0, 1.d-6, nfail, "range_f64 4.1")
	call test(v(size(v)), 1.d0, 1.d-6, nfail, "range_f64 4.2")
	call test(size(v), 7, nfail, "range_f64 4.3")

	!********
	! Steps

	v = range_f64(0.d0, 10.d0, 1.d0)
	call test(v(1), 0.d0, 1.d-6, nfail, "range_f64 5.1")
	call test(v(size(v)), 10.d0, 1.d-6, nfail, "range_f64 5.2")
	call test(v(2) - v(1), 1.d0, 1.d-6, nfail, "range_f64 5.3")

	v = range_f64(4.d0, 10.d0, 1.d0)
	call test(v(1), 4.d0, 1.d-6, nfail, "range_f64 6.1")
	call test(v(size(v)), 10.d0, 1.d-6, nfail, "range_f64 6.2")
	call test(v(2) - v(1), 1.d0, 1.d-6, nfail, "range_f64 6.3")

	v = range_f64(0.d0, 1.d0, 0.1d0)
	call test(v(1), 0.d0, 1.d-6, nfail, "range_f64 7.1")
	call test(v(size(v)), 1.d0, 1.d-6, nfail, "range_f64 7.2")
	call test(v(2) - v(1), 0.1d0, 1.d-6, nfail, "range_f64 7.3")

	v = range_f64(0.4d0, 1.d0, 0.1d0)
	call test(v(1), 0.4d0, 1.d-6, nfail, "range_f64 8.1")
	call test(v(size(v)), 1.d0, 1.d-6, nfail, "range_f64 8.2")
	call test(v(2) - v(1), 0.1d0, 1.d-6, nfail, "range_f64 8.3")

	!********
	! Steps with step sizes such that the last element does not always exactly
	! hit the stop value

	v = range_f64(0.d0, 10.d0, 3.d0)
	call test(v(1), 0.d0, 1.d-6, nfail, "range_f64 9.1")
	call test(v(size(v)), 9.d0, 1.d-6, nfail, "range_f64 9.2")
	call test(v(2) - v(1), 3.d0, 1.d-6, nfail, "range_f64 9.3")

	v = range_f64(4.d0, 10.d0, 3.d0)
	!print *, "v = ", v
	call test(v(1), 4.d0, 1.d-6, nfail, "range_f64 10.1")
	call test(v(size(v)), 10.d0, 1.d-6, nfail, "range_f64 10.2")
	call test(v(2) - v(1), 3.d0, 1.d-6, nfail, "range_f64 10.3")

	v = range_f64(0.d0, 1.d0, 0.3d0)
	call test(v(1), 0.d0, 1.d-6, nfail, "range_f64 11.1")
	call test(v(size(v)), 0.9d0, 1.d-6, nfail, "range_f64 11.2")
	call test(v(2) - v(1), 0.3d0, 1.d-6, nfail, "range_f64 11.3")

	v = range_f64(0.4d0, 1.d0, 0.3d0)
	call test(v(1), 0.4d0, 1.d-6, nfail, "range_f64 12.1")
	call test(v(size(v)), 1.d0, 1.d-6, nfail, "range_f64 12.2")
	call test(v(2) - v(1), 0.3d0, 1.d-6, nfail, "range_f64 12.3")

	!********

end function test_basics

!===============================================================================

integer function test_bad_numa_usage() result(nfail)

	character(len = *), parameter :: label = "test_bad_numa_usage"

	integer :: io

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0
	!-----------------------------------

	b_lu_factor: block
	double precision, allocatable :: a(:,:)
	integer, allocatable :: pivot(:)

	!********
	! Pivot not allocated

	!call lu_factor(a, pivot)              ! this will stop in lu_factor()
	call lu_factor(a, pivot, iostat = io)  ! this will let caller handle error
	call test_ne(io, 0, nfail, "lu_factor pivot unallocated")

	!********
	! Singular matrix `a`
	allocate(a(3, 3))
	pivot = [1, 2, 3]
	a = 0

	call lu_factor(a, pivot, iostat = io)
	call test_ne(io, 0, nfail, "lu_factor singular matrix")

	!********
	! Non-square `a`
	deallocate(a)
	allocate(a(3, 2))
	call lu_factor(a, pivot, iostat = io)
	call test_ne(io, 0, nfail, "lu_factor non-square matrix")

	!********

	deallocate(pivot)
	!call lu_factor(a, pivot)     ! this will stop in lu_factor()

	end block b_lu_factor
	!-----------------------------------

	b_lu_factor_c64: block
	double complex, allocatable :: a(:,:)
	integer, allocatable :: pivot(:)

	!********
	! Pivot not allocated

	call lu_factor(a, pivot, iostat = io)
	call test_ne(io, 0, nfail, "lu_factor pivot unallocated")

	!********
	! Singular matrix `a`
	allocate(a(3, 3))
	pivot = [1, 2, 3]
	a = 0

	call lu_factor(a, pivot, iostat = io)
	call test_ne(io, 0, nfail, "lu_factor singular matrix")

	!********
	! Non-square `a`
	deallocate(a)
	allocate(a(3, 2))
	call lu_factor(a, pivot, iostat = io)
	call test_ne(io, 0, nfail, "lu_factor non-square matrix")

	end block b_lu_factor_c64
	!-----------------------------------

	b_cholesky_factor: block
	double precision, allocatable :: a(:,:)

	!********
	! Singular matrix `a`
	allocate(a(3, 3))
	a = 0

	call cholesky_factor(a, iostat = io)
	call test_ne(io, 0, nfail, "cholesky_factor singular matrix")

	!********
	! Non-square `a`
	deallocate(a)
	allocate(a(3, 2))
	call cholesky_factor(a, iostat = io)
	call test_ne(io, 0, nfail, "cholesky_factor non-square matrix")

	end block b_cholesky_factor
	!-----------------------------------

	b_tridiag_factor: block
	double precision, allocatable :: a(:,:)

	!********
	! Non-tridiagonal a
	allocate(a(4, 5))
	a = 0

	call tridiag_factor(a, iostat = io)
	call test_ne(io, 0, nfail, "tridiag_factor non-tridiagonal")

	!********
	! Matrix a with non-zeros in bad places
	deallocate(a)
	allocate(a(3, 5))
	a = 1
	call tridiag_factor(a, iostat = io)
	call test_ne(io, 0, nfail, "tridiag_factor with bad non-zeros")

	!********
	! Singular matrix
	deallocate(a)
	allocate(a(3, 5))
	a = 0
	call tridiag_factor(a, iostat = io)
	call test_ne(io, 0, nfail, "tridiag_factor with bad non-zeros")

	end block b_tridiag_factor
	!-----------------------------------

end function test_bad_numa_usage

!===============================================================================

integer function chapter_4_lu() result(nfail)

	character(len = *), parameter :: label = "chapter_4_lu"

	double precision :: t0, t, t_lu
	double precision, allocatable :: a(:,:), a0(:,:), bx(:), x(:), kr(:), &
		aexp(:,:)
	integer :: n, irep, p0
	integer, allocatable :: pivot(:)

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Test a fixed case
	allocate(a(5, 5))

	a(:,1) = [-4.130000d+00, 6.000000d+00, 1.100000d+01, 1.600000d+01, 2.100000d+01]
	a(:,2) = [ 2.000000d+00, 7.960000d+00, 1.200000d+01, 1.700000d+01, 2.200000d+01]
	a(:,3) = [ 3.000000d+00, 8.000000d+00, 1.280000d+00, 1.800000d+01, 2.300000d+01]
	a(:,4) = [ 4.000000d+00, 9.000000d+00, 1.400000d+01, 3.400000d+00, 2.400000d+01]
	a(:,5) = [ 5.000000d+00, 1.000000d+01, 1.500000d+01, 2.000000d+01, 5.600000d+00]

	! TODO: fill the above matrix in with some zeros to compare tridiag/banded solvers

	x = [+4.0000d+00, +2.0000d+00, -1.0000d+00, +7.0000d+00, +1.8000d+01]

	!! Print this to check that the solution is correct in the transpose sense
	!print *, "a * x = ", matmul(a, x)

	bx = [102.48d0, 274.92d0, 434.72d0, 463.8d0, 373.8d0]

	!print *, "a = "
	!print "(5es18.6)", a

	call lu_invmul(a, bx)
	!print "(a,5es18.6)", "bx = ", bx
	call test(norm2(bx - x), 0.d0, 1.d-11, nfail, "lu_invmul() 5x5")

	!********

	a = reshape([ &
		0.,   0.,   0.,   1.,     1.,     0., &
		0.,   0.,   0.,   1.,     0.25,   0., &
		1.,   0.,   0.,   1.,    -1.,     0., &
		0.,   1.,   0.,  -0.25,  -1.,     0., &
		0.,   0.,   0.,   1.,     1.,    -1., &
		0.,   0.,   1.,  -1.,     1.,     0.  &
		], &
		[6, 6] &
	)
	x = [69.d0, 420.d0, 1337.d0, 690.d0, 4200.d0, 13370.d0]
	bx = matmul(a, x)
	bx = invmul(a, bx)
	call test(norm2(bx - x), 0.d0, 1.d-14, nfail, "lu_factor 4")

	!********
	! Test a pivoting bug that went undetected for a surprisingly long time
	! until it was found in revised simplex

	a = reshape([ &
		0.0000d+00,   2.0000d+00,   5.0000d+00, &
		0.0000d+00,   3.0000d+00,   2.5000d+00, &
		1.0000d+00,   8.0000d+00,   1.0000d+01  &
		], &
		[3, 3] &
	)
	a0 = a
	pivot = range_i32(3)

	aexp = transpose(reshape([ &
		0.0000d+00,    0.0000d+00,    1.0000d+00, &
		4.0000d-01,    2.0000d+00,    4.0000d+00, &
		5.0000d+00,    2.5000d+00,    1.0000d+01  &
		], &
		[3, 3] &
	))

	call lu_factor(a, pivot)

	!call print_mat(a, "lu(a) = ")
	!print *, "pivot = ", pivot

	call test(norm2(a - aexp), 0.d0, 1.d-11, nfail, "lu_factor 1")
	call test(norm2(pivot - 1.d0*[3, 2, 1]), 0.d0, 1.d-11, nfail, "lu_factor 2")

	a = a0
	x = [69.d0, 420.d0, 1337.d0]
	bx = matmul(a, x)
	!print *, "bx = ", bx
	bx = invmul(a, bx)
	!print *, "bx = ", bx
	!print *, "x  = ", x

	call test(norm2(bx - x), 0.d0, 1.d-14, nfail, "lu_factor 3")

	!********
	! Do some fuzz testing with random data.  Chances are almost 0 for randomly
	! creating a singular matrix
	!
	! TODO: add fuzz testing for some other problems, e.g. tridiagonal and
	! banded solvers.  We will need special matmul() implementations to verify
	! results

	deallocate(a, x, bx)

	call rand_seed_determ()

	t_lu = 0.d0
	do n = 2, 90, 2

		! LU decomposition is still fast for n >> 90, but maybe not fast enough
		! for multiple reps in a unit test

		allocate(a(n, n), x(n))

		if (mod(n, 10) == 0) then
			print *, "Testing lu_invmul() with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 6
			call random_number(a)  ! random matrix
			call random_number(x)  ! random expected answer

			bx = matmul(a, x)      ! calculate rhs

			!print *, "a = "
			!print "(6es15.5)", a
			!print "(a,6es15.5)", "x  = ", x
			!print "(a,6es15.5)", "b  = ", bx

			call cpu_time(t0)
			call lu_invmul(a, bx)
			call cpu_time(t)
			t_lu = t_lu + t - t0
			!print "(a,6es15.5)", "bx = ", bx

			! Rounding error grows with `n`
			call test(norm2(bx - x), 0.d0, 1.d-9 * n, nfail, "lu_invmul() fuzz n x n")

		end do

		! Fuzz test with a randomly permuted diagonal matrix.  My pivoting bug
		! went undetected because my existing fuzz-tests have *too many*
		! non-zeros. A more spares matrix poses more of a pivoting challenge

		call random_number(x)
		a = diag(x)
		a = a(:, rand_perm(n))
		!call print_mat(a, "a = ")

		call random_number(x)
		bx = matmul(a, x)      ! calculate rhs

		call lu_invmul(a, bx)
		call test(norm2(bx - x), 0.d0, 1.d-9 * n, nfail, "lu_invmul() fuzz n x n")
		!print "(a,6es15.5)", "bx = ", bx

		deallocate(a, x, bx)
	end do
	write(*,*) "LU fuzz time = ", to_str(t_lu), " s"
	write(*,*)

	!********
	! Test lu_kernel()

	p0 = -1
	do n = 2, 40, 5

		! LU decomposition is still fast for n >> 90, but maybe not fast enough
		! for multiple reps in a unit test

		allocate(a(n, n))

		if (n/10 > p0) then
			p0 = n/10
			print *, "Testing lu_kernel() with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 1
			call random_number(a)  ! random matrix

			! Make the last column a linear combination of the other columns
			a(:,n) = matmul(a(:, :n-1), a(:n-1, n))
			!print *, "a = "
			!print "("//to_str(n)//"es15.5)", a

			a0 = a

			kr = lu_kernel(a)
			!print *, "kernel = ", kernel
			!print *, "a * kernel = ", matmul(a0, kernel)

			call test(norm2(matmul(a0, kr)), 0.d0, 1.d-9 * n, nfail, "lu_kernel()")

		end do

		deallocate(a)
	end do

	!********
	print *, ""

end function chapter_4_lu

!===============================================================================

integer function chapter_4_modify() result(nfail)

	character(len = *), parameter :: label = "chapter_4_modify"

	double precision :: aij, aij0
	double precision, allocatable :: a(:,:), a0(:,:), bx(:), x(:), &
		u(:), b(:), b0(:)
	integer :: n, irep, row, col
	integer, allocatable :: pivot(:)

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Test a fixed case
	n = 5
	allocate(a(n, n))

	a(:,1) = [-4.130000d+00, 6.000000d+00, 1.100000d+01, 1.600000d+01, 2.100000d+01]
	a(:,2) = [ 2.000000d+00, 7.960000d+00, 1.200000d+01, 1.700000d+01, 2.200000d+01]
	a(:,3) = [ 3.000000d+00, 8.000000d+00, 1.280000d+00, 1.800000d+01, 2.300000d+01]
	a(:,4) = [ 4.000000d+00, 9.000000d+00, 1.400000d+01, 3.400000d+00, 2.400000d+01]
	a(:,5) = [ 5.000000d+00, 1.000000d+01, 1.500000d+01, 2.000000d+01, 5.600000d+00]
	a0 = a

	x = [+4.0000d+00, +2.0000d+00, -1.0000d+00, +7.0000d+00, +1.8000d+01]

	bx = [102.48d0, 274.92d0, 434.72d0, 463.8d0, 373.8d0]
	b = bx
	b0 = b

	!print *, "a = "
	!print "(5es18.6)", a

	pivot = range_i32(size(a,1))
	call lu_factor(a, pivot)
	call lu_solve(a, bx, pivot)

	!print "(a,5es18.6)", "bx = ", bx
	call test(norm2(bx - x), 0.d0, 1.d-11, nfail, "lu_invmul() 5x5")

	! Test a rank-1 update.  Change the element of `a` at `row` and `col` to
	! value `aij`.  Now with the modified version of `a`, solve a*x == b
	! without refactoring `a` from scratch
	row = 3
	col = 2
	aij0 = a0(row, col)
	aij = 69.d0
	a0(row, col) = aij
	print *, "a["//to_str(row)//","//to_str(col)//"]: ", aij0, " -> ", aij

	! Use the Sherman-Morrison-Woodbury formula
	!
	! Source:  https://leobouts.medium.com/rank-1-updates-in-linear-algebra-in-simple-terms-7b5032a11a3e
	!
	! Should this be a subroutine?  Idk what the interface should be

	! Construct u and v such that modified a = a - u*v' (or plus?)
	u = zeros(n)
	u(row) = aij0 - aij

	! Solve a*z == u for z, overwrite to u
	call lu_solve(a, u, pivot)

	! Solve a*y == b for y, overwrite to b
	call lu_solve(a, b, pivot)

	! Compute x = y + ((v'*y) / (1 - v'z)) * z
	x = b + (b(col) / (1 - u(col))) * u

	call test(norm2(matmul(a0, x) - b0), 0.d0, 1.d-11, nfail, "Sherman-Morrison 5x5")

	deallocate(a, x, bx)

	!********
	! Fuzz

	call rand_seed_determ()

	do n = 2, 30, 3

		! LU decomposition is still fast for n >> 90, but maybe not fast enough
		! for multiple reps in a unit test

		allocate(a(n, n), x(n))

		!if (mod(n, 10) == 0) then
		!	print *, "Testing Sherman-Morrison with n = " // to_str(n) // " ..."
		!end if

		do irep = 1, 1
			call random_number(a)  ! random matrix
			call random_number(x)  ! random expected answer (not used for Sherman-Morrison)
			a0 = a

			b = matmul(a, x)      ! calculate rhs
			b0 = b

			!********

			pivot = range_i32(size(a,1))
			call lu_factor(a, pivot)

			row = ceiling(rand_f64() * n)
			col = ceiling(rand_f64() * n)
			aij0 = a0(row, col)
			call random_number(aij)
			a0(row, col) = aij
			print *, "a["//to_str(row)//","//to_str(col)//"]: ", &
				aij0, " -> ", aij, "(n = "//to_str(n)//")"

			u = zeros(n)
			u(row) = aij0 - aij
			call lu_solve(a, u, pivot)
			call lu_solve(a, b, pivot)
			x = b + (b(col) / (1 - u(col))) * u

			call test(norm2(matmul(a0, x) - b0), 0.d0, 1.d-11, nfail, "Sherman-Morrison nxn")

		end do

		deallocate(a, x)
	end do

	!********
	print *, ""

end function chapter_4_modify

!===============================================================================

integer function chapter_4_lu_c64() result(nfail)

	character(len = *), parameter :: label = "chapter_4_lu_c64"

	double precision, allocatable :: ar(:,:), ai(:,:), xr(:), xi(:)
	double complex, allocatable :: a0(:,:), a(:,:), bx(:), x(:), kernel(:)
	integer :: n, irep, p0

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Fuzz test

	call rand_seed_determ()

	p0 = -1
	do n = 2, 30, 3

		allocate(ar(n, n), ai(n, n), xr(n), xi(n))

		if (n/10 > p0) then
			p0 = n/10
			print *, "Testing invmul_c64() with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 2

			call random_number(ar)  ! random matrix
			call random_number(ai)  ! random matrix
			call random_number(xr)
			call random_number(xi)
			a0 = ar + IMAG_ * ai
			x  = xr + IMAG_ * xi

			a = a0
			bx = matmul(a, x)      ! calculate rhs

			!print *, "a = "
			!print "(6es15.5)", a
			!print "(a,6es15.5)", "x  = ", x
			!print "(a,6es15.5)", "b  = ", bx

			bx = invmul(a, bx)
			!print "(a,6es15.5)", "bx = ", bx

			! Rounding error grows with `n`
			call test(norm2c(bx - x), 0.d0, 1.d-12 * n, nfail, "invmul_c64() fuzz n x n")

		end do

		deallocate(ar, ai, xr, xi)
	end do

	!********
	! Test lu_kernel()

	p0 = -1
	do n = 2, 30, 5

		! LU decomposition is still fast for n >> 90, but maybe not fast enough
		! for multiple reps in a unit test

		allocate(ar(n, n), ai(n, n))

		if (n/10 > p0) then
			p0 = n/10
			print *, "Testing lu_kernel() with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 2
			call random_number(ar)  ! random matrix
			call random_number(ai)  ! random matrix
			a0 = ar + IMAG_ * ai
			a = a0

			! Make the last column a (random) linear combination of the other columns
			a(:,n) = matmul(a(:, :n-1), a(:n-1, n))

			a0 = a

			kernel = lu_kernel(a)
			!print *, "kernel = ", kernel
			!print *, "a * kernel = ", matmul(a0, kernel)

			call test(norm2c(matmul(a0, kernel)), 0.d0, 1.d-12 * n, nfail, "lu_kernel_c64()")

		end do

		deallocate(ar, ai)
	end do

	!********
	print *, ""

end function chapter_4_lu_c64

!===============================================================================

integer function chapter_4_inv() result(nfail)

	character(len = *), parameter :: label = "chapter_4_inv"

	double precision :: t0, t
	double precision :: t_gj
	double precision, allocatable :: a0(:,:), a(:,:), bx(:,:), x(:,:), ainv(:,:)

	integer :: n, irep

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Test an LU solve with multiple right hand sides

	n = 5
	allocate(a(n, n))

	a(:,1) = [-4.130000d+00, 6.000000d+00, 1.100000d+01, 1.600000d+01, 2.100000d+01]
	a(:,2) = [ 2.000000d+00, 7.960000d+00, 1.200000d+01, 1.700000d+01, 2.200000d+01]
	a(:,3) = [ 3.000000d+00, 8.000000d+00, 1.280000d+00, 1.800000d+01, 2.300000d+01]
	a(:,4) = [ 4.000000d+00, 9.000000d+00, 1.400000d+01, 3.400000d+00, 2.400000d+01]
	a(:,5) = [ 5.000000d+00, 1.000000d+01, 1.500000d+01, 2.000000d+01, 5.600000d+00]

	! Copy backup
	a0 = a

	allocate(x(n, 2))
	x(:,1) = [+4.0000d+00, +2.0000d+00, -1.0000d+00,  +7.0000d+00, +1.8000d+01]
	x(:,2) = [+1.0000d+00, +3.0000d+00, -2.0000d+00, +11.0000d+00, +1.7000d+01]

	bx = matmul(a, x)

	!print *, "a * x = "
	!print "(5es18.6)", bx

	!print *, "a = "
	!print "(5es18.6)", a

	call lu_invmul(a, bx)
	!print *, "bx = "
	!print "(5es18.6)", bx
	call test(norm2(bx - x), 0.d0, 1.d-11, nfail, "lu_invmul_mat() 5x5")

	deallocate(a, x, bx)

	!********
	! Test whole matrix inversion.  This is more for convenience than serious
	! use, as you would usually want to avoid explicit inversions

	! Reset
	a = a0

	ainv = eye(n)

	!print *, "eye = "
	!print "(5es18.6)", ainv

	call lu_invmul(a, ainv)

	!print *, "ainv = "
	!print "(5es18.6)", ainv

	!print *, "a * ainv = "
	!print "(5es18.6)", matmul(a0, ainv)

	call test(norm2(matmul(a0, ainv) - eye(n)), 0.d0, 1.d-11, nfail, "matrix inverse 1")

	!********
	! `inv` is a fn that returns the inverse of `a` without modifying it, while
	! `invert` is a subroutine that replaces `a` with its inverse
	!
	! The fn and the subroutine cannot be overloaded.  In Fortran, overloads
	! must be either all fns or all subroutines

	a = a0
	ainv = inv(a)
	call test(norm2(matmul(a, ainv) - eye(n)), 0.d0, 1.d-11, nfail, "matrix inverse 2")

	a = a0
	call invert(a)
	call test(norm2(matmul(a0, a) - eye(n)), 0.d0, 1.d-11, nfail, "matrix inverse 3")

	!********
	! Test the Gauss-Jordan method directly

	a = a0
	call gauss_jordan(a)
	call test(norm2(matmul(a0, a) - eye(n)), 0.d0, 1.d-11, nfail, "matrix inverse 4")

	deallocate(a0)

	!********
	! Fuzz test

	call rand_seed_determ()

	t_gj = 0.d0

	do n = 2, 90, 4

		allocate(a0(n, n))

		if (mod(n, 10) == 0) then
			print *, "Testing Gauss-Jordan with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 3
			call random_number(a0)  ! random matrix

			!print *, "a0 = "
			!print "(6es15.5)", a0

			a = a0
			call cpu_time(t0)
			call gauss_jordan(a)
			call cpu_time(t)
			t_gj = t_gj + t - t0
			call test(norm2(matmul(a0, a) - eye(n)), 0.d0, 1.d-9 * n, nfail, "gauss_jordan() fuzz n x n")

		end do

		deallocate(a0)
	end do

	write(*,*) "Gauss-Jordan time = ", to_str(t_gj), " s"

	!********
	print *, ""

end function chapter_4_inv

!===============================================================================

integer function chapter_4_cholesky() result(nfail)

	character(len = *), parameter :: label = "chapter_4_cholesky"

	double precision, allocatable :: a0(:,:), a(:,:), bx(:), x(:), &
		tmp_mat(:,:), lmat(:,:)

	integer :: i, j, n, irep

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Fuzz test

	call rand_seed_determ()

	do n = 2, 90, 4

		allocate(tmp_mat(n, n))
		allocate(x(n))

		if (mod(n, 10) == 0) then
			print *, "Testing Cholesky with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 3

			call random_number(tmp_mat)  ! random matrix
			call random_number(x)  ! random matrix

			! Use tmp_mat to make a symmetric (and probably positive definite)
			! matrix a0
			a0 = matmul(transpose(tmp_mat), tmp_mat)
			!print *, "a0 = "
			!print "(5es15.5)", a0

			! Calculate RHS `bx` from random expected solution `x`
			bx = matmul(a0, x)

			a = a0
			call cholesky_factor(a)
			call cholesky_solve(a, bx)

			! After cholesky_factor(), the upper triangular (or lower triangular
			! in the default Fortran print ordering) part of `a` is the same as
			! it originally was.  The diagonal and lower triangle is the
			! Cholesky factorization

			!print *, "cholesky a = "
			!print "(5es15.5)", a

			lmat = a
			do i = 1, n
			do j = 1, i-1
				lmat(j,i) = 0.d0
			end do
			end do

			!print *, "lmat = "
			!print "(5es15.5)", lmat

			!print *, "lmat * lmat' = "
			!print "(5es15.5)", matmul(lmat, transpose(lmat))

			call test(norm2(matmul(lmat, transpose(lmat)) - a0), 0.d0, 1.d-9 * n, nfail, "cholesky_factor() fuzz n x n")

			!print "(a,5es15.5)", "x  = ", x
			!print "(a,5es15.5)", "bx = ", bx

			! Cholesky requires a bigger tolerance here than LU, probably
			! because it doesn't pivot
			call test(norm2(bx - x), 0.d0, 1.d-6 * n, nfail, "cholesky_solve() fuzz n x n")

			! TODO: add cholesky_invmul(), at least for single RHS

		end do

		deallocate(tmp_mat, x)
	end do

	!********
	print *, ""

end function chapter_4_cholesky

!===============================================================================

integer function chapter_4_qr() result(nfail)

	character(len = *), parameter :: label = "chapter_4_qr"

	double precision, allocatable :: a0(:,:), a(:,:), q(:,:), r(:,:), diag_(:)

	integer :: n, irep, p0

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Fuzz test

	call rand_seed_determ()

	p0 = -1
	do n = 2, 70, 5

		allocate(a0(n, n))

		if (n/10 > p0) then
			p0 = n/10
			print *, "Testing QR decomposition with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 3

			call random_number(a0)  ! random matrix
			!print *, "a0 = "
			!print "(5es15.5)", a0

			a = a0
			call qr_factor(a, diag_)

			!print *, "qr(a) = "
			!print "(5es15.5)", a

			r = triu(a)
			!print *, "r = "
			!print "(5es15.5)", r

			q = qr_get_q_expl(a, diag_)
			!print *, "q = "
			!print "(5es15.5)", q

			call test(norm2(matmul(q, r) - a0), 0.d0, 1.d-12 * n, nfail, "qr_factor() 1 fuzz n x n")
			call test(norm2(matmul(transpose(q), q) - eye(n)), 0.d0, 1.d-12 * n, nfail, "qr_factor() 1 q' * q, n x n")

		end do

		deallocate(a0)
	end do

	!********
	print *, ""

end function chapter_4_qr

!===============================================================================

integer function chapter_4_gram_schmidt() result(nfail)

	character(len = *), parameter :: label = "chapter_4_gram_schmidt"

	double precision, allocatable :: a0(:,:), a(:,:), r(:,:)

	integer :: n, irep, p0

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Fuzz test

	call rand_seed_determ()

	p0 = -1
	do n = 2, 25, 5

		allocate(a0(n, n))

		if (n/10 > p0) then
			p0 = n/10
			print *, "Testing Gram-Schmidt QR with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 1

			call random_number(a0)  ! random matrix
			!print *, "a0 = "
			!print "(5es15.5)", a0

			a = a0
			call qr_factor_gram_schmidt(a, r)

			!print *, "gram schmidt q = "
			!print "("//to_str(n)//"es15.5)", a
			!print *, "r = "
			!print "("//to_str(n)//"es15.5)", r

			call test(norm2(matmul(transpose(a), a) - eye(n)), 0.d0, 1.d-12 * n, nfail, "gram schmidt 1 q' * q, n x n")
			call test(norm2(matmul(a, r) - a0), 0.d0, 1.d-12 * n, nfail, "gram schmidt 1 fuzz n x n")

		end do

		deallocate(a0)
	end do

	!********
	print *, ""

end function chapter_4_gram_schmidt

!===============================================================================

integer function chapter_4_lls() result(nfail)
	! Linear least squares applied to polynomial data fitting

	character(len = *), parameter :: label = "chapter_4_lls"

	double precision :: xmin, xmax
	double precision, allocatable :: x(:), y(:), p(:), pk(:), xi(:), yi(:)
	integer :: i, nx, ni, fid

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0
	call rand_seed_determ()

	!********
	! polyfit_lu().  should be avoided in favor of polyfit()

	nx = 100
	allocate(x(nx), y(nx))

	call random_number(x)
	x = x * 5 - 2.5

	call random_number(y)
	y = 2*(y - 0.5d0) * 0.1d0

	pk = [3, 5, -2, 2]  ! known polynomial coefficients

	y = y + polyval(pk, x)

	p = polyfit_lu(x, y, size(pk) - 1)
	print *, "p = ", p

	call test(norm2(p - pk), 0.d0, 1.d-1, nfail, "polyfit_lu")

	deallocate(x, y)

	!********
	! polyfit()

	nx = 100
	allocate(x(nx), y(nx))

	call random_number(x)
	x = x * 5 - 2.5

	call random_number(y)
	y = 2*(y - 0.5d0) * 0.1d0

	pk = [3, 5, -2, 2]  ! known polynomial coefficients

	y = y + polyval(pk, x)

	p = polyfit(x, y, size(pk) - 1)

	print *, "p = ", p

	! Number of interpolation points
	ni = 100
	xmin = minval(x)
	xmax = maxval(x)
	xi = range_f64(xmin, xmax, ni)
	yi = polyval(p, xi)

	open(file = "plot-poly-1.txt", newunit = fid)
	write(fid, *) "# x, y"
	write(fid, "(2es18.6)") [(xi(i), yi(i), i = 1, size(xi))]
	close(fid)

	open(file = "plot-poly-data-1.txt", newunit = fid)
	write(fid, *) "# x, y"
	write(fid, "(2es18.6)") [(x(i), y(i), i = 1, size(x))]
	close(fid)

	!print *, "diff = ", norm2(p - pk)
	call test(norm2(p - pk), 0.d0, 1.d-1, nfail, "polyfit")

	!********
	print *, ""

end function chapter_4_lls

!===============================================================================

integer function chapter_4_gauss_newton() result(nfail)
	! Nonlinear least squares applied to data fitting

	character(len = *), parameter :: label = "chapter_4_gauss_newton"

	double precision :: xmin, xmax
	double precision, allocatable :: x(:), y(:), p(:), pk(:), xi(:), yi(:)
	integer :: i, ni, fid

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********

	! Data to be fit
	x = [   0.038,   0.194,   0.425,   0.626,    1.253,    2.5,      3.74  ]
	y = [   0.05,    0.127,   0.094,   0.2122,   0.2729,   0.2665,   0.3317 ]

	! Expected (known) parameters (beta)
	pk = [0.362, 0.556]

	p = gauss_newton(x, y, rate_fn, rate_dx_fn, beta0 = [0.9d0, 0.2d0], iters = 5)
	print *, "p = ", p

	! Number of interpolation points
	ni = 100
	xmin = minval(x)
	xmax = maxval(x)
	allocate(xi(ni), yi(ni))
	xi = range_f64(xmin, xmax, ni)
	do i = 1, ni
		!yi(i) = rate_fn(xi(i), pk)
		yi(i) = rate_fn(xi(i), p)
	end do

	open(file = "plot-gn-1.txt", newunit = fid)
	write(fid, *) "# x, y"
	write(fid, "(2es18.6)") [(xi(i), yi(i), i = 1, size(xi))]
	close(fid)

	open(file = "plot-gn-data-1.txt", newunit = fid)
	write(fid, *) "# x, y"
	write(fid, "(2es18.6)") [(x(i), y(i), i = 1, size(x))]
	close(fid)

	call test(norm2(p - pk), 0.d0, 1.d-3, nfail, "gauss_newton")

	!********
	print *, ""

end function chapter_4_gauss_newton

!===============================================================================

integer function chapter_4_nelder_mead() result(nfail)
	! Nonlinear least squares applied to data fitting

	character(len = *), parameter :: label = "chapter_4_nelder_mead"

	double precision :: xmin, xmax
	double precision, allocatable :: x(:), y(:), p(:), pk(:)!, xi(:), yi(:)
	integer :: nx

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Do some data fitting, as in the gauss_newton() test driver

	! Data to be fit
	x = [   0.038,   0.194,   0.425,   0.626,    1.253,    2.5,      3.74  ]
	y = [   0.05,    0.127,   0.094,   0.2122,   0.2729,   0.2665,   0.3317 ]

	! Expected (known) parameters (beta)
	pk = [0.362, 0.556]

	p = nelder_mead_fit(x, y, rate_fn, beta0 = [100.d0, 100.d0], beta_tol = 1.d-9)!, iters = 100)
	print "(a,*(es16.6))", " p = ", p

	!! Number of interpolation points
	!ni = 100
	!xmin = minval(x)
	!xmax = maxval(x)
	!allocate(xi(ni), yi(ni))
	!xi = range_f64(xmin, xmax, ni)
	!do i = 1, ni
	!	!yi(i) = rate_fn(xi(i), pk)
	!	yi(i) = rate_fn(xi(i), p)
	!end do

	!open(file = "plot-nm-1.txt", newunit = fid)
	!write(fid, *) "# x, y"
	!write(fid, "(2es18.6)") [(xi(i), yi(i), i = 1, size(xi))]
	!close(fid)

	!open(file = "plot-nm-data-1.txt", newunit = fid)
	!write(fid, *) "# x, y"
	!write(fid, "(2es18.6)") [(x(i), y(i), i = 1, size(x))]
	!close(fid)

	call test(norm2(p - pk), 0.d0, 1.d-3, nfail, "nelder_mead_fit rate_fn")

	!********
	! Abuse the data fitting nelder_mead_fit() to find the minimum of the
	! Rosenbrock banana fn

	! Data to be fit (dummy for this rosenbrock_banana case)
	x = [69.d0]
	y = [420.d0]

	! Expected (known) parameters (beta)
	pk = [1.d0, 1.d0]

	! I implemented nelder_mead_fit() for data fitting, but it can also be used for
	! general optimization/min/max.  Here this is done with a function
	! rosenbrock_banana_beta() which has an unused dummy `x` argument.  Only the beta
	! arg is used, thus it searches for the optimal point in beta space
	p = nelder_mead_fit(x, y, rosenbrock_banana_beta, beta0 = [100.d0, 100.d0], beta_tol = 1.d-9)
	print "(a,*(es16.6))", " p = ", p

	call test(norm2(p - pk), 0.d0, 1.d-6, nfail, "nelder_mead_fit rosenbrock_banana_beta")
	deallocate(x, y)

	!********
	! Use nelder_mead(), not nelder_mead_fit(), to do optimization without the
	! extra dummy args

	! Expected (known) minimum
	pk = [1.d0, 1.d0]

	p = nelder_mead(rosenbrock_banana, [100.d0, 100.d0], x_tol = 1.d-9)
	print "(a,*(es16.6))", " p = ", p

	call test(norm2(p - pk), 0.d0, 1.d-6, nfail, "nelder_mead rosenbrock_banana 2D")

	!********
	! Optimize 4D Rosenbrock banana
	!
	! In higher dimensions, it's harder to converge for this fn if you start too far away for the minimum
	pk = [1, 1, 1, 1]
	p = nelder_mead(rosenbrock_banana_nd, [20.d0, 20.d0, 20.d0, 20.d0], x_tol = 1.d-9, iters = 5000)
	!p = nelder_mead(rosenbrock_banana_nd, [100.d0, 100.d0, 100.d0, 100.d0], x_tol = 1.d-15, iters = 10000)
	print "(a,*(es16.6))", " p = ", p
	print *, "f(p) = ", rosenbrock_banana_nd(p)
	call test(norm2(p - pk), 0.d0, 1.d-6, nfail, "nelder_mead rosenbrock_banana 4D")

	!********
	! Optimize 6D Rosenbrock banana
	pk = [1, 1, 1, 1, 1, 1]
	p = nelder_mead(rosenbrock_banana_nd, [10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0], x_tol = 1.d-9, iters = 10000)
	print "(a,*(es16.6))", " p = ", p
	print *, "f(p) = ", rosenbrock_banana_nd(p)
	call test(norm2(p - pk), 0.d0, 1.d-6, nfail, "nelder_mead rosenbrock_banana 6D")

	!********
	! Fit a bimodal distribution

	nx = 200
	xmin = -10.d0
	xmax =  10.d0
	allocate(x(nx), y(nx))
	x = range_f64(xmin, xmax, nx)

	pk = [1.4d0, 1.0d0, 2.0d0, 1.0d0, 0.2d0, -1.d0]
	y = pk(1) * exp(-pk(2) * (x - pk(3))**2) &
	  + pk(4) * exp(-pk(5) * (x - pk(6))**2)

	! This doesn't converge depending on the initial guess

	!p = nelder_mead_fit(x, y, rate_fn, beta0 = [100.d0, 100.d0], beta_tol = 1.d-9)!, iters = 100)
	!p = nelder_mead_fit(x, y, bimodal_fn, beta0 = [0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0], beta_tol = 1.d-9, iters = 10000)
	p = nelder_mead_fit(x, y, bimodal_fn, beta0 = [1.0d0, 1.0d0, 3.0d0, 1.0d0, 1.0d0, -3.0d0], beta_tol = 1.d-9)
	print "(a,*(es16.6))", " p = ", p

	call test(norm2(p - pk), 0.d0, 1.d-6, nfail, "nelder_mead_fit bimodal")

	!********
	print *, ""

end function chapter_4_nelder_mead

!===============================================================================

integer function chapter_4_qr_c64() result(nfail)
	character(len = *), parameter :: label = "chapter_4_qr_c64"

	double precision, allocatable :: ar(:,:), ai(:,:)
	double complex, allocatable :: a0(:,:), a(:,:), q(:,:), r(:,:), diag_(:)

	integer :: n, irep, p0

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********

	n = 2
	allocate(a0(n,n))
	!a0(:,1) = [dcmplx(1.d0), 2*IMAG_]
	a0(:,1) = dcmplx([(1.d0, 1.d0), 2*IMAG_])
	a0(:,2) = dcmplx([-IMAG_, dcmplx(1.d0)])

	!print *, "a0 = "
	!print "("//to_str(2*n)//"es15.5)", a0

	a = a0
	call qr_factor(a, diag_)

	!print *, "qr(a) = "
	!print "("//to_str(2*n)//"es15.5)", a

	r = triu(a)
	!print *, "r = "
	!print "("//to_str(2*n)//"es15.5)", r

	q = qr_get_q_expl(a, diag_)
	!print *, "q = "
	!print "("//to_str(2*n)//"es15.5)", q

	!print *, "q * r ="
	!print "("//to_str(2*n)//"es15.5)", matmul(q, r)
	!print *, "q' * q ="
	!print "("//to_str(2*n)//"es15.5)", matmul(transpose(conjg(q)), q)

	call test(norm2(abs(matmul(q, r) - a0)), 0.d0, 1.d-12 * n, nfail, "qr_factor_c64() q*r 2x2")
	call test(norm2(abs(matmul(transpose(conjg(q)), q) - eye(n))), 0.d0, 1.d-12 * n, nfail, "qr_factor_c64() q'*q 2x2")
	deallocate(a0)

	!--------------------------------------

	n = 3
	allocate(a0(n,n))
	a0(:,1) = dcmplx([( 1.d0,  1.d0), ( 0.d0,  2.d0), ( 0.d0,  0.d0)])
	a0(:,2) = dcmplx([( 0.d0, -1.d0), ( 1.d0,  0.d0), ( 1.d0,  0.d0)])
	a0(:,3) = dcmplx([( 0.d0,  3.d0), ( 0.d0,  0.d0), ( 2.d0,  0.d0)])
	!print *, "a0 = "
	!print "("//to_str(2*n)//"es15.5)", a0

	a = a0
	call qr_factor(a, diag_)
	!stop

	!print *, "qr(a) = "
	!print "("//to_str(2*n)//"es15.5)", a

	r = triu(a)
	!print *, "r = "
	!print "("//to_str(2*n)//"es15.5)", r

	q = qr_get_q_expl(a, diag_)
	!print *, "q = "
	!print "("//to_str(2*n)//"es15.5)", q

	!print *, "q * r ="
	!print "("//to_str(2*n)//"es15.5)", matmul(q, r)
	!print *, "q' * q ="
	!print "("//to_str(2*n)//"es15.5)", matmul(transpose(conjg(q)), q)
	!print *, "q' * q  - eye ="
	!print "("//to_str(2*n)//"es15.5)", matmul(transpose(conjg(q)), q) - eye(n)
	!print *, "norm(q' * q  - eye) =", norm2(abs(matmul(transpose(conjg(q)), q) - eye(n)))

	call test(norm2(abs(matmul(q, r) - a0)), 0.d0, 1.d-12 * n, nfail, "qr_factor_c64() q*r 3x3")
	call test(norm2(abs(matmul(transpose(conjg(q)), q) - dcmplx(eye(n)))), 0.d0, 1.d-12 * n, nfail, "qr_factor_c64() q'*q 3x3")
	deallocate(a0)

	!********
	! Fuzz test

	call rand_seed_determ()

	p0 = -1
	do n = 2, 70, 3

		allocate(ai(n,n), ar(n,n))

		if (n/10 > p0) then
			p0 = n/10
			print *, "Testing complex QR decomposition with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 3

			call random_number(ar)  ! random matrix
			call random_number(ai)  ! random matrix
			!call print_mat(ar, "ar = ")
			a0 = ar + IMAG_ * ai

			!print *, "a0 = "
			!print "("//to_str(2*n)//"es15.5)", a0

			a = a0
			call qr_factor(a, diag_)

			!print *, "qr(a) = "
			!print "("//to_str(2*n)//"es15.5)", a

			r = triu(a)
			!print *, "r = "
			!print "("//to_str(2*n)//"es15.5)", r

			q = qr_get_q_expl(a, diag_)
			!print *, "q = "
			!print "("//to_str(2*n)//"es15.5)", q

			!print *, "q * r ="
			!print "("//to_str(2*n)//"es15.5)", matmul(q, r)
			!print *, "q' * q ="
			!print "("//to_str(2*n)//"es15.5)", matmul(transpose(conjg(q)), q)
			!print *, "q' * q  - eye ="
			!print "("//to_str(2*n)//"es15.5)", matmul(transpose(conjg(q)), q) - eye(n)
			!print *, "norm(q' * q  - eye) =", norm2(abs(matmul(transpose(conjg(q)), q) - eye(n)))

			call test(norm2(abs(matmul(q, r) - a0)), 0.d0, 1.d-12 * n, nfail, "qr_factor_c64() q*r nxn")
			call test(norm2(abs(matmul(transpose(conjg(q)), q) - dcmplx(eye(n)))), 0.d0, 1.d-12 * n, nfail, "qr_factor_c64() q'*q nxn")

		end do

		deallocate(ai, ar)
	end do

	!********
	print *, ""
	!stop

end function chapter_4_qr_c64

!===============================================================================

integer function chapter_4_rank() result(nfail)
	character(len = *), parameter :: label = "chapter_4_rank"

	double precision :: delta
	double precision, allocatable :: a(:,:), a0(:,:), kr(:,:)
	integer :: i, n, rank_, rank_expect, irep

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********

	! Test rank calculation

	a = reshape([ &
		1, 0, 0, &
		0, 1, 0, &
		0, 0, 1  &
		], &
		[3, 3] &
	)
	a0 = a
	rank_ = qr_rank(a)
	print *, "qr_rank 3 = ", rank_
	call test(1.d0 * rank_, 3.d0, 1.d-14, nfail, "qr_rank 1")
	call test(is_full_rank(a0), .true., nfail, "is_full_rank 1")

	a = reshape([ &
		1, 0, 0, &
		0, 1, 0, &
		0, 1, 0  &
		], &
		[3, 3] &
	)
	a0 = a
	rank_ = qr_rank(a)
	print *, "qr_rank 2 = ", rank_
	call test(1.d0 * rank_, 2.d0, 1.d-14, nfail, "qr_rank 2")
	call test(is_full_rank(a0), .false., nfail, "is_full_rank 2")

	a = reshape([ &
		1, 0, 0, &
		0, 1, 0, &
		0, 0, 0  &
		], &
		[3, 3] &
	)
	a0 = a
	rank_ = qr_rank(a)
	print *, "qr_rank 2 = ", rank_
	call test(1.d0 * rank_, 2.d0, 1.d-14, nfail, "qr_rank 3")
	call test(is_full_rank(a0), .false., nfail, "is_full_rank 3")

	a = reshape([ &
		1, 0, 0, &
		1, 0, 0, &
		1, 0, 0  &
		], &
		[3, 3] &
	)
	rank_ = qr_rank(a)
	print *, "qr_rank 1 = ", rank_
	call test(1.d0 * rank_, 1.d0, 1.d-14, nfail, "qr_rank 4")

	a = reshape([ &
		1, 0, 0, &
		0, 0, 0, &
		1, 0, 0  &
		], &
		[3, 3] &
	)
	rank_ = qr_rank(a)
	print *, "qr_rank 1 = ", rank_
	call test(1.d0 * rank_, 1.d0, 1.d-14, nfail, "qr_rank 5")

	!********

	! This one breaks without pivoting
	a = reshape([ &
		1, 0, 0, &
		0, 0, 0, &
		0, 1, 0  &
		], &
		[3, 3] &
	)
	rank_ = qr_rank(a)
	print *, "qr_rank 2 = ", rank_
	call test(1.d0 * rank_, 2.d0, 1.d-14, nfail, "qr_rank 6")

	!********

	a = reshape([ &
		0, 0, 0, &
		1, 0, 0, &
		0, 1, 0  &
		], &
		[3, 3] &
	)
	rank_ = qr_rank(a)
	print *, "qr_rank 2 = ", rank_
	call test(1.d0 * rank_, 2.d0, 1.d-14, nfail, "qr_rank 7")

	!********
	! Source:  https://nhigham.com/2021/05/19/what-is-a-rank-revealing-factorization/

	a = reshape([ &
		0, 1, 0, 0, &
		0, 0, 1, 0, &
		0, 0, 0, 1, &
		0, 0, 0, 0  &
		], &
		[4, 4] &
	)

	a0 = a
	rank_ = qr_rank(a)
	print *, "qr_rank 4x4 = ", rank_
	call test(1.d0 * rank_, 3.d0, 1.d-14, nfail, "qr_rank 8")

	a = a0
	kr = kernel(a)
	call print_mat(kr, "kr = ")
	delta = norm2(matmul(a0, kr))
	call test(delta, 0.d0, 1.d-9, nfail, "kernel 4x4 1")

	a = transpose(a0)
	a0 = a
	rank_ = qr_rank(a)
	print *, "qr_rank 4x4 = ", rank_
	call test(1.d0 * rank_, 3.d0, 1.d-14, nfail, "qr_rank 9")

	a = a0
	kr = kernel(a)
	!call print_mat(kr, "kr = ")
	delta = norm2(matmul(a0, kr))
	call test(delta, 0.d0, 1.d-9, nfail, "kernel 4x4 2")

	!********

	a = reshape([ &
		1, 0, 0, &
		0, 0, 0, &
		0, 0, 0  &
		], &
		[3, 3] &
	)

	a0 = a
	kr = kernel(a)
	call print_mat(kr, "kr = ")
	delta = norm2(matmul(a0, kr))
	call test(delta, 0.d0, 1.d-9, nfail, "kernel 3x3 3")

	!********

	a = reshape([ &
		0, 0, 0, &
		1, 0, 0, &
		0, 0, 0  &
		], &
		[3, 3] &
	)

	a0 = a
	kr = kernel(a)
	!call print_mat(kr, "kr = ")
	delta = norm2(matmul(a0, kr))
	call test(delta, 0.d0, 1.d-9, nfail, "kernel 3x3 4")

	!********

	a = reshape([ &
		1, 0, 0, &
		1, 0, 0, &
		0, 0, 0  &
		], &
		[3, 3] &
	)

	a0 = a
	kr = kernel(a)
	!call print_mat(kr, "kr = ")
	delta = norm2(matmul(a0, kr))
	call test(delta, 0.d0, 1.d-9, nfail, "kernel 3x3 5")

	!********

	a = reshape([ &
		1, 0, 0, &
		1, 0, 0, &
		1, 0, 0  &
		], &
		[3, 3] &
	)

	a0 = a
	kr = kernel(a)
	!call print_mat(kr, "kr = ")
	!print *, "kernel col norms = ", norm2(kr, 1)
	delta = norm2(matmul(a0, kr))
	call test(delta, 0.d0, 1.d-9, nfail, "kernel 3x3 6")

	!********

	!return  ! TODO
	!stop

	! Fuzz test
	!
	! Make a few random rows/cols in an n-by-n matrix, and then make the remaining
	! rows as a linear combination of the first rows.  Finally, randomly permute
	! the rows.  I might be loose with my row/col terminology, but it doesn't
	! matter because the row rank is equivalent to the col rank
	!
	! These require a non-zero tolerance in qr_rank()

	print *, "Fuzz testing rank ..."
	do n = 10, 20
	do irep = 1, 5
		deallocate(a)
		allocate(a(n,n))
		call random_number(a)  ! fully random a
		!call print_mat(a, "a = ")

		rank_expect = rand_i32(n)

		!rank_expect = n  ! edge cases work
		!rank_expect = 0

		!print *, "n            = ", n
		!print *, "rank_expect  = ", rank_expect

		! Keep the original rank_expect cols.  Make the rest linear combos of the
		! first rank_expect cols
		do i = rank_expect+1, n
			! Use the initial values of col i as the weights of the linear combo.
			! Could be any random set of numbers
			a(:,i) = matmul(a(:, 1:rank_expect), a(1:rank_expect, i))
		end do
		!call print_mat(a, "a = ")

		! Permute
		a = a(:, rand_perm(n))
		!a = a(rand_perm(n), :)
		!if (rank_expect == 0) call print_mat(a, "a = ")
		!call print_mat(a, "a = ")
		!call print_mat(a, "a = ", fmt = "es18.8")

		a0 = a

		! Finished setting `a`
		!********

		rank_ = qr_rank(a)
		!print *, "qr_rank fuzz = ", rank_

		call test(1.d0 * rank_, 1.d0 * rank_expect, 1.d-14, nfail, "qr_rank fuzz")

		!********
		a = a0

		kr = kernel(a)

		delta = norm2(matmul(a0, kr))
		call test(delta, 0.d0, 1.d-9, nfail, "kernel fuzz()")
		!if (delta > 1.d-9) then
		!	call print_mat(kr, "kr = ")
		!	stop
		!end if

		!print *, "kernel col norms = ", norm2(kr, 1)
		!call print_mat(matmul(a0, kr), "a * kr = ")
		!print *, "a * kr = ", matmul(a0, kr)
		!print *, "norm    = ", norm2(matmul(a0, kr))

		!print *, ""
	end do
	end do
	!********
	print *, ""

end function chapter_4_rank

!===============================================================================

integer function chapter_4_linprog() result(nfail)
	character(len = *), parameter :: label = "chapter_4_linprog"

	double precision :: fval, fexpect
	double precision, allocatable :: c(:), a(:,:), b(:), x(:), expect(:), &
		a_ub(:,:), b_ub(:), lb(:), ub(:), a_eq(:,:), b_eq(:)
	integer, allocatable :: methods(:)
	integer :: n, p, im, method

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Standard form, using linprog_simplex() directly

	c = [8, 6, 0, 0]
	a = transpose(reshape([ &
	        1.0, 1.0, -1.0, 0.0, &
	        -0.05, 0.07, 0.0, 1.0  &
		], [4, 2]))
	b = [1, 0]
	expect = [0.58333333d0, 0.41666667d0, 0.d0, 0.d0]

	print *, "a = "
	print "(4es13.3)", transpose(a)

	x = linprog_simplex(c, a, b)

	print *, "x = ", x

	call test(norm2(x - expect), 0.d0, 1.d-7, nfail, "linprog_simplex 1")

	!********

	methods = [LINPROG_SIMPLEX_METH, LINPROG_REVISED_SIMPLEX_METH]
	!methods = [LINPROG_REVISED_SIMPLEX_METH, LINPROG_SIMPLEX_METH]
	do im = 1, size(methods)
	method = methods(im)

	if (method == LINPROG_SIMPLEX_METH) then
		print *, "Testing method LINPROG_SIMPLEX_METH"
	else if (method == LINPROG_REVISED_SIMPLEX_METH) then
		print *, "Testing method LINPROG_REVISED_SIMPLEX_METH"
	end if

	!********

	c = [4.d0, 8.d0, 3.d0, 0.d0, 0.d0, 0.d0]
	a_eq = transpose(reshape([ &
		2.d0, 5.d0, 3.d0, -1.d0, 0.d0, 0.d0,   &
		3.d0, 2.5d0, 8.d0, 0.d0, -1.d0, 0.d0, &
		8.d0, 10.d0, 4.d0, 0.d0, 0.d0, -1.d0   &
		], &
		[6, 3] &
	))
	b_eq = [185.d0, 155.d0, 600.d0]

	if (allocated(a_ub)) deallocate(a_ub, b_ub)
	allocate(a_ub(0,6), b_ub(0))

	expect = [66.25d0, 0.d0, 17.5d0, 0.d0, 183.75d0, 0.d0]

	x = linprog(c, a_ub, b_ub, a_eq=a_eq, b_eq=b_eq, method=method)
	print *, "x = ", x

	! Revised simplex has slightly less accuracy here than than simplex
	call test(norm2(x - expect), 0.d0, 1.d-12, nfail, "linprog 13")
	!stop

	!********
	! Example from MATLAB docs:  https://www.mathworks.com/help/optim/ug/linprog.html

	c = [-1.d0, -1.d0/3]  ! `f` in MATLAB's linprog()
	a_ub = transpose(reshape([ &  ! `A` in MATLAB
			1.0, 1.0, &
			1.0, 0.25, &
			1.0, -1.0, &
			-0.25, -1.0, &
			-1.0, -1.0, &
			-1.0, 1.0  &
		] &
		, [2, 6] &
	))
	b_ub = [2, 1, 2, 1, -1, 2]

	expect = [2.d0/3, 4.d0/3]

	x = linprog(c, a_ub, b_ub, method=method)
	print *, "x = ", x

	call test(norm2(x - expect), 0.d0, 1.d-14, nfail, "linprog 6")
	!stop

	!********
	! Scipy example including non-default bounds

	c = [-1, 4]
	a_ub = transpose(reshape([-3, 1, 1, 2] , [2, 2]))
	b_ub = [6, 4]
	lb = [-inf(), -3.d0]
	!ub = [inf(), inf()]  ! still default
	print *, "lb = ", lb

	expect = [10, -3]

	x = linprog(c, a_ub, b_ub, lb = lb, fval = fval, method=method)

	print *, "x = ", x
	print *, "fval = ", fval

	call test(norm2(x - expect), 0.d0, 1.d-14, nfail, "linprog 3")

	!stop

	!********
	! General, using linprog()
	!
	! This is the example from the scipy docs but without any non-default bounds
	!
	! By default there is still a lower bound at 0, but no upper bound

	c = [-1, 4]
	a_ub = transpose(reshape([-3, 1, 1, 2] , [2, 2]))
	b_ub = [6, 4]

	expect = [4, 0]

	x = linprog(c, a_ub, b_ub, method=method)
	print *, "x = ", x

	call test(norm2(x - expect), 0.d0, 1.d-14, nfail, "linprog 2")

	!********
	! Example from https://github.com/khalibartan/simplex-method/issues/2

	c = [4, 1]
	a_ub = transpose(reshape([-4, -3, 1, 2] , [2, 2]))
	b_ub = [-6, 4]
	a_eq = reshape([3, 1] , [1, 2])
	b_eq = [3]

	expect = [0.4d0, 1.8d0]

	x = linprog(c, a_ub, b_ub, a_eq, b_eq, method=method)
	print *, "x = ", x

	call test(norm2(x - expect), 0.d0, 1.d-14, nfail, "linprog 4")

	!********
	! Example from https://en.wikipedia.org/wiki/Simplex_algorithm#Example

	c = [-2, -3, -4]
	a_ub = transpose(reshape([3, 2, 1,   2, 5, 3] , [3, 2]))
	b_ub = [10, 15]

	expect = [0, 0, 5]

	x = linprog(c, a_ub, b_ub, method=method)
	print *, "x = ", x

	call test(norm2(x - expect), 0.d0, 1.d-14, nfail, "linprog 5")

	!********
	! MATLAB example with an equality

	c = [-1.d0, -1.d0/3]  ! `f` in MATLAB's linprog()
	a_ub = transpose(reshape([ &  ! `A` in MATLAB
			1.0, 1.0, &
			1.0, 0.25, &
			1.0, -1.0, &
			-0.25, -1.0, &
			-1.0, -1.0, &
			-1.0, 1.0  &
		] &
		, [2, 6] &
	))
	b_ub = [2, 1, 2, 1, -1, 2]
	a_eq = reshape([1.0, 0.25] , [1, 2])
	b_eq = [0.5]

	expect = [0, 2]

	x = linprog(c, a_ub, b_ub, a_eq, b_eq, method=method)
	print *, "x = ", x

	call test(norm2(x - expect), 0.d0, 1.d-14, nfail, "linprog 7")

	!********
	! MATLAB example with bounds (all constraint types)

	c = [-1.d0, -1.d0/3]  ! `f` in MATLAB's linprog()
	a_ub = transpose(reshape([ &  ! `A` in MATLAB
			1.0, 1.0, &
			1.0, 0.25, &
			1.0, -1.0, &
			-0.25, -1.0, &
			-1.0, -1.0, &
			-1.0, 1.0  &
		] &
		, [2, 6] &
	))
	b_ub = [2, 1, 2, 1, -1, 2]
	a_eq = reshape([1.0, 0.25] , [1, 2])
	b_eq = [0.5]
	lb = [-1.0, -0.5]
	ub = [1.5, 1.25]

	expect = [0.1875d0, 1.25d0]

	x = linprog(c, a_ub, b_ub, a_eq, b_eq, lb, ub, method=method)
	print *, "x = ", x

	call test(norm2(x - expect), 0.d0, 1.d-14, nfail, "linprog 8")

	!********
	! Another example from MATLAB docs

	c = [-5, -4, -6]
	a_ub = transpose(reshape([ &  ! `A` in MATLAB
		1, -1, 1, &
		3,  2, 4,  &
		3,  2, 0   &
		] &
		, [3, 3] &
	))
	b_ub = [20, 42, 30]

	expect = [0, 15, 3]

	x = linprog(c, a_ub, b_ub, method=method)
	print *, "x = ", x

	call test(norm2(x - expect), 0.d0, 1.d-14, nfail, "linprog 9")

	!********
	! Example from scipy test suite:  https://github.com/scipy/scipy/blob/main/scipy/optimize/tests/test_linprog.py

	c = [-3, -2]
	a_ub = transpose(reshape([ &
		2, 1, &
		1, 1, &
		1, 0 &
		], &
		[2, 3] &
	))
	b_ub = [10, 8, 4]

	fexpect = -18
	expect = [2, 6]

	x = linprog(c, a_ub, b_ub, fval = fval, method=method)
	print *, "x = ", x

	call test(norm2(x - expect), 0.d0, 1.d-14, nfail, "linprog 10")
	call test(fval - fexpect, 0.d0, 1.d-14, nfail, "linprog fval 10")

	!********
	! A network flow problem with supply and demand at nodes and with costs
	! along directed edges:  https://www.princeton.edu/~rvdb/542/lectures/lec10.pdf

	c = [2, 4, 9, 11, 4, 3, 8, 7, 0, 15, 16, 18]
	n = -1
	p = 1
	a_eq = transpose(reshape([ &
	    n, n, p, 0, p, 0, 0, 0, 0, p, 0, 0, &
	    p, 0, 0, p, 0, p, 0, 0, 0, 0, 0, 0, &
	    0, 0, n, n, 0, 0, 0, 0, 0, 0, 0, 0, &
	    0, 0, 0, 0, 0, 0, p, p, 0, 0, p, 0, &
	    0, 0, 0, 0, n, n, n, 0, p, 0, 0, 0, &
	    0, 0, 0, 0, 0, 0, 0, n, n, 0, 0, p, &
	    0, 0, 0, 0, 0, 0, 0, 0, 0, n, n, n], &
		[12, 7] &
	))
	b_eq = [0, 19, -16, 33, 0, 0, -36]
	expect = [19, 0, 16, 0, 0, 0, 0, 0, 0, 3, 33, 0]
	fexpect = 755

	! Make [a_ub, b_ub] optional?  It's tricky now because you have to be very
	! careful to get size(a_ub,2) correct even if size 1 is 0.  This could be
	! dangerous because ranks and types of a_ub/b_ub match a_eq/b_eq, so it
	! could invite misuse if both are optional
	!
	! It's slightly less tricky now that there's a half-decent error message to
	! notify you of a size mismatch (although, not what the size is or should
	! be)

	deallocate(a_ub, b_ub)
	!allocate(a_ub(0,0), b_ub(0))
	allocate(a_ub(0,12), b_ub(0))

	!x = linprog(c, a_ub, b_ub, method=method)
	x = linprog(c, a_ub, b_ub, a_eq=a_eq, b_eq=b_eq, fval=fval, method=method)
	!print *, "x = ", x
	print *, "fval = ", fval

	call test(norm2(x - expect), 0.d0, 1.d-12, nfail, "linprog 11")
	call test(fval - fexpect, 0.d0, 1.d-12, nfail, "linprog fval 11")
	!stop

	!********
	! "nontrivial" test from scipy

	c = [-1, 8, 4, -6]
	a_ub = transpose(reshape([ &
		-7, -7, 6, 9, &
		1, -1, -3, 0, &
		10, -10, -7, 7, &
		6, -1, 3, 4 &
		], &
		[4, 4] &
	))
	b_ub = [-3, 6, -6, 6]
	a_eq = reshape([-10, 1, 1, -8] , [1, 4])
	b_eq = [-4]

	expect = [101.d0 / 1391, 1462.d0 / 1391, 0.d0, 752.d0 / 1391]
	fexpect = 7083.d0 / 1391

	x = linprog(c, a_ub, b_ub, a_eq, b_eq, fval=fval, method=method)
	print *, "x = ", x

	call test(norm2(x - expect), 0.d0, 1.d-14, nfail, "linprog 12")
	call test(fval - fexpect, 0.d0, 1.d-14, nfail, "linprog fval 12")

	!********

	! Scipy calls `_remove_redundancy_svd()` during presolve for this example,
	! which is not implemented here
	!
	!     "a_eq does not appear to be of full row rank. To improve performance,
	!     check the problem formulation for redundant equality constraints."
	!

	c = [2.8, 6.3, 10.8, -2.8, -6.3, -10.8]
	a_eq = transpose(reshape([ &
		-1, -1, -1, 0, 0, 0, &
		0, 0, 0, 1, 1, 1, &
		1, 0, 0, 1, 0, 0, &
		0, 1, 0, 0, 1, 0, &
		0, 0, 1, 0, 0, 1 &
		], &
		[6, 5] &
	))
	b_eq = [-0.5, 0.4, 0.3, 0.3, 0.3]

	expect = [0.3d0, 0.2d0, 0.0d0, 0.0d0, 0.1d0, 0.3d0]

	deallocate(a_ub, b_ub)
	allocate(a_ub(0,6), b_ub(0))

	! With the default tolerance, the "solution is infeasible" for this case.
	! It can be worked around with a looser tolerance, although the accuracy of
	! the solution is also worsened
	!
	! Surely there are other pathological cases where the tolerance workaround
	! does not help.  Those cases would actually require redundancy removal
	!
	! No difference between LINPROG_SIMPLEX and LINPROG_REVISED_SIMPLEX

	x = linprog(c, a_ub, b_ub, a_eq, b_eq, tol=1.d-5, method=method)
	print *, "x = ", x

	call test(norm2(x - expect), 0.d0, 1.d-6, nfail, "linprog 13.1")

	!********

	c = [-0.1, -0.07, 0.004, 0.004, 0.004, 0.004]
	a_ub = transpose(reshape([ &
		1.0, 0., 0., 0., 0., 0., &
		-1.0, 0., 0., 0., 0., 0., &
		0., -1.0, 0., 0., 0., 0., &
		0., 1.0, 0., 0., 0., 0., &
		1.0, 1.0, 0., 0., 0., 0. &
		], &
		[6, 5] &
	))
	b_ub = [3.0, 3.0, 3.0, 3.0, 20.0]
	a_eq = transpose(reshape([ &
		1.0, 0., -1., 1., -1., 1., &
		0., -1.0, -1., 1., -1., 1. &
		], &
		[6, 2] &
	))
	b_eq = [0, 0]

	expect = [0, 0, 0, 0, 0, 0]

	x = linprog(c, a_ub, b_ub, a_eq, b_eq, method=method)
	print *, "x = ", x

	call test(norm2(x - expect), 0.d0, 1.d-14, nfail, "linprog 14")

	!********
	! Scipy doc example, but `c` is modified so we properly test the
	! `n_unbounded` block in the post-solve stage

	c = [4, -1]
	a_ub = transpose(reshape([-3, 1, 1, 2] , [2, 2]))
	b_ub = [6, 4] - 0
	lb = [-inf(), -3.d0]

	expect = [-3, -3]

	x = linprog(c, a_ub, b_ub, lb = lb, method=method)
	print *, "x = ", x

	call test(norm2(x - expect), 0.d0, 1.d-14, nfail, "linprog 15")

	!********
	print *, ""
	end do  ! methods
	!********

end function chapter_4_linprog

!===============================================================================

integer function chapter_6_basic_qr() result(nfail)

	character(len = *), parameter :: label = "chapter_6_basic_qr"

	double precision :: diff
	double precision, allocatable :: a(:,:), d(:,:), s(:,:), eigvals(:), &
		expect(:)

	integer :: n, irep

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Fuzz test

	call rand_seed_determ()

	! There are fewer tests here and with smaller matrices compared to
	! qr_factor() because eig_basic_qr() is an expensive iterative algorithm,
	! involving a full QR decomposition and matmul at each iteration
	do n = 4, 15, 3

		allocate(s (n, n))

		!if (mod(n, 10) == 0) then
		!	print *, "Testing basic QR algorithm with n = " // to_str(n) // " ..."
		!end if

		do irep = 1, 1

			! Construct a random matrix `a` with known real eigenvalues

			! Known eigenvalues
			expect = zeros(n)
			call random_number(expect)

			d = diag(expect)
			call sort(expect)
			!print *, "expect  = ", expect
			print "(a,*(es15.5))", " expect  = ...", expect(n-3: n)

			call random_number(s)  ! random matrix
			a = matmul(matmul(s, d), inv(s))
			!print *, "a = "
			!print "(4es15.5)", a

			eigvals = eig_basic_qr(a, iters = 5 * n)
			print "(a,*(es15.5))", " eigvals = ...", eigvals(n-3: n)

			! There is a large tolerance here because the basic QR algorithm
			! converges slowly
			diff = norm2(eigvals - expect)
			call test(diff, 0.d0, 1.d-2 * n, nfail, "eig_basic_qr 1")
			print *, "diff = ", diff
			print *, ""

		end do

		deallocate(s)
	end do

	!********
	!print *, ""

end function chapter_6_basic_qr

!===============================================================================

integer function chapter_6_hessenberg_qr() result(nfail)

	character(len = *), parameter :: label = "chapter_6_hessenberg_qr"

	double precision :: diff, t, t0, t_kr, t_pq
	double precision, allocatable :: a(:,:), d(:,:), s(:,:), eigvals(:), &
		expect(:), eigvecs(:,:), a0(:,:)

	integer :: n, irep, iters

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Fuzz test

	call rand_seed_determ()

	t_pq = 0.d0
	t_kr = 0.d0

	do n = 5, 25, 7

		allocate(s (n, n))

		print *, "Testing Hessenberg with n = " // to_str(n) // " ..."

		do irep = 1, 1

			! Construct a random matrix `a` with known real eigenvalues

			! Known eigenvalues
			expect = zeros(n)
			call random_number(expect)

			!expect(3) = expect(2)

			d = diag(expect)
			call sort(expect)
			!print *, "expect  = ", expect
			print "(a,*(es18.8))", " expect  = ...", expect(n-3: n)

			call random_number(s)  ! random matrix
			a = matmul(matmul(s, d), inv(s))
			!print *, "a = "
			!print "(4es18.8)", a

			a0 = a

			! Benchmark results:  null space (kernel) algorithm is nearly twice
			! as fast as the Q product algorithm for iters = 5 * n**2:
			!
			!     Q product  fuzz time = 2.6070019999999996 s
			!     Null space fuzz time = 1.5295850000000000 s
			!
			! Q product costs more per iteration, because it has to update the
			! product.  Null space is faster per iteration, but it has a more
			! expensive step at the end to calculate the eigenvectors.  Null
			! space has to perform n QR decompositions to get the eigenvectors,
			! while the Q product method only has to do one LU solve on the full
			! n x n matrix (but many more LU solves on smaller matrices of size
			! n-1, n-2, n-3, ...)
			!
			! If you only do a very small number of iterations, Q product can
			! come out faster, but it doesn't converge very well.  For practical
			! numbers of iterations, the null space algo seems better
			!
			! Also, the `diff vec` printed below has a 10x smaller residual for
			! the null space algo
			!
			! Shifting algorithms should converge faster, so it's good to keep
			! the Q product algorithm for now.  It may be the better performer
			! with shifting

			iters = 5 * n**2
			!iters = n**2
			!iters = n
			!iters = 10*n

			!********
			print *, "Hessenberg QR with Q product:"

			! Hessenberg QR doesn't converge in any fewer iterations than basic
			! QR, but each iteration is only O(n**2) instead of O(n**3), so it's
			! cheap to just do a bunch of iterations
			call cpu_time(t0)
			eigvals = eig_hess_qr(a, iters, eigvecs)
			call cpu_time(t)
			t_pq = t_pq + t - t0

			!print *, "eigvecs = "
			!print "(5es15.5)", eigvecs

			! Check the eigenvalues
			diff = norm2(sorted(eigvals) - expect)
			call test(diff, 0.d0, 1.d-4 * n, nfail, "eig_hess_qr val 1")
			print *, "diff val = ", diff

			!print *, "a * eigvecs / eigvecs = "
			!print "(5es15.5)", matmul(a0, eigvecs) / eigvecs
			!!print *, "spread = "
			!!print "(4es15.5)", spread(eigvals, 1, n)

			! Check the eigenvectors: A * eigvecs = eigvals * eigvecs
			diff = norm2(matmul(a0, eigvecs) / eigvecs - spread(eigvals, 1, n))
			call test(diff, 0.d0, 1.d-3 * n, nfail, "eig_hess_qr vec 2")
			print *, "diff vec = ", diff

			print *, ""
			!********
			print *, "Hessenberg QR with kernel:"

			a = a0
			call cpu_time(t0)
			eigvals = eig_hess_qr_kernel(a, iters, eigvecs)
			call cpu_time(t)
			t_kr = t_kr + t - t0

			! The magnitudes of eigvecs may differ from above, but their
			! directions are the same, thus the spread test passes in all cases
			!print *, "eigvecs = "
			!print "(5es15.5)", eigvecs

			! Check the eigenvalues
			diff = norm2(sorted(eigvals) - expect)
			call test(diff, 0.d0, 1.d-4 * n, nfail, "eig_hess_qr_kernel val 1")
			print *, "diff val = ", diff

			!print *, "a * eigvecs / eigvecs = "
			!print "(5es15.5)", matmul(a0, eigvecs) / eigvecs
			!!print *, "spread = "
			!!print "(4es15.5)", spread(eigvals, 1, n)

			! Check the eigenvectors: A * eigvecs = eigvals * eigvecs
			diff = norm2(matmul(a0, eigvecs) / eigvecs - spread(eigvals, 1, n))
			call test(diff, 0.d0, 1.d-3 * n, nfail, "eig_hess_qr_kernel vec 2")
			print *, "diff vec = ", diff

			print *, ""

			!********
			!stop

		end do

		deallocate(s)
	end do

	write(*,*) "Q product  fuzz time = ", to_str(t_pq), " s"
	write(*,*) "Null space fuzz time = ", to_str(t_kr), " s"

	!********
	print *, ""

end function chapter_6_hessenberg_qr

!===============================================================================

double precision function diff_complex_vecs(a, b) result(delta)
	! Compare the "difference" between two complex vectors.  Cost is O(n**2) so
	! avoid this when possible
	double complex, intent(in) :: a(:), b(:)
	!********

	double complex :: ai, bi
	double precision :: di
	integer :: i

	delta = 1.d0
	if (size(a) /= size(b)) return

	delta = 0.d0

	! Apparently comparing complex vectors is hard.  It's not enough to sort
	! then take the norm of the difference, as is possible with doubles, because
	! small differences will completely change the lexicographical sort order

	! For each a, find the closest element in b, and vice versa.  Add up the
	! differences in a norm2 sense
	!
	! The two-way comparison is required in case the mapping is onto but not 1:1

	do i = 1, size(a)
		ai = a(i)
		di = minval(abs(b - ai))
		delta = delta + di**2
		!print *, "di = ", di
	end do

	do i = 1, size(b)
		bi = b(i)
		di = minval(abs(a - bi))
		delta = delta + di**2
		!print *, "di = ", di
	end do

	delta = sqrt(delta)

end function diff_complex_vecs

!===============================================================================

integer function chapter_6_francis_qr() result(nfail)

	character(len = *), parameter :: label = "chapter_6_francis_qr"

	double precision :: diff
	double precision, allocatable :: a(:,:), d(:,:), s(:,:), &
		expect(:), a0(:,:)!, ar(:,:), ai(:,:)
	double complex, allocatable :: eigvals(:), eigvals2(:), eigvecs(:,:)

	integer :: n, irep, p0

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	! Matrix from literature with complex eigenvalues

	n = 6
	allocate(a(n,n))
	a(:,1) = [ 7,  3,  4, -11, -9, -2]
	a(:,2) = [-6,  4, -5,   7,  1, 12]
	a(:,3) = [-1, -9,  2,   2,  9,  1]
	a(:,4) = [-8,  0, -1,   5,  0,  8]
	a(:,5) = [-4,  3, -5,   7,  2, 10]
	a(:,6) = [ 6,  1,  4, -11, -7, -1]
	print *, "a = "
	print "("//to_str(n)//"es19.9)", a
	a0 = a
	!    --> spec(a')
	!     ans  =
	!
	!       5. + 6.i
	!       5. - 6.i
	!       1. + 2.i
	!       1. - 2.i
	!       3. + 0.i
	!       4. + 0.i

	eigvals = eig_francis_qr(a, eigvecs)
	print *, "eigvals = [real, imag]"
	print "(2es15.5)", sorted(eigvals)

	! Eigenvectors are optional to save work.  Beware, rounding and ordering may be different
	a = a0
	eigvals2 = eig_francis_qr(a)
	diff = diff_complex_vecs(eigvals, eigvals2)
	print *, "diff eigvals2 = ", diff
	call test(diff, 0.d0, 1.d-8, nfail, "eig_francis_qr no vecs")

	!print *, "eigvals2 = [real, imag]"
	!print "(2es15.5)", sorted(eigvals2)

	!print *, "eigvecs = "
	!print "("//to_str(2*n)//"es15.5)", eigvecs

	!print *, "a * eigvecs / eigvecs = "
	!print "("//to_str(2*n)//"es15.5)", matmul(a0, eigvecs) / eigvecs

	!print *, "a * eigvecs = "
	!print "("//to_str(2*n)//"es15.5)", matmul(a0, eigvecs)

	!print *, "eigvals * eigvecs = "
	!print "("//to_str(2*n)//"es15.5)", spread(eigvals, 1, n) * eigvecs

	!print *, "a * eigvecs - eigvals * eigvecs = "
	!print "("//to_str(2*n)//"es15.5)", matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs

	diff = norm2(abs( &
		matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs))
	print *, "diff = ", diff
	call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_francis_qr vec 1")

	!! Some components of some eigenvectors are 0, so dividing like this is
	!! ill-conditioned.  Need to multiply instead
	!diff = norm2(abs(matmul(a0, eigvecs) / eigvecs - spread(eigvals, 1, n)))
	!call test(diff, 0.d0, 1.d-3 * n, nfail, "eig_francis_qr vec 2")

	!********
	! Fuzz test

	call rand_seed_determ()

	p0 = -1
	do n = 5, 35, 3

		allocate(s (n, n))

		if (n/10 > p0) then
			p0 = n/10
			print *, "Testing real eig_francis_qr() with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 2

			! Construct a random matrix `a` with known real eigenvalues

			! Known eigenvalues
			expect = zeros(n)
			call random_number(expect)

			!expect(3) = expect(2)

			d = diag(expect)
			call sort(expect)
			!print "(a,*(es18.8))", " expect  = ...", expect(n-3: n)

			call random_number(s)  ! random matrix
			a = matmul(matmul(s, d), inv(s))
			!print *, "a = "
			!print "("//to_str(n)//"es19.9)", a

			a0 = a

			!********

			eigvals = eig_francis_qr(a, eigvecs)

			!print *, "eigvecs = "
			!print "("//to_str(2*n)//"es15.5)", eigvecs

			! Check the eigenvalues

			diff = norm2(sorted(dble(eigvals)) - expect)

			call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_francis_qr val 1")
			!print *, "diff val = ", diff
			!print *, "expect = ", expect
			!print *, "actual = ", sorted(dble(eigvals))

			!********

			!print *, "a * eigvecs / eigvecs = "
			!print "("//to_str(2*n)//"es15.5)", matmul(a0, eigvecs) / eigvecs
			!!print *, "spread = "
			!!print "("//to_str(n)//"4es15.5)", spread(eigvals, 1, n)

			! Check the eigenvectors: A * eigvecs = eigvals * eigvecs
			!
			! idk if norm2(abs( [complex_matrix] )) is a good norm but i don't feel like making my
			! own matrix norm rn.  need to figure out what built-in norm2() is
			! doing on matrices, whether it just reshapes them to vec or does
			! some more involved matrix norm

			!print *, "eigvals = ", eigvals

			diff = norm2(abs( &
				matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs))
			call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_francis_qr vec 3")
			!print *, "diff vec = ", diff
			!print *, ""

		end do

		deallocate(s)
	end do
	deallocate(a)

	p0 = -1
	do n = 5, 65, 7
		! eig_francis_qr() can't converge over about n == 80

		allocate(a(n, n))

		if (n/10 > p0) then
		!if (n/1 > p0) then
			p0 = n/10
			print *, "Testing complex eig_francis_qr() with n = " // to_str(n) // " ..."
		end if

		!do irep = 1, 10
		do irep = 1, 1
			!print *, "irep = ", irep

			! Construct a random matrix `a`.  It can have complex eigenvalues,
			! most likely a mixture of real and complex eigenvalues
			call random_number(a)
			!print *, "a = "
			!print "("//to_str(n)//"es19.9)", a
			a0 = a

			!********

			eigvals = eig_francis_qr(a, eigvecs)

			!print *, "eigvals = [real, imag]"
			!print "(2es15.5)", eigvals

			!print *, "eigvecs = "
			!print "("//to_str(2*n)//"es15.5)", eigvecs

			! Eigenvalues cannot be verified by themselves because we don't know
			! what they should be a priori
			!
			! We can still verify that `a0 * eigvecs == eigvals * eigvecs`

			!print *, "a * eigvecs / eigvecs = "
			!print "("//to_str(2*n)//"es15.5)", matmul(a0, eigvecs) / eigvecs
			!!print *, "spread = "
			!!print "("//to_str(n)//"4es15.5)", spread(eigvals, 1, n)

			!print *, "a * eigvecs - eigvals * eigvecs = "
			!print "("//to_str(2*n)//"es15.5)", matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs

			!print *, "checking diff"
			!diff = norm2(abs(matmul(a0, eigvecs) / eigvecs - spread(eigvals, 1, n)))
			diff = norm2(abs( &
				matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs))
			call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_francis_qr vec 4")
			!print *, "diff vec = ", diff

			!print *, ""

		end do

		deallocate(a)
	end do
	!********
	print *, ""

end function chapter_6_francis_qr

!===============================================================================

integer function chapter_6_eig_lapack() result(nfail)

	character(len = *), parameter :: label = "chapter_6_eig_lapack"

	double precision :: diff
	double precision, allocatable :: a(:,:), d(:,:), s(:,:), &
		expect(:), a0(:,:)!, ar(:,:), ai(:,:)
	double complex, allocatable :: eigvals(:), eigvals2(:), eigvecs(:,:)

	integer :: n, irep, p0

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	! Matrix from literature with complex eigenvalues

	n = 6
	allocate(a(n,n))
	a(:,1) = [ 7,  3,  4, -11, -9, -2]
	a(:,2) = [-6,  4, -5,   7,  1, 12]
	a(:,3) = [-1, -9,  2,   2,  9,  1]
	a(:,4) = [-8,  0, -1,   5,  0,  8]
	a(:,5) = [-4,  3, -5,   7,  2, 10]
	a(:,6) = [ 6,  1,  4, -11, -7, -1]
	print *, "a = "
	print "("//to_str(n)//"es19.9)", a
	a0 = a
	!    --> spec(a')
	!     ans  =
	!
	!       5. + 6.i
	!       5. - 6.i
	!       1. + 2.i
	!       1. - 2.i
	!       3. + 0.i
	!       4. + 0.i

	eigvals = eig_lapack(a, eigvecs)
	print *, "eigvals = [real, imag]"
	print "(2es15.5)", sorted(eigvals)

	! Eigenvectors are optional to save work.  Beware, rounding and ordering may be different
	a = a0
	eigvals2 = eig_lapack(a)
	diff = diff_complex_vecs(eigvals, eigvals2)
	print *, "diff eigvals2 = ", diff
	call test(diff, 0.d0, 1.d-8, nfail, "eig_lapack no vecs")

	diff = norm2(abs( &
		matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs))
	print *, "diff = ", diff
	call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_lapack vec 1")

	!stop

	!! Some components of some eigenvectors are 0, so dividing like this is
	!! ill-conditioned.  Need to multiply instead
	!diff = norm2(abs(matmul(a0, eigvecs) / eigvecs - spread(eigvals, 1, n)))
	!call test(diff, 0.d0, 1.d-3 * n, nfail, "eig_lapack vec 2")

	!********
	! Fuzz test

	call rand_seed_determ()

	p0 = -1
	do n = 5, 25, 5

		allocate(s (n, n))

		if (n/10 > p0) then
			p0 = n/10
			print *, "Testing real eig_lapack() with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 2

			! Construct a random matrix `a` with known real eigenvalues

			! Known eigenvalues
			expect = zeros(n)
			call random_number(expect)

			d = diag(expect)
			call sort(expect)

			call random_number(s)  ! random matrix
			a = matmul(matmul(s, d), inv(s))
			!print *, "a = "
			!print "("//to_str(n)//"es19.9)", a
			a0 = a

			!********

			eigvals = eig_lapack(a, eigvecs)

			! Check the eigenvalues
			diff = norm2(sorted(dble(eigvals)) - expect)
			call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_lapack val 1")
			!********

			! Check the eigenvectors: A * eigvecs = eigvals * eigvecs
			!
			! idk if norm2(abs( [complex_matrix] )) is a good norm but i don't feel like making my
			! own matrix norm rn.  need to figure out what built-in norm2() is
			! doing on matrices, whether it just reshapes them to vec or does
			! some more involved matrix norm
			diff = norm2(abs( &
				matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs))
			call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_lapack vec 3")

		end do

		deallocate(s)
	end do
	deallocate(a)

	p0 = -1
	!do n = 5, 175, 5  ! eig_lapack() is robust even for larger matrices
	do n = 5, 95, 11

		allocate(a(n, n))

		if (n/10 > p0) then
			p0 = n/10
			print *, "Testing complex eig_lapack() with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 1
		!do irep = 1, 2

			!print *, "irep = ", irep

			! Construct a random matrix `a`.  It can have complex eigenvalues,
			! most likely a mixture of real and complex eigenvalues
			call random_number(a)
			!print *, "a = "
			!print "("//to_str(n)//"es19.9)", a
			a0 = a
			!********

			eigvals = eig_lapack(a, eigvecs)

			! Eigenvalues cannot be verified by themselves because we don't know
			! what they should be a priori
			!
			! We can still verify that `a0 * eigvecs == eigvals * eigvecs`
			diff = norm2(abs( &
				matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs))
			call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_lapack vec 4")

		end do

		deallocate(a)
	end do
	!********
	print *, ""

end function chapter_6_eig_lapack

!===============================================================================

end module numa__exercises

