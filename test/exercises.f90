
!> Unit tests and exercises
module numa__exercises

	use numa

	implicit none

	interface test
		procedure :: test_f64
		procedure :: test_i32
		procedure :: test_i64
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

subroutine test_i64(val, expect, nfail, msg)
	integer(kind = 8), intent(in) :: val, expect
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

end subroutine test_i64

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

	b_test_rng: block

	use numa__rng
	integer :: i
	integer(kind = 8) :: dummy

	type(rng_t) :: rng, ra, rb

	! To compare results with reference C implementation, call init_genrand(0)
	! instead of init_by_array() as in the original.  The reference C
	! implementation is here:
	!
	!     http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/CODES/mt19937ar.c
	!
	! It can be easier to print hex values instead of decimal-formatted ints
	! because of signed vs unsigned differences with Fortran

	!print *, rng%int32(), rng%int32(), rng%int32(), rng%int32(), rng%int32()

	! Test default (not explicitly seeded) seed
	call test(rng%uint32(), int(z"D091BB5C", 8), nfail, "rng 1")

	! Test explicit seed
	call rng%seed(0)
	call test(rng%uint32(), int(z"8C7F0AAC", 8), nfail, "rng 2")
	call test(rng%uint32(), int(z"97C4AA2F", 8), nfail, "rng 3")

	! Test re-seeding.  Should get same number as before
	call rng%seed(0)
	call test(rng%uint32(), int(z"8C7F0AAC", 8), nfail, "rng 4")

	! Test > 624 calls.  This will trigger another twist_mt19937() call
	call rng%seed(0)
	do i = 1, 997
		dummy = rng%uint32()
	end do
	call test(rng%uint32(), 3814118674_8, nfail, "rng 5")
	call test(rng%uint32(), 2172679577_8, nfail, "rng 6")
	call test(rng%uint32(), 3043451800_8, nfail, "rng 7")

	! Test multiple RNGs running concurrently
	call ra%seed(0)
	call rb%seed(0)
	call test(ra%uint32(), 2357136044_8, nfail, "rng 8.1")
	call test(rb%uint32(), 2357136044_8, nfail, "rng 8.2")
	call test(ra%uint32(), 2546248239_8, nfail, "rng 8.3")
	call test(rb%uint32(), 2546248239_8, nfail, "rng 8.4")
	call test(ra%uint32(), 3071714933_8, nfail, "rng 8.5")
	call test(rb%uint32(), 3071714933_8, nfail, "rng 8.6")

	! Test (signed) int32 generation
	call rng%seed(0)
	call test(rng%int32(), -1937831252, nfail, "rng 9.1")
	call test(rng%int32(), -1748719057, nfail, "rng 9.2")
	call test(rng%int32(), -1223252363, nfail, "rng 9.3")
	call test(rng%int32(), -668873536 , nfail, "rng 9.4")
	call test(rng%int32(), -1706118333, nfail, "rng 9.5")

	end block b_test_rng
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

end module numa__exercises

