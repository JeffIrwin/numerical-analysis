
!> This is a Fortran module for Mersenne twister (MT) pseudo-random number
!> generation (PRNG)
!>
!> MT is *not* cryptographically secure.  From wikipedia:
!>
!> >   [MT] is not cryptographically secure. The reason is that observing a
!> >   sufficient number of iterations (624 in the case of MT19937, since this is
!> >   the size of the state vector from which future iterations are produced)
!> >   allows one to predict all future iterations
!>
!> From https://github.com/JeffIrwin/rfng
module numa__rng

#define int32_t integer(kind = int32)
#define int64_t integer(kind = int64)

	use iso_c_binding
	use iso_fortran_env

	implicit none

	public :: rng_t
	private

	! TODO: add 64-bit version.  From wikipedia it seems that the difference is
	! just a few constants
	int32_t, parameter :: n32 = 624

	!********

	public :: rng_state_t
	type, bind(c) :: rng_state_t
		! C struct interface
		integer(kind = c_int32_t) :: mt(0: n32-1)  ! state vector
		integer(kind = c_int32_t) :: index_ = n32 + 1
	end type rng_state_t

	!********

	type :: rng_t
		! Type with bound procedures cannot bind(c), so this type just wraps the
		! C interface `rng_state_t` and adds convenience procedures for Fortran

		type(rng_state_t) :: s

		contains
			procedure :: &
				seed   => seed_mt19937, &
				int32  => int32_mt19937, &
				uint32 => uint32_mt19937

	end type rng_t

!===============================================================================

contains

!===============================================================================

function seed_array(seed_, len_) bind(c, name = "seed_array_fort") result(rng)

	! This is not public but C++ dgaf
	!
	! TODO: DRY with seed()

	integer(c_int32_t), intent(in) :: seed_(len_)
	integer(c_int32_t), intent(in) :: len_

	type(rng_state_t) :: rng

	!********

	int32_t, parameter :: &
		w = 32, &
		f = 1812433253

	int32_t :: i

	!print *, "len_ = ", len_
	!print *, "seed_(1) = ", seed_(1)
	!print *, "seed_ = ", seed_

	rng%index_ = n32

	!rng%mt(0) = seed_(1)
	rng%mt(0: min(len_, n32) - 1) = seed_(1: min(len_, n32))

	! Reference C implementation expands keys differently but YOLO
	do i = min(len_, n32), n32 - 1
		rng%mt(i)  = f * ieor(rng%mt(i-1), shiftr(rng%mt(i-1), w-2)) + i
	end do
	!print *, "rng%mt = ", rng%mt
	!print *, "rng%mt = ", rng%mt(0: 10)

end function seed_array

!===============================================================================

function seed(seed_) bind(c, name = "seed_fort") result(rng)

	integer(c_int32_t), intent(in) :: seed_
	type(rng_state_t) :: rng

	!********

	int32_t, parameter :: &
		w = 32, &
		f = 1812433253

	int32_t :: i

	!print *, "seed_ = ", seed_

	rng%index_ = n32
	rng%mt(0) = seed_

	! Honestly it's a lot easier and more compact to write bit-fiddling
	! operations like this in C instead of Fortran
	do i = 1, n32 - 1
		rng%mt(i)  = f * ieor(rng%mt(i-1), shiftr(rng%mt(i-1), w-2)) + i
	end do
	!print *, "rng%mt = ", rng%mt
	!print *, "rng%mt = ", rng%mt(0: 4)

end function seed

!===============================================================================

subroutine seed_mt19937(rng, seed_)

	! Fortran convenience interface to seed()

	class(rng_t) :: rng
	int32_t, intent(in) :: seed_

	rng%s = seed(seed_)

end subroutine seed_mt19937

!===============================================================================

function get_int32(rng) bind(c, name = "int32_fort") result(num)

	type(rng_state_t)  :: rng
	integer(c_int32_t) :: num

	!********

	int32_t, parameter :: &
		b = int(z"9d2c5680", int32), &
		c = int(z"efc60000", int32), &
		d = int(z"ffffffff", int32), &
		l = 18, &
		s =  7, &
		t = 15, &
		u = 11

	int32_t :: y

	if (rng%index_ >= n32) then
		if (rng%index_ > n32) rng = seed(5489)
		call twist_mt19937(rng)
	end if

	y = rng%mt(rng%index_)
	y = ieor(y, iand(shiftr(y, u), d))
	y = ieor(y, iand(shiftl(y, s), b))
	y = ieor(y, iand(shiftl(y, t), c))
	y = ieor(y, shiftr(y, l))
	rng%index_ = rng%index_ + 1

	num = y

end function get_int32

!===============================================================================

function get_uint32(rng) bind(c, name = "uint32_fort") result(num)

	type(rng_state_t)  :: rng
	integer(c_int64_t) :: num

	!********
	!int64_t :: num

	num = iand( &
		!int(rng%s%int32(), int64), &
		int(get_int32(rng), int64), &
		int(z"ffffffff"   , int64))

	!print *, "num = ", num

end function get_uint32

!===============================================================================

function int32_mt19937(rng) result(num)

	class(rng_t) :: rng
	int32_t :: num

	num = get_int32(rng%s)

end function int32_mt19937

!===============================================================================

function uint32_mt19937(rng) result(num)

	! Extract an unsigned 32 bit int but return a signed 64 bit int because
	! that's all that Fortran has :(

	class(rng_t) :: rng

	int64_t :: num

	num = iand( &
		int(rng%int32(), int64), &
		int(z"ffffffff", int64))

end function uint32_mt19937

!===============================================================================

subroutine twist_mt19937(rng)

	type(rng_state_t) :: rng

	!********

	int32_t, parameter :: &
		a = int(z"9908b0df", int32), &
		m = 397!, &
		!r =  31

	!! nvfortran doesn't like using shiftl() in parameter definition
	!int32_t, parameter :: &
	!	lower_mask = shiftl(1, r) - 1, &
	!	upper_mask = not(lower_mask)
	int32_t :: lower_mask, upper_mask

	int32_t :: i, x, xa

	!lower_mask = shiftl(1, r) - 1
	!lower_mask = int(int(shiftl(1, r),8) - 1, 4)
	lower_mask = int(z"7FFFFFFF")

	upper_mask = not(lower_mask)

	!print "(a,b0.32)", "lower_mask = ", lower_mask
	!print "(a,b0.32)", "upper_mask = ", upper_mask

	!! Only using error unit to debug without causing recursive IO (e.g. if
	!! caller is calling int32() within a print statement)
	!write(error_unit,*) "twist_mt19937()"

	! Fortran `mod()` is consistent with the C `%` operator. Fortran `modulo()`
	! works differently for negative args
	do i = 0, n32 - 1
		x = ior(iand(rng%mt(i)            , upper_mask), &
		        iand(rng%mt(mod(i+1, n32)), lower_mask))
		xa = shiftr(x, 1)
		if (mod(x, 2) /= 0) xa = ieor(xa, a)
		rng%mt(i) = ieor(rng%mt(mod(i + m, n32)), xa)
	end do
	rng%index_ = 0

end subroutine twist_mt19937

!===============================================================================

end module numa__rng

!!===============================================================================
!
!#ifdef RNG_FORT_TEST
!
!module rng_test_m
!
!	use numa__rng
!	implicit none
!
!	integer :: nfail_glbl = 0, npass_glbl = 0
!
!contains
!
!!********
!
!function to_str(i) result(str_)
!
!	integer, intent(in) :: i
!	character(len = :), allocatable :: str_
!
!	character :: buffer*16
!
!	write(buffer, "(i0)") i
!	str_ = trim(buffer)
!
!end function to_str
!
!!********
!
!! Might be better if we take expected vs actual values as args, then we can
!! print them instead of just having a binary pass/fail
!#define TEST_(x) call test_((x), __FILE__, __LINE__)
!subroutine test_(test, file_, line)
!
!	logical, intent(in) :: test
!
!	character(len = *), intent(in) :: file_
!
!	integer, intent(in) :: line
!
!	if (test) then
!		npass_glbl = npass_glbl + 1
!	else
!		nfail_glbl = nfail_glbl + 1
!		write(*,*) "Error: assertion failed in """//file_ &
!			//""" at line "//to_str(line)
!	end if
!
!end subroutine test_
!
!!********
!
!subroutine rng_test()
!
!	integer :: i, hist_min, hist_max, mod_, rand_
!	integer(kind = 8) :: dummy
!	!integer :: fid
!	integer, allocatable :: hist(:)
!
!	type(rng_t) :: rng, ra, rb
!
!	! To compare results with reference C implementation, call init_genrand(0)
!	! instead of init_by_array() as in the original.  The reference C
!	! implementation is here:
!	!
!	!     http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/CODES/mt19937ar.c
!	!
!	! It can be easier to print hex values instead of decimal-formatted ints
!	! because of signed vs unsigned differences with Fortran
!
!	!print *, rng%int32(), rng%int32(), rng%int32(), rng%int32(), rng%int32()
!
!	! Test default (not explicitly seeded) seed
!	TEST_(rng%uint32() == 3499211612)
!
!	! Test explicit seed
!	call rng%seed(0)
!	TEST_(rng%uint32() == 2357136044)
!	TEST_(rng%uint32() == 2546248239)
!
!	! Test re-seeding.  Should get same number as before
!	call rng%seed(0)
!	TEST_(rng%uint32() == 2357136044)
!
!	! Test > 624 calls.  This will trigger another twist_mt19937() call
!	call rng%seed(0)
!	do i = 1, 997
!		dummy = rng%uint32()
!	end do
!	TEST_(rng%uint32() == 3814118674)
!	TEST_(rng%uint32() == 2172679577)
!	TEST_(rng%uint32() == 3043451800)
!
!	! Test multiple RNGs running concurrently
!	call ra%seed(0)
!	call rb%seed(0)
!	TEST_(ra%uint32() == 2357136044)
!	TEST_(rb%uint32() == 2357136044)
!	TEST_(ra%uint32() == 2546248239)
!	TEST_(rb%uint32() == 2546248239)
!	TEST_(ra%uint32() == 3071714933)
!	TEST_(rb%uint32() == 3071714933)
!
!	! Test (signed) int32 generation
!	call rng%seed(0)
!	TEST_(rng%int32() == -1937831252)
!	TEST_(rng%int32() == -1748719057)
!	TEST_(rng%int32() == -1223252363)
!	TEST_(rng%int32() == -668873536)
!	TEST_(rng%int32() == -1706118333)
!
!	!********
!
!	write(*,*) to_str(npass_glbl)//" / "//to_str(npass_glbl + nfail_glbl) &
!		//" tests passed"
!	write(*,*) to_str(nfail_glbl)//" tests failed"
!	write(*,*)
!
!	! TODO: histogram test
!	return
!
!	do i = 1, 1000
!		!print "(a,z0.8)", "rng = ", rng%int32()
!		print *, "rng = ", rng%uint32()
!		!print "(a,z0.8)", "rng = ", rng%uint32()
!		!print *, "rng mod = ", modulo(rng%int32(), 10001)
!	end do
!
!	!********
!
!	! Check some results by making a histogram that counts the number of each
!	! random value
!	mod_ = 5
!	hist_min = 0
!	hist_max = mod_ - 1
!	allocate(hist(hist_min: hist_max))
!	hist = 0
!
!	call rng%seed(0)
!	do i = 1, 1000000
!		rand_ = modulo(rng%int32(), mod_)
!		hist(rand_) = hist(rand_) + 1
!	end do
!	write(*,*) "Histogram ="
!	write(*,*) "["
!	write(*,"(i12)") hist
!	write(*,*) "]"
!
!	!open(newunit = fid, file = "rng-out-1.csv")
!	!do i = 1, 1000
!	!	write(fid, "(i0)") modulo(rng%int32(), 10)
!	!end do
!	!close(fid)
!
!end subroutine rng_test
!
!end module rng_test_m
!
!!===============================================================================
!
!program main
!
!	use numa__rng
!	use rng_test_m
!
!	implicit none
!
!	write(*,*) "starting rng main"
!	write(*,*)
!
!	call rng_test()
!
!	write(*,*) "ending rng main"
!	write(*,*)
!
!	call exit(nfail_glbl)
!
!end program main
!
!#endif
!
!!===============================================================================
!
