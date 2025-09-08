
!> Module for high-level random number generation
module numa__rand

	use numa__rng
	implicit none

	interface rand_f64
		procedure :: rand_sca_f64
		procedure :: rand_vec_f64
		procedure :: rand_mat_f64
	end interface rand_f64

	type(rng_t) :: global_rng

contains

!===============================================================================

!> Seed the RNG with 0's for repeatable tests
subroutine rand_seed_determ()
	!integer :: i, nrng

	! numa__rng
	call global_rng%seed(0)

	!! Built-in Fortran default rng.  The built-in rng is excellent for most
	!! compilers but nvfortran creates random matrices with many repeated
	!! elements, which is the worst case scenario for my tests which depend on
	!! randomly-generated non-singular matrices
	!call random_seed(size = nrng)
	!!call random_seed(put = [1, (0, i = 2, nrng)])  ! nvfortran complains with all 0's
	!call random_seed(put = [(0, i = 1, nrng)])

end subroutine rand_seed_determ

!********

double precision function rand_sca_f64()

	integer :: i

	! numa__rng

	!rand_sca_f64 = abs(1.d0 * global_rng%int32() / huge(i))
	!rand_sca_f64 = 1.d0 * mod(global_rng%int32(), huge(i)) / huge(i)
	rand_sca_f64 = 1.d0 * modulo(global_rng%int32(), huge(i)) / huge(i)

	!print *, "rand_sca_f64 = ", rand_sca_f64

	!! Built-in Fortran default rng
	!call random_number(rand_sca_f64)

end function rand_sca_f64

function rand_vec_f64(n)
	integer, intent(in) :: n
	double precision, allocatable :: rand_vec_f64(:)
	integer :: i

	allocate(rand_vec_f64(n))
	!call random_number(rand_vec_f64)

	! nvfortran rng sucks
	do i = 1, n
		rand_vec_f64(i) = rand_f64()
	end do

end function rand_vec_f64

function rand_mat_f64(m, n)
	integer, intent(in) :: m, n
	double precision, allocatable :: rand_mat_f64(:,:)
	integer :: i, j

	allocate(rand_mat_f64(m,n))
	!call random_number(rand_mat_f64)

	! nvfortran rng sucks
	do j = 1, n
	do i = 1, m
		rand_mat_f64(i,j) = rand_f64()
	end do
	end do

end function rand_mat_f64

!********

!> Random integer in the range 0 <= rand_i32() < n
integer function rand_i32(n)
	integer, intent(in) :: n

	rand_i32 = modulo(global_rng%int32(), n)
	!rand_i32 = mod(abs(global_rng%int32()), n)

	!rand_i32 = floor(n * rand_f64())

end function rand_i32

!===============================================================================

function rand_perm(n)
	! Fisher-Yates shuffle
	integer, intent(in) :: n
	integer, allocatable :: rand_perm(:)
	!********
	integer :: i, j

	rand_perm = [(i, i = 1, n)]  ! initially identity
	do i = 1, n-1
		j = i + rand_i32(n-i+1)
		rand_perm([i, j]) = rand_perm([j, i])
	end do

end function rand_perm

!===============================================================================

end module numa__rand

