
#include "panic.F90"

!> Basic linear algebra (sub)routines with minimal dependencies
!>
!> This actually has nothing to do with "BLAS" but it's an ok name.  I could
!> call it math, but this is all math.  Even linalg is too generic.
module numa__blas

	implicit none

	interface outer_product
		procedure :: outer_product_c64
		procedure :: outer_product_f64
	end interface outer_product

contains

!===============================================================================

function outer_product_f64(a, b) result(c)
	double precision, intent(in) :: a(:), b(:)
	double precision, allocatable :: c(:,:)

	integer :: i, j, na, nb
	na = size(a)
	nb = size(b)

	allocate(c(na, nb))
	do j = 1, nb
	do i = 1, na
		c(i,j) = a(i) * b(j)
	end do
	end do

end function outer_product_f64

function outer_product_c64(a, b) result(c)
	double complex, intent(in) :: a(:), b(:)
	double complex, allocatable :: c(:,:)

	integer :: i, j, na, nb
	na = size(a)
	nb = size(b)

	allocate(c(na, nb))
	do j = 1, nb
	do i = 1, na
		!c(i,j) = a(i) * conjg(b(j))
		c(i,j) = a(i) * b(j)
	end do
	end do

end function outer_product_c64

!********

! TODO: it makes sense for zeros() to be defined in here too
function eye(n)
	! n x n identity matrix

	integer, intent(in) :: n
	double precision, allocatable :: eye(:,:)

	integer :: i, j
	allocate(eye(n, n))
	do i = 1, n
		do j = 1, n
			if (i == j) then
				eye(i,j) = 1.d0
			else
				eye(i,j) = 0.d0
			end if
		end do
	end do

end function eye

!********

double precision function rand_f64()
	call random_number(rand_f64)
end function rand_f64

!********

subroutine house_c64(alpha, x, diag_, iostat)
	! Return the Householder reflector represting `pp` such that pp * x == [1, 0, 0, ...]
	!
	! This was adapted from reference LAPACK:  https://netlib.org/lapack/explore-html/d8/d0d/group__larfg_ga11ff37151d75113678103648c1a68c11.html
	!
	! The diagonal "tau" factor would've been impossible to figure out
	! otherwise.  I had an unpacked version working which exploded Q and R into
	! separate arrays, but keeping them both in the A matrix (and diagonal
	! vector) was too hard without looking at LAPACK

	use numa__utils
	double complex, intent(inout) :: alpha
	double complex, intent(inout) :: x(:)
	double complex, intent(out) :: diag_
	integer, optional, intent(out) :: iostat
	!********

	character(len = :), allocatable :: msg
	double precision :: xnorm, beta, alphr, alphi
	integer :: n

	if (present(iostat)) iostat = 0

	n = size(x)

	xnorm = norm2c(x)
	alphr = alpha%re
	alphi = alpha%im

	if (xnorm == 0 .and. alphi == 0) then
		! For eig_lapack(), it's very important that house_f64() does not panic
		! here.  Maybe house_c64() should do the same
		diag_ = 0
		msg = "vector is singular in house_c64()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if
	beta = -sign_(alphr) * norm2([alphr, alphi, xnorm])

	! Skip checking if abs(beta) < machine_precision
	diag_ = dcmplx((beta - alphr) / beta, -alphi / beta)
	alpha = 1.d0 / (alpha - beta)
	x = x * alpha

	alpha = beta

end subroutine house_c64

!********

subroutine house_f64(alpha, x, diag_, iostat)
	! Return the Householder reflector represting `pp` such that pp * x == [1, 0, 0, ...]
	!
	! This was adapted from reference LAPACK:  https://netlib.org/lapack/explore-html/d8/d0d/group__larfg_ga11ff37151d75113678103648c1a68c11.html
	!
	! The diagonal "tau" factor would've been impossible to figure out
	! otherwise.  I had an unpacked version working which exploded Q and R into
	! separate arrays, but keeping them both in the A matrix (and diagonal
	! vector) was too hard without looking at LAPACK

	use numa__utils
	double precision, intent(inout) :: alpha
	double precision, intent(inout) :: x(:)
	double precision, intent(out) :: diag_
	integer, optional, intent(out) :: iostat
	!********

	!character(len = :), allocatable :: msg
	double precision :: xnorm, beta
	integer :: n

	if (present(iostat)) iostat = 0

	n = size(x)

	xnorm = norm2(x)

	if (xnorm <= 0) then
		! It's important for eig_lapack() to not panic here
		diag_ = 0
		return

		!diag_ = 0
		!msg = "vector is singular in house_f64()"
		!print *, "x = ", x
		!call PANIC(msg, present(iostat))
		!iostat = 1
		!return
	end if
	beta = -sign_(alpha) * norm2([alpha, xnorm])

	! Skip checking if abs(beta) < machine_precision
	diag_ = (beta - alpha) / beta

	alpha = 1.d0 / (alpha - beta)
	x = x * alpha

	alpha = beta

end subroutine house_f64

!===============================================================================

double precision function norm2c(v)

	double complex, intent(in) :: v(:)

	! It's weird that Fortran's intrinsic norm2() can't handle complex args.
	! Intel MKL has dznrm2() which is equivalent
	!
	! Note that dot_product() already conjugates one of the arguments, so there
	! is no need for additional conjugation

	!print *, "v = ", v
	norm2c = dble(sqrt(dot_product(v, v)))

end function norm2c

!===============================================================================

end module numa__blas

