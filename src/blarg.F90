
#include "panic.F90"

!> Basic Linear Algebra Routine Group
!>
!> > blarg:  (informal) Expressing frustration or disappointment.
!> >
!> > Blarg! I'm sick of this.
module numa__blarg

	implicit none

	interface outer_product
		procedure :: outer_product_c64
		procedure :: outer_product_f64
	end interface outer_product

	interface zeros
		procedure :: zeros_vec
		procedure :: zeros_mat
	end interface zeros

	interface zeros_i32
		procedure :: zeros_vec_i32
		!procedure :: zeros_mat_i32
	end interface zeros_i32

	interface diag
		procedure :: diag_set
		procedure :: diag_get
		procedure :: diag_get_c64
	end interface diag

	interface triu
		procedure :: triu_c64
		procedure :: triu_f64
	end interface triu

	interface hstack
		procedure :: hstack_mat_mat
		procedure :: hstack_mat_vec
	end interface hstack

	interface vstack
		procedure :: vstack_mat_mat
		procedure :: vstack_mat_vec
	end interface vstack

	interface mask_to_index
		procedure :: mask_to_index_vec
		procedure :: mask_to_index_mat
	end interface mask_to_index

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

function zeros_mat(m, n) result(a)
	! m x n matrix of 0
	integer, intent(in) :: m, n
	double precision, allocatable :: a(:,:)
	allocate(a(m, n))
	a = 0
end function zeros_mat

function zeros_vec(n) result(a)
	! Size n vector of 0
	integer, intent(in) :: n
	double precision, allocatable :: a(:)
	allocate(a(n))
	a = 0
end function zeros_vec

function zeros_vec_i32(n) result(a)
	! Size n vector of 0
	integer, intent(in) :: n
	integer, allocatable :: a(:)
	allocate(a(n))
	a = 0
end function zeros_vec_i32

!********

function diag_set(v) result(d)
	! Spread a diagonal vector `v` into a matrix
	double precision, intent(in) :: v(:)
	double precision, allocatable :: d(:,:)
	!********
	integer :: i, n

	n = size(v)
	d = zeros(n, n)
	do i = 1, n
		d(i,i) = v(i)
	end do

end function diag_set

function diag_get(a) result(v)
	! Get the diagonal vector `v` from a matrix `a`
	double precision, intent(in) :: a(:,:)
	double precision, allocatable :: v(:)
	!********
	integer :: i, n

	n = min(size(a,1), size(a,2))
	allocate(v(n))
	do i = 1, n
		v(i) = a(i,i)
	end do

end function diag_get

function diag_get_c64(a) result(v)
	! Get the diagonal vector `v` from a matrix `a`
	double complex, intent(in) :: a(:,:)
	double complex, allocatable :: v(:)
	!********
	integer :: i, n

	n = min(size(a,1), size(a,2))
	allocate(v(n))
	do i = 1, n
		v(i) = a(i,i)
	end do

end function diag_get_c64

!===============================================================================

function triu_f64(a) result(r)
	! Explicitly get R by zeroing lower triangle
	double precision, intent(in) :: a(:,:)
	double precision, allocatable :: r(:,:)
	!********
	integer :: i, n

	n = min(size(a,1), size(a,2))
	r = a(:n, :n)
	do i = 1, n
		r(i+1:, i) = 0
	end do

end function triu_f64

!********

function triu_c64(a) result(r)
	! Explicitly get R (U) by zeroing lower triangle
	double complex, intent(in) :: a(:,:)
	double complex, allocatable :: r(:,:)
	!********
	integer :: i

	r = a
	do i = 1, size(a, 1)
		r(i+1:, i) = 0
	end do

end function triu_c64

!===============================================================================

function vstack_mat_vec(a, b, iostat) result(c)
	use numa__utils
	double precision, intent(in) :: a(:,:), b(:)
	double precision, allocatable :: c(:,:)
	integer, optional, intent(out) :: iostat
	!********
	character(len = :), allocatable :: msg
	integer :: na
	integer, parameter :: nb = 1

	if (present(iostat)) iostat = 0

	! Sizes must match, unless at least one of the arrays is empty
	if (size(a,2) /= size(b) .and. size(a) > 0 .and. size(b) > 0) then
		msg = "size(a,2) does not match size(b) in vstack_mat_vec()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	na = size(a, 1)
	allocate(c(na+nb, size(a,2)))
	if (size(a,2) <= 0) return

	if (na > 0) c(1:na, :) = a
	if (nb > 0) c(na+1, :) = b

end function vstack_mat_vec

!********

function vstack_mat_mat(a, b, iostat) result(c)
	use numa__utils
	double precision, intent(in) :: a(:,:), b(:,:)
	double precision, allocatable :: c(:,:)
	integer, optional, intent(out) :: iostat
	!********
	character(len = :), allocatable :: msg
	integer :: na, nb

	if (present(iostat)) iostat = 0
	if (size(a,2) /= size(b,2) .and. size(a) > 0 .and. size(b) > 0) then
		msg = "size(a,2) does not match size(b,2) in vstack_mat_mat()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	na = size(a, 1)
	nb = size(b, 1)

	allocate(c(na+nb, size(a,2)))
	if (size(a,2) <= 0) return

	if (na > 0) c(1:na , :) = a
	if (nb > 0) c(na+1:, :) = b

end function vstack_mat_mat

!===============================================================================

function hstack_mat_mat(a, b, iostat) result(c)
	! `hstack` and `vstack` are named after corresponding numpy functions
	use numa__utils
	double precision, intent(in) :: a(:,:), b(:,:)
	double precision, allocatable :: c(:,:)
	integer, optional, intent(out) :: iostat
	!********
	character(len = :), allocatable :: msg
	integer :: na, nb

	if (present(iostat)) iostat = 0
	if (size(a,1) /= size(b,1) .and. size(a) > 0 .and. size(b) > 0) then
		msg = "size(a,1) does not match size(b,1) in hstack_mat_mat()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	na = size(a, 2)
	nb = size(b, 2)
	allocate(c(size(a,1), na+nb))
	if (size(a,1) <= 0) return

	if (na > 0) c(:, 1: na) = a
	if (nb > 0) c(:, na+1:) = b

end function hstack_mat_mat

!===============================================================================

function hstack_mat_vec(a, b, iostat) result(c)
	use numa__utils
	double precision, intent(in) :: a(:,:), b(:)
	double precision, allocatable :: c(:,:)
	integer, optional, intent(out) :: iostat
	!********
	character(len = :), allocatable :: msg
	integer :: na, nb

	if (present(iostat)) iostat = 0
	if (size(a,1) /= size(b) .and. size(a) > 0 .and. size(b) > 0) then
		msg = "size(a,1) does not match size(b) in hstack_mat_vec()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	na = size(a, 2)
	nb = 1

	allocate(c(size(a,1), na+nb))
	if (size(a,1) <= 0) return

	if (na > 0) c(:, 1: na) = a
	if (nb > 0) c(:, na+1 ) = b

end function hstack_mat_vec

!===============================================================================

function mask_to_index_vec(mask) result(indices)
	! Convert a logical mask array to an index array `indices`
	!
	! For example:
	!
	!     mask_to_index([.true., .false., .false, .true.])
	!     ==
	!     [1, 4]
	!
	! This is useful if you actually need the index array.  If you just want the
	! result of indexing the parent array by the returned index array, you can
	! just use built-in pack(parent, mask) instead

	logical, intent(in) :: mask(:)
	integer, allocatable :: indices(:)
	!********
	integer :: i, j

	! Overallocate and then trim.  Could make two passes, explicitly or
	! implicitly using count()
	allocate(indices(size(mask)))
	i = 0
	do j = 1, size(mask)
		if (.not. (mask(j))) cycle
		i = i + 1
		indices(i) = j
	end do
	indices = indices(:i)  ! trim

end function mask_to_index_vec

!===============================================================================

function mask_to_index_mat(mask) result(indices)
	! Convert a logical mask array to an index array `indices`
	!
	! For example:
	!
	!     mask_to_index([.true., .false., .false, .true.])
	!     ==
	!     [1, 4]
	!
	! This is useful if you actually need the index array.  If you just want the
	! result of indexing the parent array by the returned index array, you can
	! just use built-in pack(parent, mask) instead

	logical, intent(in) :: mask(:,:)
	integer, allocatable :: indices(:,:)
	!********
	integer :: i, j, k

	! Overallocate and then trim.  Could make two passes, explicitly or
	! implicitly using count()
	allocate( indices(2, size(mask)) )
	i = 0
	do j = 1, size(mask,2)
	do k = 1, size(mask,1)
		if (.not. (mask(k,j))) cycle
		i = i + 1
		indices(:,i) = [k, j]
	end do
	end do
	indices = indices(:, :i)  ! trim

end function mask_to_index_mat

!===============================================================================

function reverse(a) result(r)
	integer, intent(in) :: a(:)
	integer, allocatable :: r(:)
	!********
	integer :: i

	r = a([(i, i = size(a), 1, -1)])

end function reverse

!===============================================================================

end module numa__blarg

