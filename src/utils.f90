
!> Common utilities with minimal dependencies
module numa__utils

	implicit none

	character, parameter :: &
			NULL_CHAR       = char( 0), &
			TAB             = char( 9), &
			LINE_FEED       = char(10), &
			VERT_TAB        = char(11), &
			CARRIAGE_RETURN = char(13), &
			ESC             = char(27)

	character(len = *), parameter :: &
			BOLD               = ESC//"[;1m", &
			YELLOW             = ESC//"[33m", &
			!BRIGHT_RED         = ESC//"[91m", &
			RED                = ESC//"[91;1m", &
			BOLD_BRIGHT_YELLOW = ESC//"[93;1m", &
			GREEN              = ESC//"[92;1m", &
			BRIGHT_YELLOW      = ESC//"[93m", &
			BLUE               = ESC//"[94m", &
			MAGENTA            = ESC//"[95m", &
			CYAN               = ESC//"[96m", &
			WHITE              = ESC//"[97m", &
			COLOR_RESET        = ESC//"[0m"

	character(len = *), parameter :: &
		ERROR = RED // "Error: " // COLOR_RESET

	interface to_str
		procedure :: to_str_i32
		procedure :: to_str_i64
		procedure :: to_str_f32
		procedure :: to_str_f64
	end interface to_str

	interface sorted
		! Functions
		procedure :: sorted_f64
		procedure :: sorted_c64
	end interface

	interface sort
		! Subroutines
		procedure :: sort_f64
		procedure :: sort_c64
	end interface

	interface compare_lex
		procedure :: compare_lex  ! TODO: rename routines that clash with the overload
		procedure :: compare_lex_c64
	end interface
	interface lt_lex
		procedure :: lt_lex
		procedure :: lt_lex_c64
	end interface
	interface lte_lex
		procedure :: lte_lex
		procedure :: lte_lex_c64
	end interface

	integer, parameter :: &
		COMPARE_LESS_ = -1, &
		COMPARE_EQUAL_ = 0, &
		COMPARE_GREATER_ = 1

contains

!===============================================================================

function select_max_n(v, n) result(iv)
	! Select the `n` largest values from `v` and return their indices `iv`
	!
	! Use a shortened bubble-sort based method O(n * nv).  Better selection
	! algorithms exist

	double precision, intent(in) :: v(:)
	integer, intent(in) :: n

	integer, allocatable :: iv(:)

	!********

	integer :: i, j, nv

	nv = size(v)
	iv = [(i, i = 1, nv)]

	do i = 1, n
	do j = nv-1, i, -1
		if (v(iv(j)) < v(iv(j+1))) then
			iv([j, j+1]) = iv([j+1, j])
		end if
	end do
	end do

	! Trim
	iv = iv(1: n)

end function select_max_n

!===============================================================================

function to_str_i32(val) result(str)
	integer(kind = 4), intent(in) :: val
	character(len = :), allocatable :: str
	character :: buffer*16
	write(buffer, "(i0)") val
	str = trim(buffer)
end function to_str_i32

!********

function to_str_i64(val) result(str)
	integer(kind = 8), intent(in) :: val
	character(len = :), allocatable :: str
	character :: buffer*16
	write(buffer, "(i0)") val
	str = trim(buffer)
end function to_str_i64

!********

function to_str_f32(val) result(str)
	real, intent(in) :: val
	character(len = :), allocatable :: str
	character :: buffer*32
	write(buffer, *) val
	str = trim(adjustl(buffer))
end function to_str_f32

!********

function to_str_f64(val) result(str)
	double precision, intent(in) :: val
	character(len = :), allocatable :: str
	character :: buffer*32
	write(buffer, *) val
	str = trim(adjustl(buffer))
end function to_str_f64

!===============================================================================

integer function compare_lex(a, b) result(compare)

	! Lexicographically compare int vectors a and b.  Return COMPARE_LESS_ if a < b,
	! COMPARE_EQUAL_ if a == b, or COMPARE_GREATER_ if a > b.
	!
	! TODO: rename `compare_lex` and related routines like `lt_lex` with type in
	! name for possible future overloading

	integer, intent(in) :: a(:), b(:)
	!integer, intent(in) :: a, b

	integer :: i

	if (size(a) /= size(b)) then
		write(*,*) 'Error: incompatible sizes of args in compare_lex()'
		! TODO: panic
		stop
	end if

	do i = 1, size(a)
		if (a(i) < b(i)) then
			compare = COMPARE_LESS_
			return
		else if (a(i) > b(i)) then
			compare = COMPARE_GREATER_
			return
		end if
	end do

	compare = COMPARE_EQUAL_

end function compare_lex

integer function compare_lex_c64(a, b) result(compare)

	! Lexicographically compare int vectors a and b.  Return COMPARE_LESS_ if a < b,
	! COMPARE_EQUAL_ if a == b, or COMPARE_GREATER_ if a > b.

	double complex, intent(in) :: a, b

	if (a%re < b%re) then
		compare = COMPARE_LESS_
		return
	else if (a%re > b%re) then
		compare = COMPARE_GREATER_
		return
	end if
	if (a%im < b%im) then
		compare = COMPARE_LESS_
		return
	else if (a%im > b%im) then
		compare = COMPARE_GREATER_
		return
	end if

	compare = COMPARE_EQUAL_

end function compare_lex_c64

!===============================================================================

logical function lt_lex(a, b) result(lt)
	integer, intent(in) :: a(:), b(:)
	lt = compare_lex(a, b) == COMPARE_LESS_
end function lt_lex

logical function lt_lex_c64(a, b) result(lt)
	double complex, intent(in) :: a, b
	lt = compare_lex(a, b) == COMPARE_LESS_
end function lt_lex_c64

!===============================================================================

logical function lte_lex(a, b) result(lte)
	integer, intent(in) :: a(:), b(:)
	integer :: comp
	comp = compare_lex(a, b)
	lte = comp == COMPARE_LESS_ .or. comp == COMPARE_EQUAL_
end function lte_lex

logical function lte_lex_c64(a, b) result(lte)
	double complex, intent(in) :: a, b
	integer :: comp
	comp = compare_lex(a, b)
	lte = comp == COMPARE_LESS_ .or. comp == COMPARE_EQUAL_
end function lte_lex_c64

!===============================================================================

logical function gt_lex(a, b) result(gt)

	integer, intent(in) :: a(:), b(:)

	gt = compare_lex(a, b) == COMPARE_GREATER_

end function gt_lex

!===============================================================================

logical function eq_lex(a, b) result(eq)

	integer, intent(in) :: a(:), b(:)

	eq = compare_lex(a, b) == COMPARE_EQUAL_

end function eq_lex

!===============================================================================

recursive subroutine sortidx_i32_2(a, idx, lo_in, hi_in)

	! Quicksort a rank-2 array a along its 2nd dimension and return the
	! sort permutation idx
	!
	! See also:  https://github.com/JeffIrwin/aoc-2022/blob/381da3c2d468b4e2a9d4bb1068d84a9e8ae6bec6/2022/23/main.f90#L114

	integer, intent(in) :: a(:,:)
	integer, allocatable :: idx(:)
	integer, optional :: lo_in, hi_in

	integer :: lo, hi, p, i

	logical :: outer

	lo = 1
	hi = size(a, 2)
	outer = .true.

	if (present(lo_in)) lo = lo_in
	if (present(hi_in)) then
		hi = hi_in
		outer = .false.
	end if

	if (lo >= hi .or. lo < 1) return

	if (.not. allocated(idx)) then
		idx = [(i, i = 1, size(a, 2))]
	end if

	p = partition_i32_2(a, idx, lo, hi)

	call sortidx_i32_2(a, idx, lo, p - 1)
	call sortidx_i32_2(a, idx, p + 1, hi)

	!if (outer) then
	!	! Check result.  TODO debug only
	!	do i = 2, size(a, 2)
	!		if (lt_lex(a(:,idx(i)), a(:,idx(i-1)))) then
	!			write(*,*) 'Error in sortidx_i32_2'
	!			stop
	!		end if
	!	end do
	!end if

end subroutine sortidx_i32_2

!===============================================================================

integer function partition_i32_2(a, idx, lo, hi) result(ans)

	integer, intent(in) :: a(:,:)
	integer, allocatable :: idx(:), pivot(:)
	integer :: lo, hi, i, j, mid

	!pivot = a(:, idx(hi))

	mid = (lo + hi) / 2
	if (lt_lex(a(:, idx(mid)), a(:, idx(lo)))) then
		idx([lo, mid]) = idx([mid, lo])
	else if (lt_lex(a(:, idx(hi)), a(:, idx(lo)))) then
		idx([lo, hi]) = idx([hi, lo])
	else if (lt_lex(a(:, idx(mid)), a(:, idx(hi)))) then
		idx([mid, hi]) = idx([hi, mid])
	end if
	pivot = a(:, idx(hi))

	i = lo - 1
	do j = lo, hi - 1
	!do j = lo, hi
		if (lte_lex(a(:, idx(j)), pivot)) then
			i = i + 1
			idx([i, j]) = idx([j, i])
		end if
	end do

	i = i + 1
	idx([i, hi]) = idx([hi, i])
	ans = i

end function partition_i32_2

!===============================================================================

subroutine sort_f64(a)
	! Replace `a` with its ascending sort
	double precision, intent(inout) :: a(:)

	integer, allocatable :: idx(:)
	call sortidx_f64_1(a, idx)
	a = a(idx)

end subroutine sort_f64

subroutine sort_c64(a)
	! Replace `a` with its ascending sort
	double complex, intent(inout) :: a(:)

	integer, allocatable :: idx(:)
	call sortidx_c64_1(a, idx)
	a = a(idx)

end subroutine sort_c64

function sorted_f64(a) result(s)
	! Return `a` sorted in ascending order
	double precision, intent(in) :: a(:)
	double precision, allocatable :: s(:)

	s = a
	call sort_f64(s)

end function sorted_f64

function sorted_c64(a) result(s)
	! Return `a` sorted in ascending order with lexicographical ordering for
	! complex values
	double complex, intent(in) :: a(:)
	double complex, allocatable :: s(:)

	s = a
	call sort_c64(s)

end function sorted_c64

recursive subroutine sortidx_f64_1(a, idx, lo_in, hi_in)

	! Quicksort a rank-2 array a along its 2nd dimension and return the
	! sort permutation idx
	!
	! See also:  https://github.com/JeffIrwin/aoc-2022/blob/381da3c2d468b4e2a9d4bb1068d84a9e8ae6bec6/2022/23/main.f90#L114

	double precision, intent(in) :: a(:)
	integer, allocatable :: idx(:)
	integer, optional :: lo_in, hi_in

	integer :: lo, hi, p, i

	logical :: outer

	lo = 1
	hi = size(a, 1)
	outer = .true.

	if (present(lo_in)) lo = lo_in
	if (present(hi_in)) then
		hi = hi_in
		outer = .false.
	end if

	if (lo >= hi .or. lo < 1) return

	if (.not. allocated(idx)) then
		idx = [(i, i = 1, size(a, 1))]
	end if

	p = partition_f64_1(a, idx, lo, hi)

	call sortidx_f64_1(a, idx, lo, p - 1)
	call sortidx_f64_1(a, idx, p + 1, hi)

end subroutine sortidx_f64_1

recursive subroutine sortidx_c64_1(a, idx, lo_in, hi_in)

	! Quicksort a rank-2 array a along its 2nd dimension and return the
	! sort permutation idx
	!
	! See also:  https://github.com/JeffIrwin/aoc-2022/blob/381da3c2d468b4e2a9d4bb1068d84a9e8ae6bec6/2022/23/main.f90#L114

	double complex, intent(in) :: a(:)
	integer, allocatable :: idx(:)
	integer, optional :: lo_in, hi_in

	integer :: lo, hi, p, i

	logical :: outer

	lo = 1
	hi = size(a, 1)
	outer = .true.

	if (present(lo_in)) lo = lo_in
	if (present(hi_in)) then
		hi = hi_in
		outer = .false.
	end if

	if (lo >= hi .or. lo < 1) return

	if (.not. allocated(idx)) then
		idx = [(i, i = 1, size(a, 1))]
	end if

	p = partition_c64_1(a, idx, lo, hi)

	call sortidx_c64_1(a, idx, lo, p - 1)
	call sortidx_c64_1(a, idx, p + 1, hi)

end subroutine sortidx_c64_1

!===============================================================================

integer function partition_f64_1(a, idx, lo, hi) result(ans)

	double precision, intent(in) :: a(:)
	integer, allocatable :: idx(:)  ! TODO: intent, allocatable
	double precision :: pivot
	integer :: lo, hi, i, j, mid

	!pivot = a(:, idx(hi))

	mid = (lo + hi) / 2
	if (a(idx(mid)) < a(idx(lo))) then
		idx([lo, mid]) = idx([mid, lo])
	else if (a(idx(hi)) < a(idx(lo))) then
		idx([lo, hi]) = idx([hi, lo])
	else if (a(idx(mid)) < a(idx(hi))) then
		idx([mid, hi]) = idx([hi, mid])
	end if
	pivot = a(idx(hi))

	i = lo - 1
	do j = lo, hi - 1
	!do j = lo, hi
		if (a(idx(j)) <= pivot) then
			i = i + 1
			idx([i, j]) = idx([j, i])
		end if
	end do

	i = i + 1
	idx([i, hi]) = idx([hi, i])
	ans = i

end function partition_f64_1

integer function partition_c64_1(a, idx, lo, hi) result(ans)

	double complex, intent(in) :: a(:)
	integer, allocatable :: idx(:)  ! TODO: intent, allocatable
	double complex :: pivot
	integer :: lo, hi, i, j, mid

	!pivot = a(:, idx(hi))

	mid = (lo + hi) / 2
	if (lt_lex(a(idx(mid)), a(idx(lo)))) then
		idx([lo, mid]) = idx([mid, lo])
	else if (lt_lex(a(idx(hi)), a(idx(lo)))) then
		idx([lo, hi]) = idx([hi, lo])
	else if (lt_lex(a(idx(mid)), a(idx(hi)))) then
		idx([mid, hi]) = idx([hi, mid])
	end if
	pivot = a(idx(hi))

	i = lo - 1
	do j = lo, hi - 1
	!do j = lo, hi
		if (lte_lex(a(idx(j)), pivot)) then
			i = i + 1
			idx([i, j]) = idx([j, i])
		end if
	end do

	i = i + 1
	idx([i, hi]) = idx([hi, i])
	ans = i

end function partition_c64_1

!===============================================================================

end module numa__utils

