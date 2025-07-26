
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

end module numa__utils

