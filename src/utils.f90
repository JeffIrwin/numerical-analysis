
module utils_m

	implicit none

	character, parameter :: &
			NULL_CHAR       = char( 0), &
			TAB             = char( 9), &
			LINE_FEED       = char(10), &
			VERT_TAB        = char(11), &
			CARRIAGE_RETURN = char(13), &
			ESC             = char(27)

	! TODO: make these variables, with colors disabled if output_unit is not tty
	! and an option to --force-color
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

end module utils_m

