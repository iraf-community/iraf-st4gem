include	<error.h>
include	"../fterr.h"

# FT_ERROR -- Fourier analysis error handler.
#
# This does two things : buffers an application-specific message and then
# issues the message and returns or aborts, depending on the disaster rating.
#
# C. D. Biemesderfer, STScI, 10 Dec 87

procedure ft_error (code, severity)

int	code		# Error code (defined in fterr.h)
int	severity	# Severity level (indicates action)
#--
char	emsg[SZ_LINE]

begin
	call strcpy ("FT pkg: ", emsg, SZ_LINE)

	switch (code) {
	case FT_NOINPUT:
	    call strcat ("input data unspecified or inaccessible",
			emsg, SZ_LINE)
	case FT_NOOUTPUT:
	    call strcat ("output unspecified", emsg, SZ_LINE)
	case FT_NOFTPAIRS:
	    call strcat ("ftpairs file not found", emsg, SZ_LINE)
	case FT_BADNAXIS:
	    call strcat ("illegal number of axes - must be 1 or 2",
			emsg, SZ_LINE)
	default:
	    call strcat ("unrecognized error", emsg, SZ_LINE)
	}

	switch (severity) {
	case FT_ERROR:
	    call error (0, emsg)
	case FT_FATAL:
	    call fatal (0, emsg)
	default:
	    call eprintf ("%s\n")
		call pargstr (emsg)
	    call erract (EA_WARN)
	}
end
