###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University


#  Synopsis:	procedure getspec(specfile)
#		char	specfile[ARB]

#  Description:	GETSPEC opens the text-file IUE image specified by specfile and
#		stores the data in arrays resident in "specfit.com"
#		Intended exclusively for analysis of NGC5548 campaign data.
#		

#  Arguments:	char	specfile[ARB]	- Input file name

#  Returns:	none

#  Notes:	Shares data in "specfit.com"

#  History:	May	1989	Gerard Kriss
#		August	1989	Added Julian Date calculation for N5548 data
#		10/20/89	gak	Put JD in specfit.com

###########################################################################

include	<ctype.h>
include	"specfit.h"

procedure getspec(specfile)

char	specfile[ARB]

char	skips[80], hm[20], arg[10]
int	i, istat, fd
int	ic, ip, nc
int	year, day, hr, min
int	open(), fscan(), stridx(), ctoi()

include	"specfit.com"

begin

# Open the file
	fd = open(specfile, READ_ONLY, TEXT_FILE)

# Read header lines and get number of points
	istat = fscan(fd)	# Get the object name
	    call gargstr(specname, 20)
	    call gargstr(skips,16)
	    call gargi(year)
	    call gargstr(skips,6)
	    call gargi(day)
	    call gargstr(hm,8)
	istat = fscan(fd)	# Get the number of spectral points
	    call gargi(npts)

# Generate the Julian Date and enter into the log
	ic = stridx(":",hm)
	call strcpy(hm,arg,ic-1)
	ip = 1
	nc = ctoi(arg, ip, hr)
	nc = ctoi(hm, ic+1, min)

	jd = 0.0d+00


	call julian(year,day,hr,min,jd)	

	for ( i = 1; i <= nlogfd; i = i + 1) {
		call fprintf(logfd[i],"\nJD=%10.3f\n")
		    call pargd(jd)
	}

# Allocate space for data
	call malloc(spectrum, npts, TY_REAL)
	call malloc(lambdas, npts, TY_REAL)
	call malloc(errors, npts, TY_REAL)

# Read in the data
	i = 0
	while ( fscan(fd) != EOF ) {
		call gargr(Memr[lambdas+i])
		call gargr(Memr[spectrum+i])
		call gargr(Memr[errors+i])
		i = i + 1
	}
	if ( i != npts ) {
	    call eprintf("%d points expected, but %d points read from %s.\n")
		call pargi(npts)
	    	call pargi(i)
	    	call pargstr(specfile)
	}

	call close( fd )

end




