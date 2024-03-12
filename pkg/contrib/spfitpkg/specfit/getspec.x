###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	procedure getspec(specfile)
#		char	specfile[ARB]
#
#  Description: GETSPEC opens the IRAF or text-file IUE image specified by
#		specfile and stores the data in arrays resident in 
#		"specfit.com"
#		
#
#  Arguments:	char	specfile[ARB]	- Input file name
#
#  Returns:	none
#
#  Notes:	Shares data in "specfit.com"
#
#  History:	May	1989	Gerard Kriss
#		August	1989	Added Julian Date calculation for N5548 data
#		October 1989	Added option to open IRAF images as well as IUE
#		Oct  1  1991	gak	Updated  and simplified
#		April 	1995	grimes	Added default file extension and 
#					allowed to open hhh files
#
###########################################################################


include <ctype.h>
include "specfit.h"


procedure getspec(specfile)

char	specfile[ARB]

char 	file_extension[SZ_FNAME]
char 	tempspecfile[SZ_FNAME]
char	ferror[SZ_LINE]
int	imlen, length
int 	junk
int 	length_extension,length_file

#	imlen = the number of characters ".imh" should appear from
#	the end of an iraf file name. length = the number of characters
#	that appear in the file name. 

int 	envgets()
bool 	access()
int 	strlen()

include "specfit.com"

begin 
	length = strlen(specfile)
	imlen = length - 3

	junk = envgets("imtype",file_extension,SZ_FNAME)
	
	length_extension = strlen(file_extension)
	length_file = strlen(specfile)
	call strcpy(specfile,tempspecfile,SZ_FNAME)
	if (file_extension[1]=='.') {
		length_extension=length_extension-1
	} else {
		call strcat(".",tempspecfile,SZ_FNAME)
	}


	call strcat(file_extension,tempspecfile,SZ_FNAME)

	if (access(specfile,0,0)) {
				#if exactly as user inputted file name
		if ((specfile[length_file]=='h' || specfile[length_file]=='H')
			&& specfile[length_file-3]=='.') {
			call rdiraf(specfile) #file is iraf
		} else {
			call rdspec(specfile) #file is specfit format
		}
	} else if (access(tempspecfile,0,0)) {
				#input filename+default file extension
				#is a file that exists-will open
		if ((tempspecfile[length_extension+length_file+1]=='h' 
			|| tempspecfile[length_extension+length_file+1]=='H')
			 && length_extension==3) {
			call rdiraf(tempspecfile) #this file is iraf
		} else {
			call rdspec(tempspecfile) #this file is specfit
		}
	} else {
		call sprintf(ferror,SZ_LINE,"Can't find either a file named %s or %s!!!!")
		call pargstr(tempspecfile)
		call pargstr(specfile)
		call error(1, ferror)
	}

end

procedure rdiraf(specfile)
char	specfile[ARB]

include "specfit.com"

begin
	call getimage(specfile)
end 

procedure rdspec(specfile)

char	 specfile[ARB]

char	 skips[80], hm[20], arg[10]
int	 i, istat, fd
int	 ic, ip, nc
int	 year, day, hr, min
double	jd3
int	 open(), fscan(), stridx(), ctoi()

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
		jd3 = 0.0d+00
		istat = fscan(fd)
#		Get the number of spectral points and the integration time
		    call gargi(npts)
		    call gargr(itime)
		    call gargd(jd3)

# Generate the Julian Date and enter into the log
		ic = stridx(":",hm)
		call strcpy(hm,arg,ic-1)
		ip = 1
		nc = ctoi(arg, ip, hr)
		nc = ctoi(hm, ic+1, min)

		jd = 0.0d+00

		call julian( year, day, hr, min, jd)

		for ( i = 1; i <= nlogfd; i = i + 1) {
			call fprintf(logfd[i],"\nJD = %10.3f\n")
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
        call eprintf("%d points expected, but %d points read from %s./n")
			call pargi(npts)
    			call pargi(i)
    			call pargstr(specfile)
		}

		call close( fd )

	end 

