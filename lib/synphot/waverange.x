include <tbset.h>
include "libsynphot.h"

# WAVERANGE -- Compute wavelength range of instrument mode

procedure waverange (mxfile, ncomp, filelist, minwave, maxwave)

int	mxfile			# i: maximum length of file name
int	ncomp			# i: number of instrument components
char	filelist[mxfile,ARB]	# i: list of component file names
real	minwave			# o: minimum wavelength
real	maxwave			# o: maximum wavelength
#--
int	icomp, status
real	onelow, onehigh

string	nocomp   "No components in optical path"
string	badfile  "Cannot compute wavelength range for file"
	
errchk	tabrange

begin
	# Loop over all instrument components

	if (ncomp <= 0)
	    call synphoterr (nocomp, "waverange")

	do icomp = 1, ncomp {
	    # Compute wavelength limits for table

	    call tabrange (filelist[1,icomp], onelow, onehigh, status)
	    if (status == ERR)
		call synphoterr (badfile, filelist[1,icomp])

	    # Take intersection of range with previous ranges

	    if (icomp == 1) {
		minwave = onelow
		maxwave = onehigh
	    } else {
		minwave = max (minwave, onelow)
		maxwave = min (maxwave, onehigh)
	    }
	}

end


