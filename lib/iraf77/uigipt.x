include <imio.h>
include <imset.h>
include <iraf77.h>

# UIGIPT--Get image processing parameter of type character

procedure uigipt (im, parnum, buffer, istat)

pointer im	# image descriptor pointer returned by uimopn
int parnum	# symbolic code; image parameter to get
%	character*(*) buffer
int istat	# return status code

char outstr[SZ_IMNAME] # SPP char string
int maxch

# For the present, only the image name in imio can be
# set or gotten by the user in imio.

begin

	switch (parnum){
	case IM_IMAGENAME:
		maxch = SZ_IMNAME
		iferr (call imstats (im, parnum, outstr, maxch)){

			# Return a Fortran blank-filled character string 
			# if an error occurred
			istat = ER_IMSTATUNKPAR
%			buffer = ' '
			return
		}

		call f77pak (outstr, buffer, maxch)
		istat = ER_OK

	default:
		# Return a Fortran blank-filled character string 
		# if the user asked for an unknown parameter
%		buffer = ' '
		istat = ER_IMSTATUNKPAR
	}

	return

end
