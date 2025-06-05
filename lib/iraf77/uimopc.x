include	<iraf77.h>
include <imio.h>
include <syserr.h>

# UIMOPC -- Open a new file and copy an existing image into it.
# f77nam is the name of the new file and is of the form:
# filename[p/n] where n is the number of images in the group-file.
# or filename[*] which means the output file must contain the same 
# number of images than the input file.
#
# NOTE: will be modified when template expander fully implemented.

procedure uimopc (f77nam, temp_imid, new_imid, istat)

%	character*(*) f77nam
pointer	temp_imid		# descriptor for template image 
pointer new_imid		# descriptor for new image
int	istat

char 	imname[SZ_PATHNAME]
pointer immap()
int	ierr
int	errcode()

begin
	istat = ER_OK
	# Convert character string to SPP string
	call f77upk (f77nam, imname, SZ_PATHNAME)

	iferr ( new_imid = immap (imname, NEW_COPY, temp_imid)) {
	   ierr = errcode()
	   if (ierr == SYS_IMMAGNCPY)
		istat = ER_IMNOTIMAG
	   if (ierr == SYS_IMSECTNEWIM)
		istat = ER_IMILLSEC
	   if (ierr == SYS_IKIEXTN)
		istat = ER_IMBADEXTN
	   if (ierr == SYS_IKIOPEN)
	     	istat = ER_IMOPCOP
	   # set pointer to zero to avoid passing back an invalid pointer
	   new_imid = 0
	}
	return

end
