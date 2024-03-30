include "libsynphot.h"

# BANDFUNC -- Compute the bandpass for a specified instrument mode

procedure bandfunc (mode, graphtab, comptab, nwave, wave, band)

char	mode[ARB]	# i: instrument mode
char	graphtab[ARB]	# i: graph table name
char	comptab[ARB]	# i: component lookup table name
int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelengths at which bandpass is computed
real	band[ARB]	# o: instrument bandpass
#--
bool	verbose
int	ncomp
pointer	sp, nparam, paramlist, filelist, banderr

data	verbose	/ false /

errchk	searchgraf, getthruput

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (nparam, MAXLIST, TY_INT)
	call salloc (paramlist, MAXPARAM*MAXLIST, TY_REAL)
	call salloc (filelist, (SZ_FNAME+1)*MAXLIST, TY_CHAR)
	call salloc (banderr, nwave, TY_REAL)

	# Get list of component throughput files

#	call searchgraf (verbose, graphtab,       # calling sequence change VGL 9/28/01
        call searchgraf (verbose, graphtab, GRF_COMPID,   # Calling sequence change VGL 9/28/01
			 comptab, mode, MAXLIST, 
			 MAXPARAM, SZ_FNAME, ncomp, Memi[nparam], 
			 Memr[paramlist], Memc[filelist])

	# Multiply component throughputs together for grand throughput

	call getthruput (SZ_FNAME, ncomp, Memc[filelist], MAXPARAM, 
			 Memi[nparam], Memr[paramlist], nwave, wave,
			 band, Memr[banderr])

	call sfree (sp)
end
