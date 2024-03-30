include "libsynphot.h"

# EVALBANDX -- Calculate the passband for all components in the instrument path

procedure evalbandx (mode, nwave, wave, graphtab, comptab, thruput, thruerr)

char	mode[ARB]	# i: instrument mode
int	nwave		# i: length of wavelength and output arrays
real	wave[ARB]	# i: wavelengths at which output is computed
char	graphtab[ARB]	# i: graph table name
char	comptab[ARB]	# i: component lookup table name
real	thruput[ARB]	# o: grand throughput
real	thruerr[ARB]	# o: grand throughput error
#--
bool	verbose
int	ncomp
pointer	sp, nparam, paramlist, filelist

data	verbose	/ false /

errchk	searchgraf, getthruput

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (nparam, MAXLIST, TY_INT)
	call salloc (paramlist, MAXPARAM*MAXLIST, TY_REAL)
	call salloc (filelist, (SZ_FNAME+1)*MAXLIST, TY_CHAR)

	# Get list of component throughput files

	call inisyntab
#	call searchgraf (verbose, graphtab, 
        call searchgraf (verbose, graphtab, GRF_COMPID,   # Calling sequence change VGL 9/28/01
			 comptab, mode, MAXLIST, 
			 MAXPARAM, SZ_FNAME, ncomp, Memi[nparam], 
			 Memr[paramlist], Memc[filelist])

	# Multiply component throughputs together for grand throughput

	call getthruput (SZ_FNAME, ncomp, Memc[filelist], MAXPARAM, 
			 Memi[nparam], Memr[paramlist], nwave, wave,
			 thruput, thruerr)

	call clssyntab
	call sfree (sp)
end
