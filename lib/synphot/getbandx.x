include "libsynphot.h"

# GETBANDX -- Calculate the wavelength and passband for an instrument mode

procedure getbandx (mode, graphtab, comptab, logspace, nwave, 
		    wave, thruput, thruerr)

char	mode[ARB]	# i: instrument mode
char	graphtab[ARB]	# i: graph table name
char	comptab[ARB]	# i: component lookup table name
bool	logspace	# i: use log spacing for wavelengths?
int	nwave		# i: number of wavelengths
real	wave[ARB]	# o: wavelength array
real	thruput[ARB]	# o: grand throughput
real	thruerr[ARB]	# o: throughput error
#--
bool	verbose
int	ncomp
pointer	sp, nparam, paramlist, filelist
real	minwave, maxwave

data	verbose	 / false /

errchk	searchgraf, waverange, getthruput

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

	# Compute wavelength set

	call waverange (SZ_FNAME, ncomp, Memc[filelist], minwave, maxwave)
	call waveset (logspace, minwave, maxwave, nwave, wave)

	# Multiply component throughputs together for grand throughput

	call getthruput (SZ_FNAME, ncomp, Memc[filelist], MAXPARAM, 
			 Memi[nparam], Memr[paramlist], nwave, 
			 wave, thruput, thruerr)

	call clssyntab
	call sfree (sp)
end
