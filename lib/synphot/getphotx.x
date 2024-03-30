include "libsynphot.h"

# GETPHOTX -- Get the photometric parameters for an instrument mode

procedure getphotx (mode, graphtab, comptab, path, mxpath, phot)

char	mode[ARB]	# i: instrument mode
char	graphtab[ARB]	# i: graph table name
char	comptab[ARB]	# i: component lookup table name
char	path[ARB]	# o: list of component throughput tables
int	mxpath		# i: maximum length of path string
real	phot[4]		# o: photometric parameters
#--
bool	verbose, logspace
int	ncomp
pointer	sp, wave, thruput, thruerr, nparam, paramlist, filelist
real	minwave, maxwave

data	verbose	 / false /
data	logspace / true  /

errchk	searchgraf, waverange, getthruput

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (wave, LENWAVE, TY_REAL)
	call salloc (thruput, LENWAVE, TY_REAL)
	call salloc (thruerr, LENWAVE, TY_REAL)
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

	# Concatenate file names into path string

	call wpath (SZ_FNAME, ncomp, Memc[filelist], mxpath, path)

	# Compute wavelength set

	call waverange (SZ_FNAME, ncomp, Memc[filelist], minwave, maxwave)
	call waveset (logspace, minwave, maxwave, LENWAVE, Memr[wave])

	# Multiply component throughputs together for grand throughput

	call getthruput (SZ_FNAME, ncomp, Memc[filelist], MAXPARAM, 
			 Memi[nparam], Memr[paramlist], LENWAVE, 
			 Memr[wave], Memr[thruput], Memr[thruerr])

	# Calculate the photometric parameters

	call phopar (LENWAVE, Memr[wave], Memr[thruput], phot)

	call clssyntab
	call sfree (sp)
end
