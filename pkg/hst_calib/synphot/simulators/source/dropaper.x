include "simtwo.h"

#* HISTORY *
#* B.Simon	01-Nov-95	original

# DROPAPER -- Remove aperture thruput from combined thruput

procedure dropaper (obsmode, graphtab, comptab, nwave, wave, thruput)

char	obsmode[ARB]	# i: observation mode string
char	graphtab[ARB]	# i: graph table name
char	comptab[ARB]	# i: component lookup table name
int	nwave		# i: length of wavelength array
real	wave[ARB]	# i: wavelengths at which thruput are computed
real	thruput[ARB]	# u: total thruput
#--
string	keyword  APERKEY

real    setzero() # Zero replacement for division procedure 
extern	setzero
pointer	sp, aperthru, apertab

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (aperthru, nwave, TY_REAL)
	call salloc (apertab, SZ_FNAME, TY_CHAR)

	# Find the aperture thruput table

	call findthruput (keyword, obsmode, graphtab, comptab, 
			  Memc[apertab], SZ_FNAME)

	# Divide out the aperture thruput

	if (Memc[apertab] != EOS) {
	    call rdband (Memc[apertab], nwave, wave, Memr[aperthru])
	    call advzr (thruput, Memr[aperthru], thruput, nwave, setzero)
	}

	call sfree (sp)
end
