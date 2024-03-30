include "dband.h"

#* HISTORY *
#* B.Simon	08-Jun-94	original

# GETDBAND -- Get bandpasses to plot and save in descriptor

procedure getdband (obsmode, norm, wavtab, graftab, comptab, dband)

char	obsmode[ARB]	# i: Instrument observation mide
bool	norm		# i: Normalize passbands to one?
char	wavtab[ARB]	# i: Table containing wavelength array
char	graftab[ARB]	# i: Instrument graph table
char	comptab[ARB]	# i: Component name table
pointer	dband		# o: bandpass descriptor
#--
int	nband, nwave, degree
pointer	sp, mode, ocode, optr, wave, band
real	factor

string	notband  "Not a passband"

int	numlist(), nxtlist()
pointer	rdlist()
real	ahivr()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (mode, SZ_LINE, TY_CHAR)
	call salloc (ocode, SZ_LINE, TY_INT)

	# Initialize table cache and synphot variables

	call inisyntab
	call undefsynvar

	# Loop over each observation mode

	dband = NULL
	optr = rdlist (obsmode)

	while (nxtlist (optr, Memc[mode], SZ_LINE) != EOF) {
	    # Compile obsmode passband into pseudocode

	    call expcompile (Memc[mode], Memi[ocode], SZ_LINE)

	    # Allocate strcture to hold passband

	    if (dband == NULL) {
		call getwavelen (wavtab, graftab, comptab, Memi[ocode], 
				 1, SZ_LINE, wave, nwave)

		nband = numlist (optr)
		call malloc (band, nwave*nband, TY_REAL)
		
		call malloc (dband, LEN_BNDSTRUCT, TY_STRUCT)
		BND_NWAVE(dband) = nwave
		BND_NBAND(dband) = nband
		BND_WAVE(dband) = wave
		BND_BAND(dband) = band
	    }

	    # Calculate passband

	    call syncalc (Memi[ocode], SZ_LINE, NULL, nwave, Memr[wave], 
			  graftab, comptab, Memr[band], degree)

	    if (degree != 0)
		call printerr_str (notband, Memc[mode])

	    # Normalize to range of [0, 1] if requested by user

	    if (norm) {
		factor = ahivr (Memr[band], nwave)
		if (factor > 0.0) {
		    factor = 1.0 / factor
		    call amulkr (Memr[band], factor, Memr[band], nwave)
		}
	    }

	    band = band + nwave
	}

	# Close files and release memory

	call clssyntab
	call freelist (optr)

	call sfree (sp)
end
