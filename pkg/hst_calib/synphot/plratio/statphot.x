include "../plspec/dspec.h"
include "../plspec/dphot.h"

#* HISTORY *
#* B.Simon	14-Jul-94	original

# STATPHOT -- Compute goodness of fit statistics for photometric ratios

procedure statphot (form, tgtlist, ntgt, dspec, dphot, 
		    pchisq, pbias, prms, pstat)

char	form[ARB]	# i: form of data
pointer	tgtlist		# i: correspondence btw spectrum and photometry
int	ntgt		# i: length of tgtlist array
pointer	dspec		# i: spectrum descriptor
pointer	dphot		# i: photmetry descriptor
real	pchisq		# o: chi squared
real	pbias		# o: bias
real	prms		# o: root mean squared
int	pstat		# o: number of data points used in statistics
#--
int	itgt, jtgt
pointer	sp, rstim, rerror, sstim, pstim

begin
	# Initialize statistics to undefined

	pchisq = INDEFR
	pbias = INDEFR
	prms = INDEFR
	pstat = INDEFI

	# Return if we don't have the two arrays needed to compute statistics

	if (tgtlist == NULL || dspec == NULL || dphot == NULL)
	    return

	# Allocate arrays to hold sorted photometry data

	call smark (sp)
	call salloc (rstim, ntgt, TY_REAL)
	call salloc (rerror, ntgt, TY_REAL)

	# Fake photometry errors

	call amovkr (INDEFR, Memr[rerror], ntgt)

	sstim = SPC_STIM(dspec)
	pstim = PHT_STIM(dphot)

	# Sort photmetry data according to tgtlist

	do itgt = 1, ntgt {
	    jtgt = Memi[tgtlist+itgt-1]
	    if (jtgt == 0) {
		Memr[rstim+itgt-1] = INDEFR
	    } else {
		Memr[rstim+itgt-1] = Memr[pstim+jtgt-1]
	    }
	}

	# Compute goodness of fit statistics

	call goodfit (Memr[rstim], Memr[rerror], Memr[sstim], ntgt, form,
		      pchisq, pbias, prms, pstat)

	call sfree (sp)
end
