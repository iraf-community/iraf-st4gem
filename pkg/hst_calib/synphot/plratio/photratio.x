include "../plspec/dspec.h"
include "../plspec/dphot.h"

#* HISTORY *
#* B.Simon	14-Jul-94	original

# PHOTRATIO -- Compute ratio of effective stimulus of photometry to spectra

procedure photratio (form, tgtlist, ntgt, dspec, dphot)

char	form[ARB]	# i: form of effective stimulus
pointer	tgtlist		# i: correspondence btw spectra and photmetry
int	ntgt		# i: length of tgtlist array
pointer	dspec		# i: spectum descriptor
pointer	dphot		# u: photometric descriptor
#--
int	mag, itgt, jtgt, nstat
pointer	sstim, pstim, rstim

int	is_magunit()

begin
	if (dphot == NULL)
	    return

	# Allocate ratio array and initialize to INDEF

	pstim = PHT_STIM(dphot)
	nstat = PHT_NSTAT(dphot)

	call malloc (rstim, nstat, TY_REAL)
	call amovkr (INDEFR, Memr[rstim], nstat)

	# Take exit if target list or spectrum descriptor is undefined

	if (dspec == NULL || tgtlist == NULL) {
	    PHT_STIM(dphot) = rstim
	    call mfree (pstim, TY_REAL)
	    return
	}

	sstim = SPC_STIM(dspec)
	mag = is_magunit (form)

	do itgt = 1, ntgt {
	    jtgt = Memi[tgtlist+itgt-1]
	    if (jtgt == 0)
		next

	    if (IS_INDEFR(Memr[sstim+itgt-1]) || IS_INDEFR(Memr[pstim+itgt-1]))
		next

	    # Compute ratio of effective stimuli

	    if (mag == NO) {
		if (Memr[pstim+jtgt-1] > 0.0)
		    Memr[rstim+jtgt-1] = Memr[pstim+jtgt-1]/ Memr[sstim+itgt-1]
	    } else {
		Memr[rstim+jtgt-1] = Memr[pstim+jtgt-1] - Memr[sstim+itgt-1]
	    }
	}

	# Substitute photometric ratio for effective stimulus array

	PHT_STIM(dphot) = rstim
	call mfree (pstim, TY_REAL)

end
