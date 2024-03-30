include <tbset.h>
include <synphot.h>

define	MAXWAVE		1200

#* HISTORY *
#* B.Simon	24-May-93	original

# EVALBAND -- Compute the total throughput for an observation mode

procedure evalband ()

#--
pointer	obsmode		# observation mode
pointer	output		# output table name
pointer	wavtable	# wavelength table name
pointer	grtbl		# graph table name
pointer	cmptbl		# component lookup table name
real	area		# telescope area

bool	logspace
int	nwave
pointer	tp, sp, wave, thruput, thruerr
real	phot[5]

data	logspace  / true /

bool	isblank()
int	tbpsta()
pointer	tbtopn()
real	clgetr(), funit(), pivlam(), rmslam()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (obsmode, SZ_FNAME, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (wavtable, SZ_FNAME, TY_CHAR)
	call salloc (grtbl, SZ_FNAME, TY_CHAR)
	call salloc (cmptbl, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgstr ("obsmode", Memc[obsmode], SZ_FNAME)
	call clgstr ("output", Memc[output], SZ_FNAME)
	call clgstr ("wavetab", Memc[wavtable], SZ_FNAME)
	call clgstr ("grtbl", Memc[grtbl], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptbl], SZ_FNAME)
	area = clgetr ("area")

	# If the wave table is blank, call getbandx() which does calculates
	# the wavelength set from the throughput tables. If it is not blank,
	# read the table and call evalbandx() which uses the input wavelength
	# set.

	if (isblank (Memc[wavtable])) {
	    nwave = MAXWAVE

	    call salloc (wave, nwave, TY_REAL)
	    call salloc (thruput, nwave, TY_REAL)
	    call salloc (thruerr, nwave, TY_REAL)

	    call getbandx (Memc[obsmode], Memc[grtbl], Memc[cmptbl], logspace,
			   nwave, Memr[wave], Memr[thruput], Memr[thruerr])

	} else {
	    tp = tbtopn (Memc[wavtable], READ_ONLY, NULL)
	    nwave = tbpsta (tp, TBL_NROWS)

	    call salloc (wave, nwave, TY_REAL)
	    call salloc (thruput, nwave, TY_REAL)
	    call salloc (thruerr, nwave, TY_REAL)

	    call rdwave (tp, nwave, Memr[wave])
	    call tbtclo (tp)

	    call evalbandx (Memc[obsmode], nwave, Memr[wave], Memc[grtbl], 
			    Memc[cmptbl], Memr[thruput], Memr[thruerr])
	}

	# Calculate the photometric paraneters from the throughput

	phot[1] = funit (area, nwave, Memr[wave], Memr[thruput])
	phot[2] = STZERO
	phot[3] = pivlam (nwave, Memr[wave], Memr[thruput])
	phot[4] = rmslam (nwave, Memr[wave], Memr[thruput])
	phot[5] = area

	# Write results to the output table

	call putband (Memc[output], Memc[obsmode], Memc[grtbl], Memc[cmptbl],
		      nwave, Memr[wave], Memr[thruput], Memr[thruerr], phot)

	call sfree (sp)
end
