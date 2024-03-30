include	<tbset.h>
include	"../limit.h"

define	SZ_WAVE		1200

# LISTWAVE -- Compute a wavelength array for a list of tables

procedure listwave (wavtab, list, depcol, wave, nwave)

char	wavtab[ARB]	# i: wavelength table name
pointer	list		# i: list descriptor
char	depcol		# i: dependent column name
pointer	wave		# o: wavelength array (must be freed by caller)
int	nwave		# o: number of wavelengths
#--
bool	first, logspace
pointer	sp, tp, wv, dep, table
real	wmin, wmax, minwave, maxwave

string	wavecol   "WAVELENGTH"
data	logspace  / true /

int	tbpsta(), nxtlist()
pointer	tbtopn(), opnsyntab()

begin
	# Allocate memory for temporaray array

	call smark (sp)
	call salloc (table, SZ_FNAME, TY_CHAR)

	# If no table name, get wavelength set from the list
	# Otherwise, read it from the table

	if (wavtab[1] != EOS) {
	    tp = tbtopn (wavtab, READ_ONLY, NULL)
	    nwave = tbpsta (tp, TBL_NROWS)
	    call malloc (wave, nwave, TY_REAL)

	    call rdwave (tp, nwave, Memr[wave])
	    call tbtclo (tp)


	} else {
	    nwave = 0
	    wave = NULL

	    # Take union of wavelength ranges of each table

	    first = true
	    while (nxtlist (list, Memc[table], SZ_FNAME) != EOF) {
		tp = opnsyntab (Memc[table])
		call syncolptr (tp, wavecol, 1, wv)
		call syncolptr (tp, depcol, 2, dep)

		call wavelimits (tp, wv, dep, wmin, wmax)

		if (first) {
		    first = false
		    nwave = SZ_WAVE
		    call malloc (wave, nwave, TY_REAL)

		    minwave = wmin
		    maxwave = wmax
		} else {
		    minwave = min (wmin, minwave)
		    maxwave = max (wmax, maxwave)
		}
	    }

	    # Compute wavelength set from range

	    if (! first)
		call waveset (logspace, minwave, maxwave, nwave, Memr[wave])

	}

	call sfree (sp)
end
