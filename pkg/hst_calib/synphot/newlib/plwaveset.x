include	"../limit.h"
define	SZ_WAVE		1200

#* HISTORY *
#* B.Simon	 7-Jul-94	original

# PLWAVESET -- Get wavelength set for plotting tasks

procedure plwaveset (wavtable, grftable, cmptable, limit, 
		     pcode, ncode, maxcode,  wave, nwave)

char	wavtable[ARB]		# i: Table containing wavelength array
char	grftable[ARB]		# i: Instrument graph table
char	cmptable[ARB]		# i: Component name table
real	limit[4]		# i: Plot limits
int	pcode[maxcode,ncode]	# i: Expression pseudocode arrays
int	ncode			# i: Number of pseudocode arrays
int	maxcode			# i: dimension of pseudocode array
pointer	wave			# o: wavelength set (must be freed by caller)
int	nwave			# o: length of wavelength set
#--
bool	logspace

data	logspace  / true /

begin
	# Calculate wavelength set from plot limits if there is no 
	# wavelength set and plot limits are specified. Otehrwise,
	# compute the standard way

	if (wavtable[1] == EOS && ! IS_INDEFR (limit[LEFT]) && 
	    ! IS_INDEFR(limit[RIGHT])) {

	    nwave = SZ_WAVE
	    call malloc (wave, nwave, TY_REAL)

	    call waveset (logspace, limit[LEFT], limit[RIGHT], 
			  nwave, Memr[wave])

	} else {
	    call getwavelen (wavtable, grftable, cmptable, pcode, 
			     ncode, maxcode, wave, nwave)
	}

end

