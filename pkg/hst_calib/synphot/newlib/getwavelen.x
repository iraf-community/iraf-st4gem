include <mach.h>
include	<tbset.h>
include	<synphot.h>

define  VEGA            "st4gem$lib/synphot/vega.dat"

#* HISTORY *
#* B.Simon	29-Apr-94	original
#* B.Simon	26-May-94	support multiple pseudocode arrays
#* B.Simon	01-Aug-94	moved wavemin,wavemax to synphot library
#* B.Simon	20-Jul-95	added call to getnwave
#* B.Simon	25-Jun-96	added nmax calculation
#* B.Simon	27-Jun-96	intersect wavelength set with Vega's
#* B.Simon	06-Aug-96	removed nmax, added call to crosswave
#* B.Simon	20-Sep-99	added calculation to assure min spacing
#* B.Simon	21-Aug-00	modified to call synwave

# GETWAVELEN -- Get the wavelength set from a table or synphot expression

procedure getwavelen (wavtable, grftable, cmptable, pcode, ncode, maxcode, 
		      wave, nwave)

char	wavtable[ARB]		# i: Table containing wavelength array
char	grftable[ARB]		# i: Instrument graph table
char	cmptable[ARB]		# i: Component name table
int	pcode[maxcode,ncode]	# i: Expression pseudocode arrays
int	ncode			# i: Number of pseudocode arrays
int	maxcode			# i: dimension of pseudocode array
pointer	wave			# o: wavelength set (must be freed by caller)
int	nwave			# o: length of wavelength set
#--
int	intersect, status
pointer	tp
real	minwave, maxwave

data	intersect  / NO /
string	badvega  "Error in Vega spectrum"
string	nocross  "Wavelength table does not intersect Vega spectrum"

int	tbpsta()
pointer	tbtopn()

begin
	# If no table name, get wavelength range from synphot expression
	# Otherwise, read it from the table

	if (wavtable[1] == EOS) {
	    call synwave (intersect, pcode, ncode, maxcode, grftable, 
			  cmptable, wave, nwave)
	} else {
	    # Read spectrum from table

	    tp = tbtopn (wavtable, READ_ONLY, NULL)
	    nwave = tbpsta (tp, TBL_NROWS)
	    call malloc (wave, nwave, TY_REAL)

	    call rdwave (tp, nwave, Memr[wave])
	    call tbtclo (tp)

	    # Take intersection of range with range of Vega's spectrum

	    call tabrange (VEGA, minwave, maxwave, status)
	    if (status == ERR)
		call printerr_str (badvega, VEGA)

	    call crosswave (minwave, maxwave, nwave, Memr[wave])
	}

	if (nwave < 1)
	    call printerr_str (nocross, wavtable)
end
