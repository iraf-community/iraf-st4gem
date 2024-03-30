define	HPLANCK	6.62620E-27	# erg s
define	CLIGHT 2.997925E+18	# A/s
include <synphot.h>

# RDVEGA -- Read the spectrum of Vega in from an ascii file.  This is the
# spectrum used to produce vegamags.  This routine assumes that the vega
# spectrum in vega.h has units of mjy, it converts the spectrum to stmag

procedure rdvega( vwave, vega, nvega, vform )

real	vwave[ARB]	# o: vega wavelength set
real	vega[ARB]	# o: vega flux in mjy
int	nvega		# o: number of vega data points
char	vform[SZ_LINE]	# o: form of vegas spectrum that is output

# --

pointer	form, sp
int	ic
int	strsearch()
bool	status
real	f_one, f_two

string	name	"synphot$lib/vega.dat"
string	notfound "Cannot find vega spectrum file vega.dat in synphot$lib"
string	wrongform "Spectrum in vega.dat is not mjy, you better change it"

begin

	call smark( sp )
	call salloc( form, SZ_FNAME, TY_CHAR )

	nvega = 0
	# Blank the form string
	Memc[form] = EOS

	# Read in the vega spectrum form the vega.h file
	call ascspec( nvega, vwave, vega, name, Memc[form], status )

	if ( !status )
	   call error( 1, notfound )
	if (strsearch( Memc[form], "mjy" ) == 0 ) 
	   call error( 1, wrongform )

	# Convert mjy to photlam and photlam to stmag
	f_one = 1.e-26 / HPLANCK
	f_two = HPLANCK * CLIGHT
	do ic = 1, nvega {
	   if ( !IS_INDEFR (vega[ic]) && vega[ic] != 0. && vwave[ic] > 0. ) {
	      # To photlam
	      vega[ic] = f_one * vega[ic] / vwave[ic]
	      # To stmag
	      vega[ic] = -2.5 * alog10 (f_two * vega[ic]/vwave[ic]) + STZERO
	   } else
	      vega[ic] = INDEFR
	}
	call strcpy ( "stmag", vform, SZ_FNAME )

	call sfree( sp )

end
