include	<tbset.h>

# ANGTOANY -- Convert wavelength units from angstroms

int procedure angtoany (units, wave, nwave)

char	units[ARB]	# i: wavelength units
real	wave[ARB]	# u: wavelength array
int	nwave		# i: number of wavelengths
#--
int	wtype, done, invert[13]
pointer	sp, uni
real	factor[13]

string	badunits "Unknown wavelength units"
string	unidic   "|angstroms|nanometers|microns|millimeters\
|centimeters|meters|hertz|kilohertz|megahertz|gigahertz|ev|kev|mev|"

data	factor	/ 1.0, 0.1, 1.0e-4, 1.0e-7, 1.0e-8, 1.0e-10, 2.9979e18,
		 2.9979e15, 2.9979e12, 2.9979e9, 1.2396e4, 1.2396e1,
	         1.2396e-2 /
data	invert	/ 6*NO, 7*YES /

int	strdic()
errchk	amulkr, arcpr

begin
	# Allocate memory for temporary string

	call smark (sp)
	call salloc (uni, SZ_COLUNITS, TY_CHAR)

	# Remove blanks and lower case units for case insensitive match

	call strcpy (units, Memc[uni], SZ_COLUNITS)
	call strfix (Memc[uni])

	# If no units found, set return flag to indicate this
	# Otherwise, find proper units for conversion and  do it
	
	if (Memc[uni] == EOS) {
	    done = NO

	} else {
	    wtype = strdic (Memc[uni], Memc[uni], SZ_COLUNITS, unidic)
	    if (wtype == 0)
		call synphoterr (badunits, units)

	    if (invert[wtype] == NO) {
		call amulkr (factor[wtype], wave, wave, nwave)
	    } else {
		call arcpr (factor[wtype], wave, wave, nwave)
	    }
	}

	call sfree (sp)
	return (done)
end
