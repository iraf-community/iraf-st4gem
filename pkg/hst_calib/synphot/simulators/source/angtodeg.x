include	<math.h>
include	<tbset.h>

#* HISTORY *
#* B.Simon	17-Feb-95	original
#* B.Simon	26-Feb-95	added angtodegr

# ANGTODEGR -- Convert single precision angle to degrees

procedure angtodegr (units, angle)

char	units[ARB]	# i: units of angle
real	angle		# u: angle to be converted
#--
double	dangle

begin
	dangle = angle
	call angtodeg (units, dangle)
	angle = dangle

end

# ANGTODEG -- Convert angular measure to degrees

procedure angtodeg (units, angle)

char	units[ARB]	# i: units of angle
double	angle		# u: angle to be converted
#--
double	factor[5]
int	utype, invert[5]
pointer	sp, uni

string	badunits "Unknown angular units"
string	unidic   "|seconds|minutes|degrees|hours|radians|"

data	factor	/ 3600.0, 60.0, 1.0, 15.0, RADIAN /
data	invert  / YES,    YES,  NO,  NO,   NO     /

int	strdic()

begin
	# Allocate memory for temporary string

	call smark (sp)
	call salloc (uni, SZ_COLUNITS, TY_CHAR)

	# Remove blanks and lower case units for case insensitive match

	call strcpy (units, Memc[uni], SZ_COLUNITS)
	call strfix (Memc[uni])

	# If units found, find proper units for conversion and do it
	
	if (Memc[uni] != EOS) {
	    utype = strdic (Memc[uni], Memc[uni], SZ_COLUNITS, unidic)
	    if (utype == 0)
		call synphoterr (badunits, units)

	    if (invert[utype] == NO) {
		angle = angle * factor[utype]
	    } else {
		angle = angle / factor[utype]
	    }
	}

	call sfree (sp)
end

