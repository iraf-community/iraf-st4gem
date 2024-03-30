# FINDSPHOT -- Given a list of spectrophotometry files and a targetid, return
# the index of the spectrum that corresponds to the photometric obs.

procedure findsphot( photid, targetid, nsphot, isp )

char	photid[ARB]		# i: ID string of photometric observation
char	targetid[SZ_FNAME,ARB]	# i: List of spectra available
int	nsphot			# i: Number of spectra in targetid
int	isp			# o: Index of spectrum in targetid

# Dave Bazell First Code April, 1990
# Put error msg inside routine, May '90, DB
#--

int	ic, nchar
int	strsearch(), nowhite()
char	pat[SZ_LINE], str[SZ_LINE], errmsg[SZ_LINE]
string	photnotfound	"Targetid %s not found"

begin

	# Remove whitespace in pattern string
	nchar = nowhite( photid, str, SZ_LINE )

	# Loop through spectra and return index if found
	do ic = 1, nsphot {

	   # Remove whitespace in string being searched
	   nchar = nowhite( targetid[1,ic], pat, SZ_LINE )

	   # Look for targetid in spectrum name and vice versa.  This should
	   # be replaced with a more sophisticated algorithm
	   if ( strsearch( str, pat ) > 0 ) {
	      isp = ic
	      return
	   }
	}

	# targetid not found so signal error
	call sprintf( errmsg, SZ_LINE, photnotfound )
	    call pargstr( pat )
	call error( 1, errmsg )

end
