###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	procedure decode_plot ( ranges )
#		char	ranges[ARB]
#
#
#  Description:	DECODE_PLOT decodes the plot numbers
#
#
#
#  Arguments:	char	ranges[ARB]	- String describing data ranges
#					  to be used in the fit
#
#  Returns:	none
#
#  Notes:	Shares data in "specfit.com"
#
#  History:	June	1994	Grimes
#
###########################################################################

include	"specfit.h"

define	COMMA	','
define	DASH	'-'
define	STAR	'*'

procedure decode_plot (ranges)
char	ranges[ARB]

char	subrange[SZ_FNAME], numstr[SZ_FNAME]
int	i, j, nreg, slen, tlen, ip, temp, temp2, k
int	cstart, idash, nc
int	comma[MAXSAMPLE]

int	strlen(), ctoi()

include "specfit.com"

begin
	for ( i = 1; i<=ncomp; i = i + 1)
		inplot[i] = NO

	tlen = strlen(ranges)

# If any character is STAR, then all data included.
	for ( i = 1; i <= tlen; i = i + 1) {
		if ( ranges[i] == STAR ) {
			for ( j = 1; j <= ncomp; j = j + 1)
				{ inplot[j] = YES }
			return
		}
	}
	
# Count number of commas to get number of sample regions to decode.
# Last entry in comma should always be the EOS location.
	nreg = 1
	for ( i = 1; i <= tlen; i = i + 1) {
		if ( ranges[i] == COMMA ) {
			comma[nreg] = i
			nreg = nreg + 1
		}
	}
	comma[nreg] = tlen + 1

# Loop through subregions to decode each one
	cstart = 1
	for ( i = 1; i <= nreg; i = i + 1) {

		slen = comma[i] - cstart
		call strcpy(ranges[cstart], subrange, slen)
		cstart = comma[i] + 1

#Look for a DASH inside the subrange.  idash = 0 means no dash.
		idash = 0
		for ( j = 1; j < slen; j = j + 1)  {
			if (  subrange[j] == DASH )  idash = j
		}
# idash = 0 means no dash.  Decode the whole string and define the
# sample region to be redundantly inclusive.
		if ( idash == 0 ) {
			ip = 1
			nc = ctoi( subrange, ip, temp )
			inplot[temp] = YES
		} else {

# Handle the more typical case with a dash
			call strcpy(subrange, numstr, idash-1)
			ip = 1
			nc = ctoi( numstr, ip, temp )
			call strcpy(subrange[idash+1], numstr, slen-idash)
			ip = 1
			nc = ctoi( numstr, ip, temp2 )
			for ( k = temp; k <= temp2; k = k + 1) 
				{ inplot[k] = YES }
		}
	}

# Not consolidating adjacent regions for now, so nsamp = nreg
	nsamp = nreg

# Call setranges to re-define the "infit" array
	call setranges()

end

