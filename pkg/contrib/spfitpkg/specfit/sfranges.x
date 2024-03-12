###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	procedure decode_ranges( ranges )
#		char	ranges[ARB]
#
#		procedure setranges()
#
#  Description:	DECODE_RANGES decodes the sample ranges defined by the ASCII
#		character string input from the parameter file, and then
#		calls setranges to implement them.
#
#		SETRANGES uses the sample ranges defined in the common blocks
#		to set flags of 1 or 0 in the array infit declaring whether a
#		particular data point is to be included in the fit or not.
#
#  Arguments:	char	ranges[ARB]	- String describing data ranges
#					  to be used in the fit
#
#  Returns:	none
#
#  Notes:	Shares data in "specfit.com"
#
#  History:	June	1989	Gerard Kriss
#		7/4/95		J Grimes  	changed sample range to
#						be real valued not integral
#               7/20/00         Ed Colbert      changed definition for ctor()
#                               for Knox Long   from real to int
###########################################################################

include	"specfit.h"

define	COMMA	','
define	DASH	'-'
define	STAR	'*'

procedure decode_ranges(ranges)
char	ranges[ARB]

char	subrange[SZ_FNAME], numstr[SZ_FNAME]
int	i, j, nreg, slen, tlen, ip
int	cstart, idash, nc
int	comma[MAXSAMPLE]

int	strlen()
int	ctor()

include "specfit.com"

begin
	tlen = strlen(ranges)

# If any character is STAR, then all data included.
	for ( i = 1; i <= tlen; i = i + 1) {
		if ( ranges[i] == STAR ) {
			nsamp = 1
			sample[1, 1] = Memr[lambdas]
			sample[2, 1] = Memr[lambdas+npts-1] + 1.
			call setranges()
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
			nc = ctor( subrange, ip, sample[1,i] )
			sample[2, i] = sample[1, i]
		} else {

# Handle the more typical case with a dash
			call strcpy(subrange, numstr, idash-1)
			ip = 1
			nc = ctor( numstr, ip, sample[1,i] )
			call strcpy(subrange[idash+1], numstr, slen-idash)
			ip = 1
			nc = ctor( numstr, ip, sample[2,i] )
		}

	}

# Not consolidating adjacent regions for now, so nsamp = nreg
	nsamp = nreg


# Call setranges to re-define the "infit" array
	call setranges()

end


procedure setranges()
int	i, j

include "specfit.com"

begin
	nfitpts = 0
	for ( i = 0; i < npts; i = i + 1) {
	  Memi[infit+i] = 0
	  for ( j = 1; j <= nsamp; j = j + 1) {
	    if((Memr[lambdas+i]>=sample[1,j])&&(Memr[lambdas+i]<=sample[2,j])) {
		Memi[infit+i] = 1
		nfitpts = nfitpts + 1
	    }
	  }
	}
end


