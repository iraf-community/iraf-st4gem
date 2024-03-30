# VALIDPT -- Searches for a valid data point for PLTRANS.

procedure validpt( starname, obsmode, plmode, colordiff, thisstar,
                   ridx, nrow, ic, validx, npt)

char	starname[SZ_FNAME,ARB]	# i: Array of starnames
char	obsmode[SZ_LINE,ARB]	# i: Array of corresponding obsmodes
char	plmode[ARB]		# i: Obsmode for plotting in pltrans
bool	colordiff		# i: Flag indicating if plmode is col.diff.
char	thisstar[ARB]		# i: Name of current star
int	ridx[ARB]		# i: Row indexs from sort routine
int	nrow			# i: Number of rows in ridx
int	ic			# i: Index of 1st occurence of thisstar in 
				#    starname
int	validx[2,ARB]		# o: Index array of valid points
int	npt			# o: Number of valid points

int	jc, ipos1, ipos2
int	strsearch()

bool	streq()

begin

	# Loop over first stars
	while( streq( obsmode[1,ridx[ic]], thisstar ) ) {

	   # Search plmode for obsmode of first star

	   ipos1 = strsearch( plmode, obsmode[ridx[ic],1] )
	   if( ipos1 > 0 ) {

	      # If plmode is a color difference...
	      if( colordiff ) {

	         # Loop over same starnames to find a valid second star
	         jc = ic + 1
	         while ( streq(starname[ridx[jc],1], thisstar) ) {

	            # Search xmode for mode of second star
	            ipos2 = strsearch( plmode, obsmode[ridx[jc],1] )

	            if ( ipos2 > 0 ) {

	               if ( ipos2 > ipos1 ) {		# ic - jc
	                  validx[1,npt] = ridx[ic]
	                  validx[2,npt] = ridx[jc]
	               } else if ( ipos2 < ipos1 ) {	# jc - ic
	                  validx[1,npt] = ridx[jc]
	                  validx[2,npt] = ridx[ic]
	               }

	            }
	            npt = npt + 1
	            jc = jc + 1

	         }

	      } else {

	         validx[1,npt] = ridx[ic]
	         validx[2,npt] = 0
	         npt = npt + 1

	      }

	   # Shouldn't get here since only valid obsmodes are supposed to
	   # be sent in
	   } else
	      call printf("Invalid obsmode for star1.\n")

	   ic = ic + 1
	}

end
