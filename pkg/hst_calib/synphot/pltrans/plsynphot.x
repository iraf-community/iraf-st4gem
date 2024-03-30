include	<gset.h>
include	"pltrans.h"

# PLSYNPHOT -- Plot synthetic photometry from pltrans

procedure plsynphot( gp, xval, yval, nstar, symbol, maxpt)

pointer	gp		# i: Graphics pointer
real	xval[ARB]	# i: x-values of points to plot
real	yval[ARB]	# i: y-values of points to plot
int	nstar[ARB]	# i: Plot nstar points of a given symbol type
pointer	symbol[ARB]	# i: Pointers to symbols for point sets (char)
int	maxpt		# i: Maximum number of points that can be plotted

int	ic, ,jc, firstpt, npt, nchar
int	mark2int(), strdic(), unbang(), stridx()

char	sym[SZ_FNAME]

# Dictionarys or allowed lines and marks for plotting
string	lines "|solid|dashed|dotted|dotdash|"
string	marks "|point|fill|box|plus|cross|diamond|hline|vline|hebar|vebar|circle|"

begin

	npt = 0
	for (ic = 1; nstar[ic] > 0 && ic <= MAXLIST; ic = ic + 1 ) {

	   firstpt = npt + 1
	   npt = min( npt + nstar[ic], maxpt)

	   # Return if we have plotted maxpt points
	   if ( firstpt > maxpt )
	      return

	   call strcpy( Memc[symbol[ic]], sym, SZ_FNAME)

	   # Plot an alphanumeric if escaped with a bang !.
	   if ( stridx("!",Memc[symbol[ic]]) > 0 ) {
	      nchar = unbang(Memc[symbol[ic]],  Memc[symbol[ic]], SZ_FNAME)
	      do jc = firstpt, npt
	         call gtext( gp, xval[jc], yval[jc], Memc[symbol[ic]],
	                    "h=c;v=c")
	      next
	   }

	   call strlwr( Memc[symbol[ic]] )

	   # Plot set as line if symbol found in line_dict
	   if( strdic( Memc[symbol[ic]], Memc[symbol[ic]], 
	              SZ_FNAME, lines) > 0 )  {

	      call gseti( gp, G_PLTYPE, mark2int(Memc[symbol[ic]]) )
	      call gpline( gp, xval[firstpt], yval[firstpt], npt )

	   }

	   # Plot points if symbol found in mark_dict
	   else if ( strdic( Memc[symbol[ic]], Memc[symbol[ic]],
	             SZ_FNAME, marks ) > 0 )
	      call gpmark( gp, xval[firstpt], yval[firstpt], npt,
	                   mark2int(Memc[symbol[ic]]), 1., 1.)

	   # Signal error and plot as point
	   else
	      do jc = firstpt, npt
	         call gtext( gp, xval[jc], yval[jc], sym, "h=c;v=c")

	}

end
