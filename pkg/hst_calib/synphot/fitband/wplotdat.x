include	"../plspec/plspec.h"

# WPLOTDAT -- Write a fitband residuals file for use by pltrans

procedure wplotdat( npar, par )

int	npar		# i: number of parameters in model
real	par[ARB]	# i: model parameter values
#--

int	nbpt, nband, nwave, nphot, nsphot, ndat, ip, ic, isp, npt
int	offset
pointer	sp, wavetab, wave, colunits, grtbl, cmptbl, model, xmode, xform
pointer	mode1, mode2, band, dat, sig, xdat, bias, chi2, plotfile, passband
real	piv
char	photlist[SZ_FNAME,MAXPHOT], pform[SZ_FNAME,MAXPHOT]
char	pmode[SZ_LINE,MAXPHOT], pstar[SZ_FNAME,MAXPHOT]
char	sphotlist[SZ_FNAME,MAXSP]
char	diff[4]
bool	streq()

string	nulcmd	""
string	toomuchphot	"Warning: Too much photometry data in %s\n"

include "../lib/amoebafit.h"
include "targets.h"

begin

	# Allocate memory
	call smark( sp )
	call salloc( wavetab, SZ_FNAME, TY_CHAR )
	call salloc( grtbl, SZ_FNAME, TY_CHAR )
	call salloc( cmptbl, SZ_FNAME, TY_CHAR )
	call salloc( model, SZ_LINE, TY_CHAR )
	call salloc( passband, SZ_LINE, TY_CHAR )
	call salloc( xmode, SZ_FNAME, TY_CHAR )
	call salloc( xform, SZ_FNAME, TY_CHAR )
	call salloc( mode1, SZ_LINE, TY_CHAR )
	call salloc( mode2, SZ_LINE, TY_CHAR )
	call salloc( band, 1, TY_REAL )
	call salloc( dat, MAXPHOT, TY_REAL )
	call salloc( sig, MAXPHOT, TY_REAL )
	call salloc( xdat, MAXPHOT, TY_REAL )
	call salloc( bias, MAXPHOT, TY_REAL )
	call salloc( chi2, MAXPHOT, TY_REAL )
	call salloc( plotfile, SZ_FNAME, TY_CHAR )

	# Get wavelength grid and reference table names for synphot
	call clgstr( "wavetab", Memc[wavetab], SZ_FNAME )
	call getwave( Memc[wavetab], nwave, wave, colunits)
	call clgstr( "grtbl", Memc[grtbl], SZ_FNAME)
	call clgstr( "cmptbl", Memc[cmptbl], SZ_FNAME )

	# Fill the model string with current parameters for evaluation
	call strcpy ( fitmodel, Memc[model], SZ_LINE )
	call insertpar( par, Memc[model] )

	# Evaluate the model and return bias for plotting
	call bandchi2( Memc[model], fitdata, Memr[chi2], Memr[bias], nbpt )

	# Get xmode and xform for evaluation of x data points
	call clgstr("xmode", Memc[xmode], SZ_LINE)
	call clgstr("xform", Memc[xform], SZ_FNAME)
	call clgstr("plotfile", Memc[plotfile], SZ_FNAME )

	call splitmode( Memc[xmode], nband, Memc[mode1], Memc[mode2] )
	if ( nband == 1 ) {
	   call synphot( nwave, wave, Memr[band], nulcmd, "pivlam", Memc[xmode],
	                    Memc[grtbl], Memc[cmptbl], piv)
	}

	# Parse the data string containing the photometry file name
	call parsedat( fitdata, sphotlist, nsphot, MAXSP, 
		       photlist, nphot, MAXPHOT )

	# Loop over photometry data and calculate xmode.  offset is the
	# position of the next empty element in the pstar array.  This lets
	# us put all the star names in one array
	offset = 1
	npt = 0
	do ip = 1, nphot {
	   ndat = MAXPHOT
	   call loadphot( ndat, photlist[1,ip], Memr[dat], Memr[sig], pform,
	                  pmode, pstar[1,offset] )
	   offset = offset + ndat
	   if ( ndat >= MAXPHOT ) {
	      call printf( toomuchphot)
	         call pargstr( photlist[1,ip] )
	   }

	   # Loop over each data point and calculate the photometry
	   # in the desired mode
	   do ic = 1, ndat {

	      npt = npt + 1
	      # If the targetid = "none" then assume that the id in pstar[]
	      # is the file name
	      if ( streq( targetid[1,1], "none" ) ) {

	         call findsphot( pstar[1,ic], pstar, ndat, isp )
	         call strcpy( pstar[1,isp], specfile[1,isp], SZ_FNAME )

	      } else
	         call findsphot( pstar[1,ic], targetid, ntarget, isp )

	      # Replace 'obsmode' with current mode if necessary
	      call insertmode( pmode[1,ic], Memc[model], Memc[passband] )
	      call synphot(nwave, wave, Memr[band], specfile[1,isp],
	                   pform[1,ic], Memc[passband], Memc[grtbl],
	                   Memc[cmptbl], Memr[xdat+npt-1])

	      # Use same sigma array for both x and y data
	      Memr[sig+npt-1] = INDEFR

	      if ( nband == 1 ) {
	         call sphotform( 1, piv, Memr[xdat+npt-1], Memr[sig+npt-1],
	                         pform[1,ic], Memr[xdat+npt-1], 
	                         Memr[sig+npt-1], Memc[xform] )

	      }
	   }
	}

	call strcpy( "diff", diff, 4 )
	call savesf(Memc[plotfile], npt, Memr[xdat], Memr[sig], 
	            Memc[xform], Memc[xmode], Memr[bias], Memr[sig], 
	            pform[1,1], pmode[1,1], pstar, diff)

	call sfree( sp )	
end
