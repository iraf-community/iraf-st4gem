include	<gset.h>
include <tbset.h>
include "../plspec/plspec.h"

# RATIOPLOT -- Plot the ratio of synthetic spectra to observed spectra.
#
# Sep 1993 H.A.Bushouse: Resample sphot data to same waveset as spec data

procedure ratioplot( gp, nwave, wave, speclist, nspec, sphotlist, nsphot,
	             phratio, sigma, piv, fwhm, nphot, ebmv1, ebmv2, nebmv,
	             form )

int	gp			# i: pointer to graphics device
int	nwave			# i: number of wavelenths
pointer	wave			# i: pointer to wavelength array
char	speclist[SZ_FNAME,ARB]	# i: list of synthetic spectrum files
int	nspec			# i: number of spectra
char	sphotlist[SZ_FNAME,ARB] # i: list of spectrophotometry files
int	nsphot			# i: number of spectrophotometry spectra
real	phratio[ARB]		# i: ratio of photometry/synphot
real	sigma[ARB]		# i: 1-sigma error bars for ratio
real	piv[ARB]		# i: pivot wavelengths for ratio
real	fwhm[ARB]		# i: fwhms for ratio
int	nphot			# i: number of ratio values
real	ebmv1			# i: first ebmv value
real	ebmv2			# i: last ebmv value
int	nebmv			# i: number of ebmv values
char	form[ARB]		# i: form of the data


int	isphot, ispec, iebmv, iw, ic, nspwave
int	strsearch()
real	debmv, ebmv, temp
pointer sp, spec, spdat, spsig, spfwhm, spform, spmode, spstar, cmd
pointer sptwv, sptdt, sptsg
pointer	inform, xdat, ydat, xfwhm, ysig, grtbl, cmptbl

bool	status, mag

# Jan 1991  DB  Add plotting of error in ratio
# Nov 1993 Bernie Simon; build ebmv string with applyebmv()

begin
	# Allocate memory
	call smark( sp )
	call salloc( spec, nwave, TY_REAL)
	call salloc( spdat, nwave, TY_REAL)
	call salloc( spsig, nwave, TY_REAL)
	call salloc( xdat, nwave, TY_REAL )
	call salloc( ydat, nwave, TY_REAL )
	call salloc( xfwhm, nwave, TY_REAL )
	call salloc( ysig, nwave, TY_REAL )
	call salloc( sptwv, MAXWAVE, TY_REAL)
	call salloc( sptdt, MAXWAVE, TY_REAL)
	call salloc( sptsg, MAXWAVE, TY_REAL)
	call salloc( spfwhm, MAXWAVE, TY_REAL)
	call salloc( spform, SZ_COLUNITS, TY_CHAR)
	call salloc( spmode, SZ_LINE, TY_CHAR)
	call salloc( spstar, SZ_FNAME, TY_CHAR)
	call salloc( cmd, SZ_LINE, TY_CHAR)
	call salloc( inform, SZ_COLUNITS, TY_CHAR)
	call salloc( grtbl, SZ_FNAME, TY_CHAR)
	call salloc( cmptbl, SZ_FNAME, TY_CHAR)

	# Get graph and component table names
	call clgstr( "grtbl", Memc[grtbl], SZ_FNAME)
	call clgstr( "cmptbl", Memc[cmptbl], SZ_FNAME)

	# Check for magnitudes in form string
	mag = false
	if (strsearch(form,"mag") > 0 || strsearch(form,"MAG") > 0 )
	   mag = true

	debmv = (ebmv2 - ebmv1) / max(1,nebmv-1)
	do isphot = 1, nsphot {
	   nspwave = MAXWAVE
	   call loadsphot( nspwave, sphotlist[1,isphot], Memr[sptwv], 
			   Memr[sptdt], Memr[sptsg], Memr[spfwhm], Memc[spform],
	                   Memc[spmode], Memc[spstar] )

	   # Resample the spectrum and sigmas
	   call linterp (Memc[spform], nspwave, Memr[sptwv], Memr[sptdt],
			 nwave, Memr[wave], Memr[spdat])
	   call linterp (Memc[spform], nspwave, Memr[sptwv], Memr[sptsg],
			 nwave, Memr[wave], Memr[spsig])

	   call sphotform( nwave, Memr[wave], Memr[spdat], Memr[spsig], 
	                   Memc[spform], Memr[spdat], Memr[spsig], form )

	   # Loop over spectra
	   do ispec = 1, nspec {

	      # Loop over E(B-V) values
	      do iebmv = 1, max(1,nebmv) {
	         ebmv = ebmv1 + (iebmv-1)*debmv
		 call applyebmv (speclist[1,ispec], Memc[cmd], ebmv, SZ_LINE)

	         iw = 1
	         call compspec ( Memc[cmd], iw, Memc[grtbl], Memc[cmptbl],
	                         nwave, wave, spec, Memc[inform] )
	         call specform( nwave, Memr[wave], Memr[spec], Memc[inform],
	                        Memr[spec], form, status )

	         # Compute ratio for plotting
	         if ( mag ) {

	            do ic = 1, nwave {
	               if ( Memr[spec+ic-1] < INDEFR 
	                                    && Memr[spdat+ic-1] < INDEFR )
	                  Memr[spec+ic-1] = Memr[spdat+ic-1] - Memr[spec+ic-1]
	               else
	                  Memr[spec+ic-1] = INDEFR
	               Memr[ysig+ic-1] = Memr[spsig+ic-1]
	               Memr[xdat+ic-1] = Memr[wave+ic-1]
	               #Memr[xfwhm+ic-1] = Memr[spfwhm+ic-1]
	            }        

	         } else {
	            do ic = 1,nwave {

	               temp = Memr[spec+ic-1]
	               if ( Memr[spec+ic-1] < INDEFR
	                                    && Memr[spdat+ic-1] < INDEFR
	                                    && Memr[spec+ic-1] != 0. )
	                  Memr[spec+ic-1] = Memr[spdat+ic-1] / Memr[spec+ic-1]
	               else
	                  Memr[spec+ic-1] = INDEFR

	               if ( Memr[spsig+ic-1] < INDEFR
	                                     && temp < INDEFR && temp != 0. )
	                  Memr[ysig+ic-1] = Memr[spsig+ic-1] / temp
	               else
	                  Memr[ysig+ic-1] = INDEFR
	               Memr[xdat+ic-1] = Memr[wave+ic-1]
	               #Memr[xfwhm+ic-1] = Memr[spfwhm+ic-1]
	            }
	         }

	         # Plot the ratio
	         call gseti( gp, G_PLTYPE, GL_SOLID )
	         call gpline( gp, Memr[wave], Memr[spec], nwave )
	      }
	   }
	}

	if (nphot > 0 ) {
	   do ic = 1, nphot {
	      call gamove( gp, piv[ic] - 0.5*fwhm[ic], phratio[ic])
	      call gadraw( gp, piv[ic] + 0.5*fwhm[ic], phratio[ic])
	      if ( !IS_INDEFR (sigma[ic]) ) {
	         call gamove( gp, piv[ic], phratio[ic] + sigma[ic] )
	         call gadraw( gp, piv[ic], phratio[ic] - sigma[ic] )
	      }
	      call gmark( gp, piv[ic], phratio[ic], GM_CIRCLE, 1., 1. )
	   }
	}

	call sfree (sp)

end
