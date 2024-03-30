include	<gset.h>
include "plspec.h"

# PHOTPLOT -- Plot photometry for PLSPEC
#
# Sep 1993 Howard Bushouse; Add warning for 2 mode photometry data

procedure photplot( gp, wave, nwave, form, photname, grtbl, cmptbl, ymin, ymax)

pointer	gp		# i: Graphics pointer
pointer	wave		# i: pointer to wavelength array
int	nwave		# i: number of wavelength points
char	form[ARB]	# i: Form of data
char	photname[ARB]	# i: Name of photometry file
char	grtbl[ARB]	# i: Graph table
char	cmptbl[ARB]	# i: Component table
real	ymin		# i: Min y value for scaling
real	ymax		# i: Max y value for scaling

int	npdat, idat, nband, ib, i, iw
real	piv, fwhm, dat, sig, peak, part
pointer	sp, pdat, psig, pstar, filt, band
char	pform[SZ_FNAME,MAXPDAT], pmode[SZ_LINE,MAXPDAT], mode[SZ_LINE]
char	mode1[SZ_LINE], mode2[SZ_LINE] 
real	bhivr()
string	nulcmd		""

begin

	# Allocate Memory
	call smark( sp )
	call salloc( pdat, nwave, TY_REAL)
	call salloc( psig, nwave, TY_REAL)
	call salloc( band, nwave, TY_REAL)
	call salloc( pstar, SZ_FNAME, TY_CHAR)

	# Plot photometry

	npdat = MAXPDAT
	call loadphot( npdat, photname, Memr[pdat], Memr[psig], pform, 
	               pmode, Memc[pstar] )
	do idat=1,npdat {
	   call strcpy( pmode[1,idat], mode, SZ_LINE)
	   call splitmode( mode, nband, mode1, mode2 )
	   if( nband == 1 )  {
	      call synphot( nwave, wave, Memr[band], nulcmd, "pivlam", mode1,
	                    grtbl, cmptbl, piv)
	      call synphot( nwave, wave, Memr[band], nulcmd, "fwhmlam", mode1,
	                    grtbl, cmptbl, fwhm)
	      dat = Memr[pdat+idat-1]
	      sig = Memr[psig+idat-1]
	      call photform( mode1, 1, nwave, wave, dat, sig, pform[1,idat],
	                                            dat, sig, form )
#	      call sphotform( 1, piv, dat, sig, pform[1,idat], 
#	                      dat, sig, form )

	      if ( !IS_INDEFR (dat) ) {
	         call gamove(gp, piv - 0.5*fwhm, dat )
	         call gadraw(gp, piv + 0.5*fwhm, dat )
	         if ( !IS_INDEFR (sig) ) {
	            call gamove(gp, piv, dat - sig )
	            call gadraw(gp, piv, dat + sig )
	         }
	         call gmark(gp, piv, dat, GM_POINT, 1)
	      }
           } else {
	      call printf ("Cannot plot 2 mode photometry data.\n")
	   }

	   # plot passbands of the photometry

	   do ib=1,nband {
	      iw = 1
	      call compband( mode1, iw, grtbl, cmptbl, nwave, wave, filt )
	      peak = bhivr( Memr[filt], nwave )
              if( peak > 0. )  {
	         do i=1,nwave {
	            part = 0.2 * Memr[filt+i-1] / peak
	            Memr[filt+i-1] = ymin * ( 1.-part ) + ymax * part
	         }
	         call gseti(gp, G_PLTYPE, GL_DOTTED)
	         call gpline( gp, Memr[wave], Memr[filt], nwave)
	         call gseti(gp, G_PLTYPE, GL_SOLID)
	      }
	      call strcpy( mode2, mode1, SZ_LINE)
	   }
	}

	call sfree( sp )
	call mfree( filt, TY_REAL)
end
