include	<gset.h>

# PLSPHOT -- Plot spectrophotometry data
#
# Sep 1993 Howard Bushouse; Fixed bug in plotting ysig in C[ontinuous] mode

procedure plsphot( gp, npt, xdat, ydat, xfwhm, ysig, errtyp )

pointer	gp		# i: pointer to graphics device
int	npt		# i: number of points
real	xdat[ARB]	# i: wavelength array
real	ydat[ARB]	# i: data array
real	xfwhm[ARB]	# i: FWHM of wavelength points
real	ysig[ARB]	# i: sigmas for y data
char	errtyp[ARB]	# i: array of error flags

real	x1, x2, dx
int	ip
int	stridx()

# Oct 1989  Dave Bazell - SPP version
# Jan 1990  DB - rewrite using different error types

begin

	call strlwr( errtyp )

	# If P[oint] is specified, plot desired set of error bars (V[ertical]
	# and/or H[orizontal]) for each point

	if( stridx("p", errtyp) > 0 ) {

	   # Vertical and Horizontal error bars
	   if( stridx("v",errtyp) > 0 && stridx("h",errtyp) > 0 ) {

	      do ip = 1,npt {
	         if ( !IS_INDEFR (xfwhm[ip]) )
	            call gmark( gp, xdat[ip], ydat[ip], GM_HLINE,-xfwhm[ip],1.)
	         else
	            call gmark( gp, xdat[ip], ydat[ip], GM_POINT, 1, 1)

	         if ( !IS_INDEFR (ysig[ip]) )
	            call gmark( gp, xdat[ip], ydat[ip], GM_VLINE,1.,-ysig[ip])
	         else
	            call gmark( gp, xdat[ip], ydat[ip], GM_POINT, 1, 1)

	      }

	   # Vertical error bars
	   } else if (stridx("v", errtyp) > 0 ) {

	      do ip = 1,npt {

	         if ( !IS_INDEFR (ysig[ip]) )
	            call gmark( gp, xdat[ip], ydat[ip], GM_VLINE,1.,-ysig[ip])
	         else
	            call gmark( gp, xdat[ip], ydat[ip], GM_POINT, 1, 1)

	      }

	   # Horizontal error bars
	   } else if (stridx("h",errtyp) > 0) {

	      do ip = 1,npt {

	         if (!IS_INDEFR (xfwhm[ip]) )
	            call gmark( gp, xdat[ip], ydat[ip], GM_HLINE,-xfwhm[ip], 1.)
	         else
	            call gmark( gp, xdat[ip], ydat[ip], GM_POINT, 1, 1)


	      }

	   # no error bars, just points
	   } else
	      call gpmark( gp, xdat, ydat, npt, GM_POINT, 1., 1.)


	# Plot C[ontinuous] lines +- 1 sigma	
	} else if ( stridx("c",errtyp) > 0 ) {

	   do ip = 1,npt-1 {
	      if ( !IS_INDEFR (ysig[ip]) && !IS_INDEFR (ysig[ip+1]) )
	         call gline( gp, xdat[ip], ydat[ip] + ysig[ip],
	                         xdat[ip+1], ydat[ip+1] + ysig[ip+1] )
	   }

	   do ip = 1,npt-1 {
	      if ( !IS_INDEFR (ysig[ip]) && !IS_INDEFR (ysig[ip+1]) )
	         call gline( gp, xdat[ip], ydat[ip] - ysig[ip],
	                         xdat[ip+1], ydat[ip+1] - ysig[ip+1] )
	   }
     
	   # Horizontal error bars with +- 1 sigma lines
	   if ( stridx("h",errtyp) > 0 ) {

	      do ip = 1,npt {

	         if (!IS_INDEFR (xfwhm[ip]) )
	            call gmark( gp, xdat[ip], ydat[ip], GM_HLINE,-xfwhm[ip], 1.)
	         else
	            call gmark( gp, xdat[ip], ydat[ip], GM_POINT, 1., 1.)

	      }
	   }

	# Plot data as bins (histogram)
	} else if ( stridx( "b",errtyp) > 0 ) {

	   dx = (xdat[npt] - xdat[1])/(npt-1)
	   x1 = min( xdat[1], xdat[npt]) - dx/2.
	   x2 = max( xdat[1], xdat[npt]) + dx/2.
	   call hgline( gp, xdat, ydat, npt)

	   # Vertical error bars too
	   if (stridx("v", errtyp) > 0 ) {

	      do ip = 1,npt {

	         if ( !IS_INDEFR (ysig[ip]) )
	            call gmark( gp, xdat[ip], ydat[ip], GM_VLINE,1.,-ysig[ip])
	         else
	            call gmark( gp, xdat[ip], ydat[ip], GM_POINT, 1, 1)

	      }

	   }

	# Just plot polyline
	} else
	   call gpline( gp, xdat, ydat, npt)

end
