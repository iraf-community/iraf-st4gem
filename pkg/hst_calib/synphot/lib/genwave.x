## GENWAVE -- generate a wavelength set interactively
#
# 	     Generate a wavelength set to be used in synthetic photometry
#
# Input:
#	NWAVE	= Maximum number of wavelengths
# Output:
#	NWAVE	= Number of wavelengths
#	WAVE	= Wavelengths (Angstroms)
#
# Nov 1985 KDH @ STScI	
# Jan 1989 E. Medeiros -- convertion to SPP
#
include <tbset.h>
procedure genwave ( nwave, wave )

int	nwave				# maximum wavelengths requested
int	npix				# actual number of wavelengths returned
int	maxpix				# maximum wavelengths requested
int	status				# STDIN input status word
int	i				# loop counter

real	wave[ARB]			# array of generated wavelengths
real	clight				# speed of light in kilometers
real	dwave				# wavelength range in Angstroms
real	meanang				# mean number of Angstroms per pixel
real	dvel				# velocity width Km/s/pixel
real	meankm				# mean velocity
real    wmin                            # minimum wavelength
real    wmax                            # maximum wavelength
real    clgetr()

double	sum				# running count of pixels
double	add				# numerical hand-job

bool	velocity			# switch for setting the sample
                                        # spacing by line width

begin
   
   # intialize constants and limits
   clight = 2.997925e5
   npix = 1101
   velocity = false
   maxpix = nwave
   status = 0
   nwave = 0
   
   # Get the cl parameter values
   wmin = clgetr( "wmin" )
   wmax = clgetr( "wmax" )
   dwave = clgetr( "dwave" )
   dvel = clgetr( "dvelocity" )
   
   # If a velocity interval was set then set "velocity = true"
   if ( dvel > 0 )
      velocity = true
   
   # set the number of wavelengths returned to reflect minimum value
   npix = MIN( npix, maxpix )
   
   if ( velocity ) { 
      
      meanang = wmin * ( EXP( (npix-1) * dvel / clight) - 1.) / npix
      call printf (
		   "km/s/pixel = %0.3g Mean Angstroms/pixel = %0.3g\n")
      call pargr ( dvel )
      call pargr ( meanang )
      call flush ( STDOUT )
      
   }else{
      
      meankm = clight * ( ALOG( wmin+dwave * (npix-1)) - 
			  ALOG( wmin ) ) / npix
      call printf ( "A/pixel = %0.3g Mean km/s/pixel = %0.3g\n")
      call pargr ( dwave )
      call pargr ( meankm )
      call flush ( STDOUT )
   }
   if ( npix > maxpix ) {

      call printf ( "** WARNING: Too many pixels. %0d\n" )
      call pargi ( maxpix )
      call flush ( STDOUT )
   }                 
   
   if ( status == 0 ) {
      if ( velocity ) {
	 npix = 1 + NINT( clight * ALOG( wmax / wmin ) / dvel )
	 wmax = wmin * EXP( (npix-1) * dvel / clight )
      }else{
	 npix = 1 + NINT( (wmax - wmin) / dwave )
	 wmax = wmin + dwave * (npix-1)
      }
   }

   if ( status == 0 ) {
      if ( dwave> 0.0 ) {
	 velocity = false
	 npix = 1 + NINT( (wmax - wmin) / dwave )
	 wmax = wmin + dwave * (npix-1)
      }
   }

   if ( status == 0 ) {
      if ( dvel > 0.0 ) {
	 velocity = true
	 npix = 1 + NINT( clight * ALOG( wmax / wmin ) / dvel )
	 wmax = wmin * EXP( (npix-1) * dvel / clight )
      }
   }

   if ( status == 0 ) {
      if ( npix > 0 ) {
	 npix = MIN( maxpix, npix )
	 if ( velocity ) {
	    wmax = wmin * EXP( (npix-1) * dvel / clight )
	 }else{
	    wmax = wmin + (npix-1) * dwave
	 }
      }
   }

   if ( npix <= maxpix ) {
      nwave = npix
      if ( velocity ) {
	 sum = ALOG( wmin )
	 add = dvel / clight
	 do i = 1, nwave {
	    wave[i] = DEXP( sum )
	    sum = sum + add
	 }
      }else{
	 sum = wmin
	 add = dwave
	 do i = 1, nwave {
	    wave[i] = sum
	    sum = sum + add
	 }
      }
      call printf ( 
		   "Range = %0.4g %0.4g Number of Pixels = %0.4d\n" )
      call pargr ( wave[1] )
      call pargr ( wave[nwave] )
      call pargi ( nwave )
      call flush ( STDOUT )
      return
   }
end
