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
include "../plspec/plspec.h"
procedure genwave ( nwave, wave )

int	nwave				# maximum wavelengths requested
int	npix				# actual number of wavelengths returned
int	ncom				# number of interactive commands
int	maxpix				# maximum wavelengths requested
int	scan				# input line scan function value
int	status				# STDIN input status word
int	icom				# command number
int	new				# interactive pixel count adjustment
int	i				# loop counter

char	command[SZ_LINE,4]		# a command string
char	token[SZ_LINE]			# the GENWAVE command token

real	wave[ARB]			# array of generated wavelengths
real	clight				# speed of light in kilometers
real	ang1				# starting wavelength in Angstoms
real	ang2				# ending wavelength in Angstroms
real	dang				# wavelength range in Angstroms
real	meanang				# mean number of Angstroms per pixel
real	dvel				# velocity width Km/s/pixel
real	meankm				# mean velocity
real	val1				# lower limit input Angstrom range 
real	val2				# upper limit input Angstrom range
real	value				# input Angstrom range value

double	sum				# running count of pixels
double	add				# numerical hand-job

bool	velocity			# switch for setting the sample
					# spacing by line width

begin

	# intialize constants and limits
	clight = 2.997925e5
	ang1 = 1000.
	ang2 = MAXWAVE.
	dang = 10.
	npix = 1101
	velocity = false
	ncom = 4
	maxpix = nwave
	nwave = 0

	# set the number of wavelengths returned to reflect minimum value
	npix = MIN( npix, maxpix )

	# initialize the command string request
	call strcpy ("RANGE ... set wavelength range", command[1,1], SZ_LINE)
	call strcpy ("ANG   ... set wavelength interval", command[1,2],
								     SZ_LINE)
	call strcpy ("VEL   ... set velocity interval", command[1,3], SZ_LINE)
	call strcpy ("PIXEL ... set number of pixels", command[1,4], SZ_LINE)
	call strcpy ("GENWAVE", token[1], 7 )

	repeat {

	   # print the current state of the wavelength sampling array
	   # the Angstrom range and the number of pixels
	   call printf ( "Angstrom range: %0.3g %0.3g  Pixels: %10d\n" )
	   call pargr ( ang1 )                   
	   call pargr ( ang2 )
	   call pargi ( npix )
	   call flush ( STDOUT )
	   
	   if ( velocity ) { 

	      meanang = ang1 * ( EXP( (npix-1) * dvel / clight) - 1.) / npix
	      call printf (
		"km/s/pixel = %0.3g Mean Angstroms/pixel = %0.3g\n")
	      call pargr ( dvel )
	      call pargr ( meanang )
	      call flush ( STDOUT )

	   }else{

	      meankm = clight * ( ALOG( ang1+dang * (npix-1)) - 
				ALOG( ang1 ) ) / npix
	      call printf ( "A/pixel = %0.3g Mean km/s/pixel = %0.3g\n")
	      call pargr ( dang )
	      call pargr ( meankm )
	      call flush ( STDOUT )
	   }
	   if ( npix > maxpix ) {

	      call printf ( "** WARNING: Too many pixels. %0d\n" )
	      call pargi ( maxpix )
	      call flush ( STDOUT )
	   }                 

	   call selectc ( token, ncom, command, icom )

	   switch (icom) {

           case 1: # Set wavelength range

	      call printf ( "Set Angstrom range \n" )
	      call flush ( STDOUT )
	      status = scan()
	      call gargr ( val1 )
  	      call gargr ( val2 )
	      if ( status == 0 ) {
	         if (val1 > 0 && val2 >= val1 ) {
		    ang1 = val1
		    ang2 = val2
		 }
		 if ( velocity ) {
		    npix = 1 + NINT( clight * ALOG( ang2 / ang1 ) / dvel )
		    ang2 = ang1 * EXP( (npix-1) * dvel / clight )
		 }else{
		    npix = 1 + NINT( (ang2 - ang1) / dang )
		    ang2 = ang1 + dang * (npix-1)
		 }
	      }

	   case 2: # Set Angstrom interval

	      call printf ( "Set Angstroms/pixel \n" )
	      call flush ( STDOUT )
	      status = scan()
	      call gargr ( value )
	      if ( status == 0 ) {
	         if ( value > 0.0 ) {
		    velocity = false
		    dang = value
		    npix = 1 + NINT( (ang2 - ang1) / dang )
		    ang2 = ang1 + dang * (npix-1)
		 }
	      }

	   case 3: # Set veloctiy interval

	      call printf ( "Set Km/s/pixel \n" )
	      call flush ( STDOUT )
	      status = scan()
	      call gargr ( value )
	      if ( status == 0 ) {
	         if ( value > 0.0 ) {
		    velocity = true
		    dvel = value
		    npix = 1 + NINT( clight * ALOG( ang2 / ang1 ) / dvel )
		    ang2 = ang1 * EXP( (npix-1) * dvel / clight )
		 }
              }
	   case 4: # Set number of pixels

	      call printf ( "Set number of Pixels \n" )
	      call flush ( STDOUT )
	      status = scan()
	      call gargr ( new )
	      if ( status == 0 ) {
	         if ( npix > 0 ) {
		    npix = MIN( maxpix, new )
		    if ( velocity ) {
		       ang2 = ang1 * EXP( (npix-1) * dvel / clight )
		    }else{
		       ang2 = ang1 + (npix-1) * dang
		    }
                 }
              }
	   default: # Generate wavelength set

	      if ( npix <= maxpix ) {
		 nwave = npix
		 if ( velocity ) {
		    sum = ALOG( ang1 )
		    add = dvel / clight
		    do i = 1, nwave {
		       wave[i] = DEXP( sum )
		       sum = sum + add
		    }
                 }else{
		    sum = ang1
		    add = dang
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
	   }
	}
end
