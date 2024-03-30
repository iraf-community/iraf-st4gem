# BOX -- Create a rectangular passband

procedure box( script, nwave, wave, iw, band)

char	script[ARB]	# i: Command script
int	nwave		# i: Number of wavelengths in wave
real	wave[ARB]	# i: Wavelength array
int	iw		# io: Index of position in script
real	band[ARB]	# o: Output passband

int	ib, nchar
int	ctor()
pointer	errmsg
real	center, fwhm, wmin, wmax
bool	setone

string	negwidth	"BOX:  FWHM must be greater than zero, %r"

begin

	# Get center wavelength

	nchar = ctor(script, iw, center)

	# Get FWHM

	nchar = ctor(script, iw, fwhm)

	if ( fwhm <= 0. ) {
	   call malloc( errmsg, SZ_LINE, TY_CHAR )
	   call sprintf( Memc[errmsg], SZ_LINE, negwidth )
	      call pargr( fwhm )
	   call error( 1, Memc[errmsg] )
	}

	wmin = center - fwhm/2.
	wmax = center + fwhm/2.

	setone = false

	do ib = 1, nwave {

	   if ( wave[ib] < wmin )
	      band[ib] = 0.
	   else if ( wave[ib] >= wmax && setone )
	      band[ib] = 0.
	   else {
	      setone = true
	      band[ib] = 1.
	   }

	}

end
	
