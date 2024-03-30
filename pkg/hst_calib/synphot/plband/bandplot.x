include	<tbset.h>

# BANDPLOT -- Plot real or synthetic passbands

procedure bandplot( gp, nwave, wave, command, grtbl, cmptbl, norm, ylog )

pointer	gp		# i: pointer to graphics device
int	nwave		# i: number of wavelength points
pointer	wave		# i: pointer to wavelength array
char	command[ARB]	# i: Synthetic photometry command
char	grtbl[ARB]	# i: Instrument graph table
char	cmptbl[ARB]	# i: Instrument component table
bool	norm		# i: Flag to normalize curves to  ymax
bool	ylog		# i: Logarithmic scaling?

real	high
real	bhivr()
int	iw, NMODE
pointer	band, sp, band1, band2

# Oct 1989  Dave Bazell - SPP version

begin
	# Allocate temporary memory
	call smark( sp )
	call salloc( band1, SZ_LINE, TY_CHAR)
	call salloc( band2, SZ_LINE, TY_CHAR)

	# Split mode into two if possible
	call splitmode( command, nmode, Memc[band1], Memc[band2] )

	# Evaluate first band
	iw = 1
	call compband(Memc[band1], iw, grtbl, cmptbl, nwave, wave, band)

	# Rescale plot if flag is set
	if ( norm ) {
	   high = bhivr( Memr[band], nwave )
	   call bmulkr( Memr[band], 1./high, Memr[band], nwave)
	}

	# set data less than or equal to zero to INDEF if log scaling
	if ( ylog )
	   call toindef( Memr[band], nwave )

	# Plot it
	call gpline(gp, Memr[wave], Memr[band], nwave )

	if ( nmode == 2 ) {

	   # Evaluate second band
	   iw = 1
	   call compband(Memc[band2], iw, grtbl, cmptbl, nwave, wave, band)

	   # Rescale plot if flag is set
	   if ( norm ) {
	      high = bhivr( Memr[band], nwave )
	      call bmulkr( Memr[band], 1./high, Memr[band], nwave)
	   }

	   # set data less than or equal to zero to INDEF if log scaling
	   if ( ylog )
	      call toindef( Memr[band], nwave )

	   # Plot it
	   call gpline(gp, Memr[wave], Memr[band], nwave )

	}

	call sfree(sp)
end
