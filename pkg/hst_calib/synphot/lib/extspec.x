define ISM 1
define LMC 2

## EXTSPEC -- apply interstellar extinction to a spectrum.
#
# Applies interstellar extinction to a spectrum.
#
# Input:
#	NWAVE	= number of wavelengths
#	WAVE	= wavelength set
#	INSPEC	= spectrum
#	FORM	= form of spectrum data
#	EBMV	= extinction parameter (mags) (may be negative)
# Output:
#	OUTSPEC	= extincted spectrum
#
# March 1986 Keith Horne @ STScI
# Mar 1988 KDH @ STScI - add FORM parameter to distinguish magnitudes
# Dec 1988 E. Medeiros -- SPP version
# May 1992 D. Bazell -- Check input spectrum for INDEF
# Jul 1992 D. Bazell -- Add provision for LMC extinction

include <tbset.h>

procedure extspec ( nwave, wave, inspec, form, extlaw, ebmv, outspec )

int nwave			# number of spectrum flux samples
int strmatch			# string match function value
int i 				# loop counter

real wave[ARB]			# spectrum flux wavelength samples
real inspec[ARB]		# spectrum flux data
real ebmv	      		# the extiction in magnitudes
int  extlaw                     # Choose which extinction law to use
real outspec[ARB]		# obsfuscated spectrum fluxes
real extmag()			# extinction function for std ISM extinction
real extlmag()                  # extinction function for LMC extinction

char form[SZ_LINE]		# spectrum flux data representation
char mag_pat[SZ_LINE,2]		# magnitude representation pattern

begin

   #initialize the magnitude pattern string
   call strcpy ( "MAG", mag_pat[1,1], 3 )
   call strcpy ( "mag", mag_pat[1,2], 3 )

   # a nul spectrum flux array has been input
   if ( nwave <= 0 ) {
      call printf(" No spectrum flux data has been input\n")
      return
   }

   # spectrum flux not extinguished
   if ( ebmv == 0.0 ) {
      do i = 1, nwave 
	 outspec[i] = inspec[i]
      call printf(" Spectrum flux has not been attenuated\n")
      return
   }

   # the spectrum flux representation is in magnitudes
   if ( strmatch ( form, mag_pat[1,1] ) > 0 ||
	strmatch ( form, mag_pat[1,2] ) > 0 ) { 

      # Switch on type of extinction
      switch (extlaw) {
	 case ISM:
	    # set spectrum flux to reflect attenuation in magnitudes
	    do i =1, nwave {
	       if (!IS_INDEFR (inspec[i]))
		  outspec[i] = inspec[i] + extmag( wave[i], ebmv )
	       else
		  outspec[i] = INDEF
	    }

	 case LMC:
	    # set spectrum flux to reflect antenuation in magnitudes
	    do i =1, nwave {
	       if (!IS_INDEFR (inspec[i]))
		  outspec[i] = inspec[i] + extlmag( wave[i], ebmv )
	       else
		  outspec[i] = INDEF
	    }
      }

   # spectrum flux representation is not in magnitudes
   }else{
      switch (extlaw) {
	 case ISM:
	    # set spectrum flux to reflect antenuation in intensity
	    do i = 1, nwave {
	       if (!IS_INDEFR (inspec[i]))
		  outspec[i] = inspec[i] * 
			       10.0 ** ( -0.4 * extmag( wave[i], ebmv ))
	       else
		  outspec[i] = INDEF
	    }

	 case LMC:
	    # set spectrum flux to reflect antenuation in intensity
	    do i = 1, nwave {
	       if (!IS_INDEFR (inspec[i]))
		  outspec[i] = inspec[i] * 
			       10.0 ** ( -0.4 * extlmag( wave[i], ebmv ))
	       else
		  outspec[i] = INDEF
	    }
      }
   }
   return
   
end
