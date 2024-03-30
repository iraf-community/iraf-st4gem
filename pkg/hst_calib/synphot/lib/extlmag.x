define R 3.1

## EXTLMAG -- Evaluation of Howarth's LMC extinction curve
#
# Values of LMC extinction are from Howarth's 1983 paper MNRAS, 203, 301.
#
# Extinction:
# x = 1/lambda(microns) = 10000/lambda(angstroms)
# X(x) = A(lambda)/E(B-V)
# R = A(V)/E(B-V) = 3.1
#
# Infrared
#     x <= 1.83
#     X(x) = [(1.86 - 0.48*x)*x - 0.1]*x
#
# Visible (Optical)
#     1.83 <= x <= 2.75
#     X(x) = R + 2.04*(x - 1.83) + 0.094*(x - 1.83)**2
#
# Ultra Violet
#     2.75 <= x <= 9.0
#     X(x) = R - 0.236 + 0.462*x + 0.105*x*x + 0.454/[(x - 4.557)**2 + 0.293]
#
# Input:
#	WAVE	= wavelength (Angstroms)
#	EBMV	= extinction parameter (mags) (can be negative)
#
# Output:
#	EXTLMAG	= extinction in magnitudes
#
#-- Dave Bazell, Jul 92, First Code

real procedure extlmag ( wave, ebmv )

real wave			# wavelength for Seaton's evaluation
real ebmv			# extinction parameter in Vmag
real x				# wavelength in microns truncated on 1000A

real extl

begin
   
   if ( wave <= 0. )
      return
   
   # Convert wave in angstroms to 1/microns
   x = 10000.0 / wave

   # Infrared - extend optical results linearly to 0 at 1/lam = 0
   if ( x <= 1.83)
      extl = ((1.86 - 0.48*x)*x - 0.1)*x

   else if ( x <= 2.75 ) { 
      extl = R + 2.04*(x - 1.83) + 0.094*(x-1.83)*(x-1.83)
      
   } else if ( x <= 10.96 ) { # continue out to lambda = 912 A
      extl = R - 0.236 + 0.462*x + 0.105*x*x +
	     0.454/( (x-4.557)*(x-4.557) + 0.293 )
   }else if ( x > 10.96 ) {
      x = 10.96
      extl = R - 0.236 + 0.462*x + 0.105*x*x +
	     0.454/( (x-4.557)*(x-4.557) + 0.293 )
      
   }

   # Rescale the magnitude extinction to the Vmag extinction
   return (extl * ebmv)

end
