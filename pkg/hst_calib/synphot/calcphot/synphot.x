include		"../lib/mac.h"

# SYNPHOT -- synthetic photometer

# jul 1988 keith horne @ stsci
# may 1989 kdh @ stsci - improve comments at top, add reset
# jul 1989 kdh @ stsci - change photpix to counts, eliminate type
# jul 1989 kdh @ stsci - use splitmode
# sep 1989 Dave Bazell - SPP Version
# Dec 1989 DB - Add ability to handle sythetic passbands
# Aug 1993 H.A.Bushouse - Fixed bugs for case nband=2 and form=obmag

procedure synphot( nwave, wave, band, command, form, mode, 
	           grtbl, cmptbl, result)

int 	nwave		# i: Number of wavelengths
pointer wave		# i: Pointer to wavelength array
real	band[ARB]	# i: Synthetic passband if mode='synthetic'
char	command[ARB]	# i: Command string used in compspec or compband
char	form[ARB]	# i: form of output spectrum
char	mode[ARB]	# i: Instrument mode string
char	grtbl[ARB]	# i: Instrument graph table
char	cmptbl[ARB]	# i: Instrument component table
real	result		# o: Returned result

char	mode1[SZ_LINE], mode2[SZ_LINE]
pointer	sp, errmsg, inform, filt, filt2, spec
int	nband, junk, iw
real	result2
bool	alldone, status

bool	streq() 
int	strsearch(), strlen(), is_magunit(), fillnull()
real	pivlam(), rmslam(), fwhmlam(), barlam(), avglam(), efflam(), rate()
real	effstim()

string	nomode		"Mode string is empty."
string	divbyzero	"Divide by zero error; Form = %s"
string	logzero		"Logarithm of zero or negative number; Form = %s"

begin

	if ( strlen(mode) == 0 )
	   call error( 1, nomode )

	# Allocate memory
	call smark( sp )
	call salloc( errmsg, SZ_LINE, TY_CHAR)
	call salloc( inform, SIZE_STR, TY_CHAR)
	alldone = false

	# Check if mode is photometric or synthetic

	if ( strsearch(mode, "synthetic") > 0 || 
	     strsearch(mode, "SYNTHETIC") > 0 ) {

	   # Copy synthetic passband to working array
	   call malloc(filt, nwave, TY_REAL)
	   call amovr(band, Memr[filt], nwave)

	} else {

	   # Evaluate photometric passband

	   call splitmode( mode, nband, mode1, mode2 )

	   # mode1 starts at index iw
	   if( nband == 1) {
	      iw = 1
	      call compband( mode1, iw, grtbl, cmptbl, nwave, wave, filt)

	   # evaluate 2 passbands
	   } else {
	      iw = 1
	      call compband( mode1, iw, grtbl, cmptbl, nwave, wave, filt)
	      iw = 1
	      call compband( mode2, iw, grtbl, cmptbl, nwave, wave, filt2)
	   }

	}

	# pivot wavelength
	if( streq(form, "pivlam") || streq(form, "PIVLAM")) {
	   result = pivlam( nwave, Memr[wave], Memr[filt] )
	   alldone = true

	# rms bandwidth
	} else if( streq(form,"rmslam") || streq(form, "RMSLAM") ) {
	   result = rmslam( nwave, Memr[wave], Memr[filt] )
	   alldone = true

	# fwhm bandwidth
	} else if(streq( form, "fwhmlam") || streq(form, "FWHMLAM") ) {
	   result = fwhmlam( nwave, Memr[wave], Memr[filt] )
	   alldone = true

	# bar wavelength
	} else if(streq( form, "barlam") || streq(form, "BARLAM") ) {
	   result = barlam( nwave, Memr[wave], Memr[filt] )
	   alldone = true

	# avg wavelength
	} else if(streq( form, "avglam") || streq(form, "AVGLAM") ) { 
	   result = avglam( nwave, Memr[wave], Memr[filt] )
	   alldone = true
	}

	# If function performed then free memory and return
	if (alldone) {
	   call sfree(sp)
	   call mfree( filt, TY_REAL)
	   if (nband > 1) {
	      call mfree( filt2, TY_REAL)
	   }
	   return
	}
	   
	# evaluate the spectrum, 'spec' is allocated in comspec
	iw = 1
	call compspec( command, iw, grtbl, cmptbl, nwave, wave, 
	               spec, Memc[inform])

	# effective wavelength
	if(streq( form,"efflam") || streq(form,"EFFLAM") ) {
	   result = efflam( nwave, Memr[wave], Memr[filt], Memr[spec],
	                    Memc[inform] )
	   alldone = true

	# count rate
	} else if(streq( form, "counts") || streq(form, "COUNTS") ||
	          streq(form,"obmag") || streq(form,"OBMAG") ) {
	   result = rate( nwave, Memr[wave], Memr[filt], Memr[spec],
	                  Memc[inform] )

	   if ( streq(form,"obmag") || streq(form,"OBMAG") ) {
		if (result > 0.) {
		    result = -2.5 * alog10(result)
		} else {
	            call printf(" result: %e; mode1: %s\n")
	               call pargr(result)
	               call pargstr(mode1)
	            call sprintf( Memc[errmsg], SZ_LINE, logzero)
	               call pargstr(form)
	            call error(1, Memc[errmsg])
		}
	   }

	   # count rate for second passband
	   if( nband > 1 ) {
	      result2 = rate( nwave, Memr[wave], Memr[filt2], Memr[spec],
	                      Memc[inform] )

	      # count rate ratio for COUNTS
	      if ( streq(form,"counts") || streq(form,"COUNTS" ) ) {
	         if( result2 > 0 ) {
	            result = result / result2
	         # trap divide by zero
	         } else {
	            call printf(" result: %e; mode1: %s\n")
	               call pargr(result)
	               call pargstr(mode1)
	            call printf(" result2: %e; mode2: %s\n")
	               call pargr(result2)
	               call pargstr(mode2)
	            call sprintf( Memc[errmsg], SZ_LINE, divbyzero)
	               call pargstr(form)
	            call error(1, Memc[errmsg])
	         }

	     # count rate difference for obmag
	      } else if ( streq(form,"obmag") || streq(form,"OBMAG") ) {
	         if (result2 > 0) {
	            result = result + 2.5 *alog10(result2)
	         } else {
	            call printf(" result: %e; mode1: %s\n")
	               call pargr(result)
	               call pargstr(mode1)
	            call printf(" result2: %e; mode2: %s\n")
	               call pargr(result2)
	               call pargstr(mode2)
	            call sprintf( Memc[errmsg], SZ_LINE, logzero)
	               call pargstr(form)
	            call error(1, Memc[errmsg])
	         }
	      }
	   }
	   alldone = true
	}

	# If function performed then free memory and return
	if (alldone) {
	   call sfree(sp)
	   call mfree( filt, TY_REAL)
	   if (nband > 1)
	      call mfree( filt2, TY_REAL)

	   return
	}

	# transform the spectrum
	call specform( nwave, Memr[wave], Memr[spec], Memc[inform], 
	               Memr[spec], form, status )

	if (is_magunit (form) == NO) {
	    junk = fillnull (0.0, nwave, Memr[wave], Memr[spec])
	} else 	{
	    junk = fillnull (100.0, nwave, Memr[wave], Memr[spec])
	}

	# calculate effective stimulus
	result = effstim( nwave, Memr[wave], Memr[filt], Memr[spec], form )

	# calculate effective stimulus for second passband
	if( nband > 1 ) {
	   result2 = effstim( nwave, Memr[wave], Memr[filt2], Memr[spec], form )

	   # color index
	   if( strsearch(form, "mag" ) > 0 || strsearch(form,"MAG") > 0 ) 
	      result = result - result2

	   # flux ratio
	   else if( result2 > 0 ) {
	      result = result / result2

	      # trap divide by zero
	   } else {
	      call printf(" result: %e; mode1: %s\n")
	         call pargr(result)
	         call pargstr(mode1)
	      call printf(" result2: %e; mode2: %s\n")
	         call pargr(result2)
	         call pargstr(mode2)
	      call sprintf( Memc[errmsg], SZ_LINE, divbyzero)
	         call pargstr(form)
	      call error(1, Memc[errmsg])
	   }
	}

	call sfree(sp)
	call mfree( filt, TY_REAL)
	if (nband > 1)
	   call mfree( filt2, TY_REAL)

end
