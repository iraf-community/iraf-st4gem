include "../plspec/plspec.h"

# EVALCHI2 -- Evaluate chisquare for a fit of spectrophotometry or
# photometry to simulated data.

procedure evalchi2( refspec, data, chisq, nfree )

char	refspec[ARB]	# i: Reference spectrum
char	data[ARB]	# i: Data string: "s= sfile p= pfile"
real	chisq		# o: chi squared
int	nfree		# o: number of data values

int	nsphot, nphot, nsum, nreject, isp, ip, ic, newrej, word, ndat
int	nwave, iw, npt
int	ctowrd(), strsearch()

pointer	sp, dat, sig, fwhm, form, mode, star, grtbl, cmptbl, inform, wavetab
pointer	colunits, spfile, pfile, spec, wave

real	band[1], result, chi2, bias, rms

double	sum

char	sphotlist[SZ_LINE,MAXSPEC], photlist[SZ_LINE,MAXPHOT]
char	pform[SZ_FNAME,MAXPHOT], pmode[SZ_LINE,MAXPHOT]
char	pstar[SZ_FNAME,MAXPHOT], junk[SZ_FNAME]

bool	status

string	nosphot		"Warning: No spectrophotometry data in %s"
string	rejectpts	"Chi2 for %s calculated using %d out of %d points"
string	toomuchphot	"Warning: Too much photometry data in %s"

# Nov 1993 Bernie Simon - Fixed allocation of wave array to avoid overflow

begin

	# Allocate memory
	call smark ( sp )
	call salloc( pfile, SZ_FNAME, TY_CHAR )
	call salloc( spfile, SZ_FNAME, TY_CHAR )

	call strlwr( data )

	# Parse the data string
	ip = 1
	for(word=1; ctowrd( data, ip, junk, SZ_FNAME ) > 0; word=word+1 ) {
	   if( word == 2 )
	      call strcpy( junk, Memc[spfile], SZ_FNAME)
	   else if (word == 4 )
	      call strcpy( junk, Memc[pfile], SZ_FNAME)
	   else if (word > 4 )
	      call printf( "Too many words in DATA string, continuing\n" )
	}

	if ( strsearch( Memc[spfile], "none" ) <= 0 ) {
	   call loadlist( Memc[spfile], nsphot, sphotlist, MAXSPEC )
	}

	if (strsearch(pfile, "none") <= 0) {
	   call loadlist( Memc[pfile], nphot, photlist, MAXPHOT )

	}

	if ( nphot <= 0 && nsphot <= 0 ) {
	   call printf("No data in %s\n")
	      call pargstr( data )
	   call sfree( sp )
	   return
	}

	call salloc( dat, MAXWAVE, TY_REAL)
	call salloc( sig, MAXWAVE, TY_REAL)
	call salloc( fwhm, MAXWAVE, TY_REAL)
	call salloc( form, SZ_FNAME, TY_CHAR)
	call salloc( mode, SZ_LINE, TY_CHAR)
	call salloc( star, SZ_FNAME, TY_CHAR)
	call salloc( grtbl, SZ_FNAME, TY_CHAR)
	call salloc( cmptbl, SZ_FNAME, TY_CHAR)
	call salloc( inform, SZ_FNAME, TY_CHAR)
	call salloc( wavetab, SZ_FNAME, TY_CHAR)

	# Get wavelength grid for synphot (Bug fix 10.93 BPS)

	if (nsphot == 0) {
	   call clgstr( "wavetab", Memc[wavetab], SZ_FNAME )
	   call getwave( Memc[wavetab], nwave, wave, colunits)
	} else {
	   nwave = MAXWAVE
	   call malloc( wave, nwave, TY_REAL )
	}

	call clgstr( "grtbl", Memc[grtbl], SZ_FNAME )
	call clgstr( "cmptbl", Memc[cmptbl], SZ_FNAME )

	nreject = 0
	nsum = 0
	sum = 0.d0

	# Loop through spectrophotometry

	do isp = 1, nsphot {
	   call loadsphot( nwave, sphotlist[1,isp], Memr[wave], Memr[dat],
	                   Memr[sig], Memr[fwhm], Memc[form], Memc[mode],
	                   Memc[star] )

	   if ( nwave <= 0 ) {
	      call printf( nosphot )
	         call pargstr( sphotlist[1,isp] )

	   } else {

	      iw = 1
	      call compspec( refspec, iw, Memc[grtbl], Memc[cmptbl], nwave, 
	                     wave, spec, Memc[inform])
	      call specform( nwave, Memr[wave], Memr[spec], Memc[inform],
	                     Memr[spec], Memc[form], status )

	      # Calculate chisq, which is normaized to nsum.  Therefore
	      # we multiply chisq by nsum.  sum accumulates the chisqs from
	      # all the spectra
	      call stats( Memr[dat], Memr[sig], Memr[spec], nwave, 
	                     Memc[form], chi2, bias, rms, npt )
	      sum = sum + chi2 * npt
	      nsum = nsum + npt
	   }
	}

	# Photometry loop
	do ip = 1, nphot {
	   ndat = MAXPHOT
	   call loadphot( ndat, photlist[1,ip], Memr[dat], Memr[sig], pform,
	                  pmode, pstar )
	   if ( ndat >= MAXPHOT ) {
	      call printf( toomuchphot)
	         call pargstr( photlist[1,ip] )
	   }

	   # Check for bad pixels
	   newrej = 0
	   do ic = 1, ndat {
	      if (Memr[sig+ic-1] <= 0 || Memr[sig+ic-1] == INDEFR ) {
	         nreject = nreject + 1
	         newrej = newrej + 1
	      } else {
	         call synphot(nwave, wave, band, refspec, pform[1,ic],
	                      pmode[1,ic], Memc[grtbl], Memc[cmptbl], result)
	         call stats( Memr[dat+ic-1], Memr[sig+ic-1], result, 1, 
	                        pform[1,ic], chi2, bias, rms, npt )
	         sum = sum + chi2
	         nsum = nsum + 1
	      }
	   }
	   # If no points have error bars then calculate statistics using
	   # equal weighting for all points
	   if ( newrej == ndat ) {
	      do ic = 1, ndat {
	         call synphot(nwave, wave, band, refspec, pform[1,ic],
	                      pmode[1,ic], Memc[grtbl], Memc[cmptbl], result)
	         call stats( Memr[dat+ic-1], Memr[sig+ic-1], result, 1, 
	                        pform[1,ic], chi2, bias, rms, npt )
	         sum = sum + chi2
	         nsum = nsum + 1
	      }

	   } else if ( newrej > 0 ) {
	      call printf( rejectpts )
	         call pargstr( photlist[1,ip] )
	         call pargi( nsum )
	         call pargi( ndat )
	   }
	}

	chisq = sum
	nfree = nsum

	call sfree ( sp )
	call mfree( wave, TY_REAL )
	call mfree( colunits, TY_CHAR )
end
