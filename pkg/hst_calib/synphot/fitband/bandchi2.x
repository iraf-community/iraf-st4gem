include "../plspec/plspec.h"

# BANDCHI2 -- Evaluate chisquare for a fit of spectrophotometry or
# photometry to simulated data.

procedure bandchi2( bandmodel, data, chisq, bias, nfree )

char	bandmodel[ARB]	# i: Passband model to be fit
char	data[ARB]	# i: Data string: "s= sfile p= pfile"
real	chisq		# o: chi squared
real	bias[ARB]	# o: Obs - Predicted for each data point
int	nfree		# o: number of data values

int	nsphot, nphot, nsum, nreject, isp, ip, ic, newrej, ndat
int	nwave, npt, nchar
int	nowhite()

bool	streq(), strne()

pointer	sp, dat, sig, fwhm, form, mode, star, grtbl, cmptbl, wavetab
pointer	colunits, wave, yform, test, passband

real	band[1], result, chi2, rms, piv

double	sum

char	sphotlist[SZ_LINE,MAXSPEC], photlist[SZ_LINE,MAXPHOT]
char	pform[SZ_FNAME,MAXPHOT], pmode[SZ_LINE,MAXPHOT]
char	pstar[SZ_FNAME,MAXPHOT], oldmode[SZ_LINE]

# Include common block containing spectrum files corresponding to target ids
include	"../lib/targets.h"

string	nodata		"No data to fit in BANDCHI2"
string	rejectpts	"Chi2 for %s calculated using %d out of %d points\n"
string	toomuchphot	"Warning: Too much photometry data in %s\n"
string	nulcmd		""

begin

	call parsedat( data, sphotlist, nsphot, MAXSPEC, 
		       photlist, nphot, MAXPHOT )

	# Allocate memory
	call smark ( sp )
	call salloc( dat, MAXWAVE, TY_REAL)
	call salloc( sig, MAXWAVE, TY_REAL)
	call salloc( fwhm, MAXWAVE, TY_REAL)
	call salloc( form, SZ_FNAME, TY_CHAR)
	call salloc( mode, SZ_LINE, TY_CHAR)
	call salloc( star, SZ_FNAME, TY_CHAR)
	call salloc( yform, SZ_FNAME, TY_CHAR)
	call salloc( test, SZ_LINE, TY_CHAR )
	call salloc( grtbl, SZ_FNAME, TY_CHAR)
	call salloc( cmptbl, SZ_FNAME, TY_CHAR)
	call salloc( wavetab, SZ_FNAME, TY_CHAR)
	call salloc( passband, SZ_LINE, TY_CHAR)

	# Get wavelength grid for synphot
	call clgstr( "wavetab", Memc[wavetab], SZ_FNAME )
	call getwave( Memc[wavetab], nwave, wave, colunits)
	call clgstr( "yform", Memc[yform], SZ_FNAME )
	call clgstr( "grtbl", Memc[grtbl], SZ_FNAME)
	call clgstr( "cmptbl", Memc[cmptbl], SZ_FNAME )

	nreject = 0
	nsum = 0
	sum = 0.d0

	# Photometry loop
	do ip = 1, nphot {
	   ndat = MAXPHOT
	   call loadphot( ndat, photlist[1,ip], Memr[dat], Memr[sig], pform,
	                  pmode, pstar )

	   if ( ndat >= MAXPHOT ) {
	      call printf( toomuchphot)
	         call pargstr( photlist[1,ip] )
	   }

	   # Change all data to consistent form 
	   do ic = 1,ndat {

	      # Calculate pivot wavelength only if mode has changed.  Strip
	      # blanks from mode strings for comparison
	      nchar = nowhite( pmode[1,ic],Memc[test],SZ_LINE )
	      nchar = nowhite( oldmode, oldmode, SZ_LINE )

	      if ( strne( Memc[test], oldmode ) ) {
	         call synphot( nwave, wave, band, nulcmd, "pivlam", 
	                       pmode[1,ic], Memc[grtbl], Memc[cmptbl], piv)
	         call strcpy( pmode[1,ic], oldmode, SZ_LINE )
	      }

	      call photform( pmode[1,ic], 1, nwave, wave, 
	                     Memr[dat+ic-1], Memr[sig+ic-1], pform[1,ic],
	                     Memr[dat+ic-1], Memr[sig+ic-1], Memc[yform] )

	      call strcpy ( Memc[yform], pform[1,ic], SZ_FNAME )
	      call strlwr( pform[1,ic] )
	   }

	   # Check for bad pixels
	   newrej = 0
	   do ic = 1, ndat {

	      if (Memr[sig+ic-1] <= 0 || Memr[sig+ic-1] == INDEFR ) {
	         nreject = nreject + 1
	         newrej = newrej + 1


	      } else {

	         # targetid and specfile are declared in the include file
	         # targets.h.  findsphot makes the correspondence between
	         # pstar[1,ic], which is the targetid of the current phot.obs,
	         # and targetid, and returns the value isp which points
	         # to the isp-th spectrum in specfile

	         # If the targetid = "none" then assume that the id in pstar[]
	         # is the file name
	         if ( streq( targetid[1,1], "none" ) ) {
	            call findsphot( pstar[1,ic], pstar, ndat, isp )
	            call strcpy( pstar[1,isp], specfile[1,isp], SZ_FNAME )
	         } else
	            call findsphot( pstar[1,ic], targetid, ntarget, isp )

	         # Insert the OBSMODE of the current observation if
	         # BANDMODEL has the variable 'obsmode' in it
	         call insertmode( pmode[1,ic], bandmodel, Memc[passband] )
	         call synphot(nwave, wave, band, specfile[1,isp], pform[1,ic],
	                      Memc[passband],Memc[grtbl],Memc[cmptbl],result)

	         call stats( Memr[dat+ic-1], Memr[sig+ic-1], result, 1, 
	                        pform[1,ic], chi2, bias[nsum+1], rms, npt )
	         sum = sum + chi2
	         nsum = nsum + 1
	      }
	   }

	   # If no points have error bars then calculate statistics using
	   # equal weighting for all points
	   if ( newrej == ndat ) {

	      do ic = 1, ndat {

	         # targetid and specfile are declared in the include file
	         # targets.h.  findsphot makes the correspondence between
	         # pstar[1,ic], which is the targetid of the current phot.obs,
	         # and targetid, and returns the value isp which points
	         # to the isp-th spectrum in specfile

	         # If the targetid = "none" then assume that the id in pstar[]
	         # is the file name
	         if ( streq( targetid[1,1], "none" ) ) {
	            call findsphot( pstar[1,ic], pstar, ndat, isp )
	            call strcpy( pstar[1,isp], specfile[1,isp], SZ_FNAME )
	         } else
	            call findsphot( pstar[1,ic], targetid, ntarget, isp )

	         # Insert the OBSMODE of the current observation if
	         # BANDMODEL has the variable 'obsmode' in it
	         call insertmode( pmode[1,ic], bandmodel, Memc[passband] )
	         call synphot(nwave, wave, band, specfile[1,isp], pform[1,ic],
	                      Memc[passband],Memc[grtbl],Memc[cmptbl],result)

	         call stats( Memr[dat+ic-1], Memr[sig+ic-1], result, 1, 
	                        pform[1,ic], chi2, bias[nsum+1], rms, npt )
	         sum = sum + chi2
	         nsum = nsum + 1
	      }

	   # Only some points had errors so report
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
