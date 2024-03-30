include	<gset.h>
include	"plspec.h"

# SPECPLOT -- Plot spectra and synthetic photometry

procedure specplot( gp, nwave, wave, form, modelist, nmode, ebmv, specname, 
	            grtbl, cmptbl, hstarea)

pointer	gp			# i: Graphics pointer
int	nwave			# i: Number of wavelength points in wave
pointer	wave			# i: Pointer to wavelength array
char	form[ARB]		# i: Form of spectral data
char	modelist[SZ_LINE, ARB]	# i: List of instrument modes
int	nmode			# i: Number of modes
real	ebmv			# i: Value of extinction for synth. Phot.
char	specname[ARB]		# i: Current spectrum name
char	grtbl[ARB]		# i: Graph table 
char	cmptbl[ARB]		# i: Component table
real	hstarea			# i: HST area

int	nband, imode, ib, iw
real	dat, piv, fwhm
char	mode[SZ_LINE], mode1[SZ_LINE], mode2[SZ_LINE], command[SZ_LINE]
pointer	sp, inform, spec, histi, histr, band
bool	status

bool	is_simple()
int	strsearch()

string	bandfmt	    "(%s)*band(%s)"
string	expfmt	    "(%s)*(%s)"
string	nulcmd	    ""

# Bug Fixes
# Jan 1990 DB; Insert 'ebmv' command before '*'
# Nov 1993 Bernie Simon; build ebmv string with applyebmv()

begin

	# Allocate memory
	call smark( sp )
	call salloc( inform, SZ_FNAME, TY_CHAR )
	call salloc( histi, nwave, TY_INT)
	call salloc( histr, nwave, TY_REAL)
	call salloc( band, nwave, TY_REAL)

	# Plot detected counts ( spectrum times passband ) as histogram
	if( (strsearch(form,"counts") > 0 || strsearch(form,"photpix") > 0 ||
	     strsearch(form,"obmag") > 0 ) &&  nmode > 0 )  {
	   do imode=1,nmode {
	      call strcpy( modelist[1,imode], mode, SZ_LINE)
	      nband = 2
	      call splitmode( mode, nband, mode1, mode2 )
	      do ib=1,nband {
		 if (is_simple (mode1)) {
		     call sprintf(command, SZ_LINE, bandfmt)
		 } else {
		     call sprintf(command, SZ_LINE, expfmt)
		 }
		 call pargstr(specname)
		 call pargstr(mode1)

		 call applyebmv (command, command, ebmv, SZ_LINE)

	         iw = 1
	         call compspec( command, iw, grtbl, cmptbl,
	                        nwave, wave, spec, Memc[inform])
	         call specform( nwave, Memr[wave], Memr[spec], 
	                        Memc[inform], Memr[spec], form, status )

	         call hgline( gp, Memr[wave], Memr[spec], nwave )
	         call strcpy( mode2, mode1, SZ_LINE)
	      }
	   }

	# Not counts so plot the spectrum as curve not histogram
	} else {
	   call applyebmv (specname, command, ebmv, SZ_LINE)

	   iw = 1
	   call compspec( command, iw, grtbl, cmptbl,
	                  nwave, wave, spec, Memc[inform])

	   call specform( nwave, Memr[wave], Memr[spec], Memc[inform],
	                  Memr[spec], form, status )

	   call gpline( gp, Memr[wave], Memr[spec], nwave)

	   # Plot synthetic photometry
 	   do imode=1,nmode {
	      call strcpy(modelist[1,imode], mode, SZ_LINE)
	      nband = 2
	      call splitmode( mode, nband, mode1, mode2 )
	      do ib=1,nband {
	         call synphot( nwave, wave, Memr[band], nulcmd, "pivlam", mode1,
	                       grtbl, cmptbl, piv)
	         call synphot( nwave, wave, Memr[band], nulcmd, "fwhmlam", mode1,
	                       grtbl, cmptbl, fwhm)
		 call applyebmv (specname, command, ebmv, SZ_LINE)

	         call synphot( nwave, wave, Memr[band], command, form,
	                       mode1, grtbl, cmptbl, dat)

	         call gamove(gp, piv - 0.5*fwhm, dat )
	         call gadraw(gp, piv + 0.5*fwhm, dat )
	         call gmark( gp, piv, dat, GM_CIRCLE, 2., 2.)
	         call strcpy( mode2, mode1, SZ_LINE)
	      }
	   }	# next mode
	}

	call sfree( sp )
	call mfree( spec, TY_REAL)
end
