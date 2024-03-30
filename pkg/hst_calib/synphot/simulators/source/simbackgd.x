include	<imhdr.h>
include	"simtwo.h"

#* HISTORY *
#* B.Simon	21-Aug-95	extracted from simimg

# SIMBACKGD -- Compute simulated background image

procedure simbackgd ()

#--
pointer	obsmode		# telescope observation mode
pointer	output		# output image name
real	exptime		# exposure time
double	det_ra		# telescope right ascension
double	det_dec		# telescope declination
double	det_ang		# telescope position angle
bool	verbose		# print output message?
pointer	wavtable	# wavelength array table
real	eshine		# fraction of maximum earthlight
pointer	time		# time of observation
long	seed		# random number seed
pointer	detcat		# catalog of detector information
pointer	zodtab		# table of zodiacal flux
pointer	earthtab	# earthlight spectrum (flux / sec^2)
pointer	thermtab	# thermal background (flux / sec^2)
pointer	grftable	# instrument graph table
pointer	cmptable	# component name table
real	hstarea		# total telescope area
real	counts		# background counts

double	apscale, jd
int	apx, apy, nwave, degree
pointer	sp, mode, pcode, ocode, zcode 
pointer	wave, thruput, out, im
real	icounts, ocounts

string	raform	 RA_UNITS
string	notband  "Observation mode is not a bandpass"

bool	clgetb()
double	clgetd(), daytime()
int	clgeti()
pointer	immap(), imps2r()
real	clgetr()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (obsmode, SZ_LINE, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (wavtable, SZ_FNAME, TY_CHAR) 
	call salloc (time, SZ_FNAME, TY_CHAR)

	call salloc (detcat, SZ_FNAME, TY_CHAR) 
	call salloc (zodtab, SZ_FNAME, TY_CHAR)
	call salloc (earthtab, SZ_FNAME, TY_CHAR)
	call salloc (thermtab, SZ_FNAME, TY_CHAR)
	call salloc (grftable, SZ_FNAME, TY_CHAR)
	call salloc (cmptable, SZ_FNAME, TY_CHAR)

	call salloc (mode, SZ_LINE, TY_CHAR)
	call salloc (pcode, 2*SZ_LINE, TY_INT)
	ocode = pcode
	zcode = pcode + SZ_LINE

	# Read task parameters from main parameter file

	call clgstr ("obsmode", Memc[obsmode], SZ_LINE)
	call clgnone ("output", Memc[output], SZ_FNAME)

	exptime = clgetr ("exptime")
	det_ra = clgetd ("det_ra")
	det_dec = clgetd ("det_dec")
	det_ang = clgetd ("det_ang")
	verbose = clgetb ("verbose")

	call clgnone ("wavetab", Memc[wavtable], SZ_FNAME)

	# Read background parameter file

	eshine = clgetr ("eshine")
	call clgstr ("time", Memc[time], SZ_FNAME)
	seed = clgeti ("seed")

	# Read catalog and file parameter file

	call clgstr ("detcat", Memc[detcat], SZ_FNAME)
	call clgnone ("zodtab", Memc[zodtab], SZ_FNAME)
	call clgnone ("earthtab", Memc[earthtab], SZ_FNAME)
	call clgnone ("thermtab", Memc[thermtab], SZ_FNAME)

	# Read reference data parameter file

	call clgstr ("grtbl", Memc[grftable], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptable], SZ_FNAME)
	hstarea = clgetr( "area" )

	# Process the task parameters

	call put_hstarea (hstarea)

	call getnaked (Memc[obsmode], Memc[mode], SZ_LINE)
	call angtodeg (raform, det_ra)
	jd = daytime (Memc[time])

	# Initialize table cache and synphot variables

	call inisyntab
	call undefsynvar

	# Compile the observation mode and magnitude passband expressions

	call setmodewarn (NO)	# Turn off msg for unused keywords

	call expcompile (Memc[mode], Memi[ocode], SZ_LINE)
	call expcompile (ZOD_BAND, Memi[zcode], SZ_LINE)

	# Compute the wavelength and thruput arrays

	call getwavelen (Memc[wavtable], Memc[grftable], Memc[cmptable],
			 Memi[pcode], 2, SZ_LINE, wave, nwave)

	call salloc (thruput, nwave, TY_REAL)
	call syncalc (Memi[ocode], SZ_LINE, NULL, nwave, Memr[wave],
		      Memc[grftable], Memc[cmptable], Memr[thruput], degree)

	call setmodewarn (YES)
	if (degree != 0)
	    call printerr_str (notband, Memc[mode])

	# Retrieve the detector size and scale

	call getscale (Memc[detcat], Memc[mode], apscale, apx, apy)

	# Compute background counts

	call outbackgd (Memc[zodtab], det_ra, det_dec, jd, Memc[earthtab], 
			eshine, exptime, apscale, Memc[grftable], 
			Memc[cmptable], nwave, Memr[wave], Memr[thruput], 
			ocounts)

	call inbackgd (Memc[thermtab], exptime, apscale, nwave, 
		       Memr[wave], Memr[thruput], icounts)

	counts = icounts + ocounts

	# Open / create output image

	if (Memc[output] != EOS) {
	    im = immap (Memc[output], NEW_FILE, NULL)

	    IM_NDIM(im) = 2
	    IM_LEN(im,1) = apx
	    IM_LEN(im,2) = apy
	    IM_PIXTYPE(im) = TY_REAL

	    # Get outout image buffer

	    out = imps2r (im, 1, apx, 1, apy)

	    # Set output to number of counts

	    call amovkr (counts, Memr[out], apx*apy)
	    call imunmap (im)
	}

	# Write results

	if (verbose) {
	    call printf ("Background counts: %10.3g\n\n")
	    call pargr (counts)
	}

	call clputr ("counts", counts)
	# Close files and free memory

	call clssyntab

	call mfree (wave, TY_REAL)
	call sfree (sp)
end
