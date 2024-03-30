include	<imhdr.h>
include	"simtwo.h"

#* HISTORY *
#* B.Simon	21-Aug-95	extracted from simimg

# SIMNOISE -- Calculate simulated noise model

procedure simnoise ()

#--
pointer	obsmode		# telescope observation mode
pointer	output		# output image name
pointer	input		# input image name
pointer	noise		# noise expression
real	exptime		# exposure time
int	nread		# number of readouts
pointer	wavtable	# wavelength array table
long	seed		# random number seed
pointer	detcat		# catalog of detector information
pointer	grftable	# instrument graph table
pointer	cmptable	# component name table
real	hstarea		# total telescope area

bool	quant
double	apscale
int	apx, apy, nwave, degree
pointer	sp, mode, pcode, ncode, wave, thruput
pointer	im, imi, out, in

string	notband  "Observation mode is not a bandpass"
string	badsize  "Input image size does not match detector size"

int	clgeti()
pointer	immap(), imps2r(), imgs2r()
real	clgetr()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (obsmode, SZ_LINE, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (input, SZ_FNAME, TY_CHAR)
	call salloc (noise, SZ_LINE, TY_CHAR)
	call salloc (wavtable, SZ_FNAME, TY_CHAR) 
	call salloc (detcat, SZ_FNAME, TY_CHAR) 
	call salloc (grftable, SZ_FNAME, TY_CHAR)
	call salloc (cmptable, SZ_FNAME, TY_CHAR)

	call salloc (mode, SZ_LINE, TY_CHAR)
	call salloc (pcode, SZ_LINE, TY_INT)
	call salloc (ncode, SZ_LINE, TY_INT)

	# Read task parameters from main parameter file

	call clgstr ("obsmode", Memc[obsmode], SZ_LINE)
	call clgstr ("output", Memc[output], SZ_FNAME)
	call clgnone ("input", Memc[input], SZ_FNAME)
	call clgnone ("noise", Memc[noise], SZ_FNAME)

	exptime = clgetr ("exptime")
	nread = clgeti ("nread")
	quant = false

	call clgnone ("wavetab", Memc[wavtable], SZ_FNAME)

	# Read catalog and file parameter file

	call clgstr ("detcat", Memc[detcat], SZ_FNAME)

	# Read reference data parameter file

	call clgstr ("grtbl", Memc[grftable], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptable], SZ_FNAME)
	hstarea = clgetr( "area" )

	# Process the task parameters

	call put_hstarea (hstarea)
	call getnaked (Memc[obsmode], Memc[mode], SZ_LINE)

	# Initialize table cache and synphot variables

	call inisyntab
	call undefsynvar

	# Compile the observation mode and magnitude passband expressions

	call setmodewarn (NO)	# Turn off msg for unused keywords

	call expcompile (Memc[mode], Memi[pcode], SZ_LINE)

	# Compile the noise expression

	if (Memc[noise] == EOS) {
	    call rdnoise (Memc[mode], Memc[grftable], Memc[cmptable], 
			  exptime, nread, Memi[ncode], SZ_LINE)

	} else {
	    call compnoise (Memc[noise], exptime, nread, Memi[ncode], SZ_LINE)
	}

	# Compute the wavelength and thruput arrays

	call getwavelen (Memc[wavtable], Memc[grftable], Memc[cmptable],
			 Memi[pcode], 1, SZ_LINE, wave, nwave)

	call salloc (thruput, nwave, TY_REAL)
	call syncalc (Memi[pcode], SZ_LINE, NULL, nwave, Memr[wave],
		      Memc[grftable], Memc[cmptable], Memr[thruput], degree)

	call setmodewarn (YES)
	if (degree != 0)
	    call printerr_str (notband, Memc[mode])

	# Retrieve the detector size and scale

	call getscale (Memc[detcat], Memc[mode], apscale, apx, apy)

	# Open / create output image

	if (Memc[input] == NULL) {
	    # Open new output image, set flux to zero

	    im = immap (Memc[output], NEW_FILE, NULL)

	    IM_NDIM(im) = 2
	    IM_LEN(im,1) = apx
	    IM_LEN(im,2) = apy
	    IM_PIXTYPE(im) = TY_REAL

	    out = imps2r (im, 1, apx, 1, apy)
	    call aclrr (Memr[out], apx*apy)

	} else {
	    # Open output image as copy of input image, set flux to input

	    imi = immap (Memc[input], READ_ONLY, NULL)
	    if (IM_LEN(imi,1) != apx || IM_LEN(imi,2) != apy)
		call printerr_str (badsize, Memc[input])

	    im = immap (Memc[output], NEW_COPY, imi)

	    in = imgs2r (imi, 1, apx, 1, apy)
	    out = imps2r (im, 1, apx, 1, apy)
	    call amovr (Memr[in], Memr[out], apx*apy)

	    call imunmap (imi)
	}


	# Add calculated noise to output image

	call addnoise (Memr[out], apx, apy, quant, seed, Memi[ncode])

	# Close files and free memory

	call imunmap (im)
	call clssyntab

	call mfree (wave, TY_REAL)
	call sfree (sp)
end
