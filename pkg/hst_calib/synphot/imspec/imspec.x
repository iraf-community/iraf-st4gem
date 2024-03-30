#* HISTORY *
#* B.Simon	03-Mar-1993	original
#* B.Simon	08-Apr-1993	revised to call getftype
#* B.Simon	07-Jan-1994	added input form
#* Phil Hodge,	19-Sep-2002	Pass the complete input file name to getftype,
#				and delete the call to imgcluster.

include	"imspec.h"

# IMSPEC -- Convert an image to a synphot spectrum or vice versa

procedure imspec

#--
pointer	input		# list of input files
pointer	output		# list of output files
pointer	wave		# optional wavelength images
pointer	inform		# units of input spectrum
pointer	outform		# units of output spectrum
int	olength		# length of output file
real	badpix		# value to substitute for INDEF in output image
real	hstarea		# telescope area in cm^2

int	junk, ilen
pointer	sp, ifile, ofile, wfile, errmsg
pointer	inlist, outlist, wavlist, iwave, ispec

string	badoutlen  "Number of input files does not match output files"
string	badwavlen  "Number of input files does not match wavelength files"
string	badinput   "Cannot access input file (%s)"
string	badtype    "Input file is not image or table (%s)"

bool	isblank()
int	clgeti(), imtlen(), imtgetim(), getftype()
pointer	imtopen()
real	clgetr()

begin
	# Allocate dynamic meory for strings

	call smark (sp)
	call salloc (input, SZ_FNAME, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (wave, SZ_FNAME, TY_CHAR)
	call salloc (ifile, SZ_FNAME, TY_CHAR)
	call salloc (ofile, SZ_FNAME, TY_CHAR)
	call salloc (wfile, SZ_FNAME, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)
	call salloc (inform, SZ_FNAME, TY_CHAR)
	call salloc (outform, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgstr ("input", Memc[input], SZ_FNAME)
	call clgstr ("output", Memc[output], SZ_FNAME)
	call clgstr ("wave", Memc[wave], SZ_FNAME)
	call clgstr ("inform", Memc[inform], SZ_FNAME)
	call clgstr ("outform", Memc[outform], SZ_FNAME)
	olength = clgeti ("olength")
	badpix = clgetr ("badpix")
	hstarea = clgetr ("area") 

	call inisyntab
	call put_hstarea (hstarea)

	# Open file templates

	inlist = imtopen (Memc[input])
	outlist = imtopen (Memc[output])

	if (imtlen (inlist) != imtlen (outlist))
	    call error (1, badoutlen)

	if (isblank (Memc[wave])) {
	    wavlist = NULL
	} else {
	    wavlist = imtopen (Memc[wave])

	    if (imtlen (inlist) != imtlen (wavlist))
		call error (1, badwavlen)
	}

	# Convert each input file

	while (imtgetim (inlist, Memc[ifile], SZ_FNAME) != EOF &&
	       imtgetim (outlist, Memc[ofile], SZ_FNAME) != EOF  ) {

	    # Read wavelength file if specified

	    if (wavlist == NULL) {
		Memc[wfile] = EOS
	    } else {
		junk = imtgetim (wavlist, Memc[wfile], SZ_FNAME)
	    }

	    # Check to see if input file is image or table
	    # Convert it to its opposite

	    switch (getftype (Memc[ifile])) {
	    case NONE:
		call sprintf (Memc[errmsg], SZ_LINE, badinput)
		call pargstr (Memc[ifile])
		call error (1, Memc[errmsg])

	    case IMAGE:
		call rdimspec (Memc[ifile], Memc[wfile], Memc[inform], badpix,
			        ilen, iwave, ispec)

		call wtabspec (Memc[ofile], Memc[outform], olength, 
			       ilen, iwave, ispec)

	    case TABLE:
		call rdtabspec (Memc[ifile], Memc[inform], badpix, 
				ilen, iwave, ispec)

		call wimspec (Memc[ofile], Memc[wfile], Memc[outform], 
			      olength, ilen, iwave, ispec)

	    case UNKNOWN:
		call sprintf (Memc[errmsg], SZ_LINE, badtype)
		call pargstr (Memc[ifile])
		call error (1, Memc[errmsg])
	    }

	    # Free arrays allocated by input routines

	    call mfree (iwave, TY_REAL)
	    call mfree (ispec, TY_REAL)
	}

	# Close file templates and free memory

	call imtclose (inlist)
	call imtclose (outlist)
	if (wavlist != NULL)
	    call imtclose (wavlist)

	call clssyntab
	call sfree (sp)
end
