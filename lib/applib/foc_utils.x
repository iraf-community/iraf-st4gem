include	<syserr.h>
include <mach.h>
include	<imio.h>
include	<imhdr.h>
include	<ctype.h>
include	<time.h>
include	<foc.h>
define	SZ_PXFORMT	9	# for the string NORMAL or ZOOM

.help foc_utils
.nf______________________________________________________________________
FOC_UTILS -- a set out routines for dealing with peculiarities of FOC 
headers, in particular the photocathode coordinate system.

pointer = x_foc_immap( image, foc_coord, accmode, foc_template)
		opens FOC image and sets up extra pointers

call x_foc_out_label( foc_coord)  label output FOC image

int 	= x_foc_im_area( coord1, coord2, coordout, print_info)
 		checks requested coordinates against image size
		and generates new coordinates to ensure overlap

call xgfmodel( modelfile, rm , max_models)
		get model data from text file

call foc_date( str, maxchar)   get string with date for FOC use

bool 	= in_imagex( image, xcoord, ycoord, ldx, udx, ldy, udy)
		check if point and surroundings are within image

bool 	= foc_incheck( x, y, delx1, delx2, dely1, dely2, foc_coord)
		check if point is sufficiently within image


D. Giaretta, 1-Aug-1987  Original SPP version
Phil Hodge, 20-Sep-1989  In foc_incheck, add/subtract the borders.
Phil Hodge, 24-May-1990  Delete foc_agp and defines taken from stf.h.
Phil Hodge, 19-Dec-1990  Include add_headname in this file.
Phil Hodge, 10-Jul-1991  Don't get sampbeg or linebeg from image header.
Phil Hodge, 18-Jul-1991  x_foc_immap:  "flip" sampbeg but not sampoff.
Phil Hodge, 15-Feb-1993  x_foc_immap:  declare streq as bool, not int.
.endhelp_________________________________________________________________


# X_FOC_IMMAP -- open FOC image and get info about foc coordinates

pointer procedure x_foc_immap( image, foc_coord, accmode, foc_template)

char	image[SZ_FNAME]		# i: image name
pointer	foc_coord		# i: pointer to foc coord structure
int	accmode			# i: access mode
pointer	foc_template		# i: FOC structure template pointer
#--

pointer	immap(), im, template
char	pxformt[SZ_PXFORMT]	# NORMAL or ZOOM
int	samppln			# can be either naxis1 or naxis1/2
int	sampbeg
int	i_overlap
bool	zoom			# true if zoom mode
real	imgetr()
int	imaccf(), imgeti(), x_foc_im_area()
bool	streq()

errchk	immap, x_foc_im_area

begin
	if (foc_template != NULL) {
	    if (accmode == NEW_COPY)
	 	template = FOC_IMPT(foc_template)
	    else 
		template = NULL
	} else
	    template = NULL

	im = immap( image, accmode, template)

	FOC_IMPT(foc_coord)	= im
	FOC_ACCMODE(foc_coord)	= accmode
	FOC_NPIX1(foc_coord)	= IM_LEN(im, 1)
	FOC_NPIX2(foc_coord)	= IM_LEN(im, 2)
	
	FOC_INC1(foc_coord)	= IM_VSTEP(im, 1)
	FOC_INC2(foc_coord)	= IM_VSTEP(im, 2)
	if (accmode != NEW_COPY) {
	    FOC_SOFF1(foc_coord)	= IM_VOFF(im, 1)+FOC_INC1(foc_coord) - 1
	    FOC_SOFF2(foc_coord)	= IM_VOFF(im, 2)+FOC_INC2(foc_coord) - 1
	    iferr (FOC_SAMPOFF(foc_coord) = imgetr(im, "SAMPOFF"))
	        FOC_SAMPOFF(foc_coord) = 0.0
	    iferr (FOC_LINEOFF(foc_coord) = imgetr(im, "LINEOFF"))
	        FOC_LINEOFF(foc_coord) = 0.0
	    sampbeg = nint (FOC_SAMPOFF(foc_coord)) + 1
	    FOC_LINEBEG(foc_coord) = nint (FOC_LINEOFF(foc_coord)) + 1
	    iferr (samppln = imgeti (im, "samppln"))
		samppln = IM_LEN(im,1)
	    if (imaccf (im, "pxformt") == YES) {
		call imgstr (im, "pxformt", pxformt, SZ_PXFORMT)
		call strupr (pxformt)
		zoom = streq (pxformt, "ZOOM")
	    } else {
		zoom = false
	    }
	    if (zoom)
		FOC_SAMPBEG(foc_coord) = 1026 - 2 * samppln - sampbeg
	    else
		FOC_SAMPBEG(foc_coord) = 1026 - samppln - sampbeg

	} else {
	    FOC_SOFF1(foc_coord)	= 0
	    FOC_SOFF2(foc_coord)	= 0
	    i_overlap = x_foc_im_area( foc_template, NULL, foc_coord, true)
	    if (i_overlap == FOC_NO_OVERLAP || i_overlap == FOC_ERROR) {
		call error (0, " no overlap in images")
	    }
	}

	return (im)
end


# X_FOC_OUT_LABEL -- label output FOC image

procedure x_foc_out_label( foc_coord)

pointer	foc_coord		# i: foc coord structure
#--

begin
	call imaddi( FOC_IMPT(foc_coord), "SAMPBEG", FOC_SAMPBEG(foc_coord))
	call imaddi( FOC_IMPT(foc_coord), "LINEBEG", FOC_LINEBEG(foc_coord))
end

# X_FOC_IM_AREA -- checks requested coordinates against image size
# and generates new coordinates to ensure overlap

int procedure x_foc_im_area( coord1, coord2, coordout, print_info)

pointer	coord1		# i: pointer to input image coords structure
pointer coord2		# i: pointer to requested coord structure - or
                        #    pointer to second image coords structure
pointer	coordout	# i: pointer to output coord structure
bool	print_info	# i: should any info messages be printed?
#--

begin
	FOC_LINEBEG(coordout)	= FOC_LINEBEG(coord1) -
					FOC_SOFF2(coordout) + FOC_SOFF2(coord1)
	FOC_SAMPBEG(coordout)	= FOC_SAMPBEG(coord1) -
					FOC_SOFF1(coordout) + FOC_SOFF1(coord1)
	FOC_LINEOFF(coordout)	= FOC_LINEOFF(coord1) -
					FOC_SOFF2(coordout) + FOC_SOFF2(coord1)
	FOC_SAMPOFF(coordout)	= FOC_SAMPOFF(coord1) -
					FOC_SOFF1(coordout) + FOC_SOFF1(coord1)

	FOC_COORD_INIT(coordout)= FOC_COORD_INIT(coord1)
	FOC_PC_COORDS(coordout)	= FOC_PC_COORDS(coord1)

	return (FOC_COORD_OK)

end

# XGFMODEL -- get model data from text file

procedure xgfmodel( modelfile, rm , max_models)

char	modelfile[SZ_FNAME]	# i: model file
pointer	rm			# i: pointer to model data
int	max_models		# i: maximum number of models allowed
#--

pointer	open()
pointer	fdmodel
char	line[SZ_LINE]
long	ip, i_in_line
real	dval
int	ctor()
int	tot_models, nrows, ncols
int	getline()
pointer	currpt
bool	in_model
int	mod_dim, mod_count

errchk	getline, open, ctor

begin

	tot_models = 0
	ip = 0
	
	# read in the models from the reseau file if present
	if (modelfile[1] != EOS) {
	    # read in the text data into buffer 
	    i_in_line = 1
	    in_model  = false
	    mod_dim   = 1
	    fdmodel = open(modelfile, READ_ONLY, TEXT_FILE)
	    while ( getline(fdmodel, line) != EOF && tot_models <= max_models) {
	        i_in_line = 1
	        while ( ctor(line, i_in_line, dval) >= 1 && line[1] != '#' ) {
		    if (!in_model) {
			if (mod_dim == 1) {
			    tot_models = tot_models + 1
			    ncols = dval
			    MODEL_NAXIS1(rm, tot_models) = ncols
			    mod_dim = 2
			} else {
			    nrows = dval
			    MODEL_NAXIS2(rm,tot_models) = nrows
			    mod_dim   = 1
			    in_model  = true
			    mod_count = 1
			    call salloc ( currpt, nrows*ncols, TY_REAL)
			    MODEL_PTR(rm, tot_models) = currpt
			    MODEL_IN_RES(rm, tot_models) = false
			}
		    } else {
		        Memr[currpt+mod_count-1]=dval
		        ip = ip + 1
			mod_count = mod_count + 1
			if ( mod_count > nrows*ncols )
			    in_model = false
		    }
	        }
	    }
	    call close( fdmodel )

	    MODEL_NUM(rm) = tot_models 
	}


end


# ADD_HEADNAME -- add keywords for HEADNAME and DATANAME and DATE
# datamin and max are added only if the keyword is present

procedure add_headname (out, output, datamin, datamax)

pointer	out			# i: output file pointer
char	output[SZ_FNAME]	# i: output file name
real	datamin, datamax	# i: datamin, dtamax
#--

int	len, strlen()
char	name[SZ_FNAME]
char	temp[SZ_FNAME]
int	imaccf()

errchk	foc_date

begin
	call strcpy( output, name, SZ_FNAME)
	call strupr( name)

	len = strlen( name)
	# add default extension is one not given
	if ( name[len-3] != '.' ) 
	   call strcat( ".HHH", name, SZ_FNAME)
	# get new length
	len = strlen( name)

	# headname and dataname
	name[len] = 'H'
	call imastr( out, "HEADNAME", name)
	name[len] = 'D'	
	call imastr( out, "DATANAME", name)

	# date
	call foc_date( temp, SZ_FNAME)
	call imastr( out, "DATE", temp)

	# min and max, if keywords present
	if ((imaccf (out, "DATAMIN") == YES) &&
	    (imaccf (out, "DATAMAX") == YES)) {
	    call imputr (out, "i_minpixval", datamin)
	    call imputr (out, "i_maxpixval", datamax)
	    IM_LIMTIME(out) = IM_MTIME(out) + 1
	}
end

# FOC_DATE -- creates date string for FOC data file use

procedure foc_date( str, maxchar)

char	str[maxchar]
int	maxchar
#--

long	clktime(), ltime
long	ptime[LEN_TMSTRUCT]

begin
	ltime = clktime( 0)

	call brktime( ltime, ptime)

	call sprintf( str, maxchar, "%02d/%02d/%02d")
	    call pargi( TM_MDAY( ptime) )
	    call pargi( TM_MONTH( ptime) )
	    call pargi( TM_YEAR( ptime) -1900)
end

# IN_IMAGEX -- check if point and surroundings are within image

bool procedure in_imagex( image, xcoord, ycoord, ldx, udx, ldy, udy)

pointer	image		# i: image to be checked
long	xcoord		# i: x coordinate of center
long	ycoord		# i: y coordinate of center	
long	ldx, udx	# i: lower and upper deltas to xcoord
long	ldy, udy	# i: lower and upper deltas to ycoord
#--

begin
	if ( xcoord - ldx >= 1 &&
	     xcoord + udx <= IM_LEN(image, 1) &&
	     ycoord - ldy >= 1 &&
	     ycoord + udy <= IM_LEN(image, 2)     )
	    return ( true )
	else
	    return ( false )

end

# FOC_INCHECK -- check if point is sufficiently within image

bool procedure foc_incheck( x, y, delx1, delx2, dely1, dely2, foc_coord)

real	x, y			# i: point coords
real	delx1, delx2		# i: borders left and right, x
real	dely1, dely2		# i: borders top and bottom, y
pointer	foc_coord		# i: foc coord structure
#--

begin
	if (x-delx1 < real (FOC_SAMPBEG(foc_coord) ) ||
	    x+delx2 > real (FOC_NPIX1(foc_coord) + FOC_SAMPBEG(foc_coord)-1) ||
	    y-dely1 < real (FOC_LINEBEG(foc_coord)) ||
	    y+dely2 > real (FOC_NPIX2(foc_coord) + FOC_LINEBEG(foc_coord)-1))

	    return (false)

	else

	    return (true)
end
