include	<imhdr.h>
include	<error.h>
include <mwset.h>

# T_IMPROJECT -- project multidimensional image down one axis
#
# Giaretta,		   Task created.
# Phil Hodge, 13-Feb-1990  Rename from project to improject.
# Phil Hodge, 26-May-1994  Delete projlabel, and open output using section.
# Phil Hodge, 15-Apr-2011  Input image is now no longer open when proj_cre is
#			   called; in projectimn, use imggsd instead of imggsc;
#			   use double precision instead of single precision.
# Phil Hodge, 21-Mar-2012  In proj_cre, initialize the entire output image.

procedure t_improject ()

char	image1[SZ_FNAME], image2[SZ_FNAME]
int	projaxis
bool	average
bool	verbose
double	high, low

char	list[SZ_LINE]

pointer inlist			# Input image list
pointer outlist			# Output image list

pointer	in, out

pointer	immap()
pointer imtopen()
int	imtgetim(), clgeti()
double	clgetd()
int	clgwrd()
bool	clgetb()
bool	cut
int	pixtype
char	dtype[SZ_FNAME]
errchk	immap, clgetd, clgetb, clgeti, clgwrd

begin
	# Get task parameters.
	call clgstr ("input", list, SZ_LINE)
	inlist = imtopen (list)
	call clgstr ("output", list, SZ_LINE)
	outlist = imtopen (list)

	projaxis = clgeti("projaxis")
	average  = clgetb("average")
	high     = clgetd("highcut")
	low      = clgetd("lowcut")

	switch ( clgwrd( "pixtype", dtype, SZ_FNAME,
			",real,long,short,int,double") ) {
	case 1:
	    pixtype = TY_REAL
	case 2:
	    pixtype = TY_LONG
	case 3:
	    pixtype = TY_SHORT
	case 4:
	    pixtype = TY_INT
	case 5:
	    pixtype = TY_DOUBLE
	}

	verbose  = clgetb("verbose")
	if ( high == low)
	    cut = false
	else
	    cut = true

	# Loop through input images.

	while ((imtgetim (inlist,  image1, SZ_FNAME) != EOF) &&
	       (imtgetim (outlist, image2, SZ_FNAME) != EOF)) {

	    # Create output image using as template the input image name with
	    # a section.  This is to get appropriate coordinate parameters.
	    iferr (call proj_cre (image1, image2, projaxis, pixtype)) {
		call erract (EA_WARN)
		next
	    }

	    # Reopen the input and output images.
	    iferr (in = immap (image1, READ_ONLY, NULL)) {
		call erract (EA_WARN)
		next
	    }
	    if (IM_NDIM(in) < projaxis) {
		call printf(
			"image %s has dimension %d; cannot project axis %d\n")
		    call pargstr( image1 )
		    call pargi( IM_NDIM(in) )
		    call pargi( projaxis )
		call imunmap(in)
		next
	    }
	    out = immap (image2, READ_WRITE, NULL)

	    if (IM_NDIM(in) == 1) {
	        call projectim1( in, out, average, cut, low, high)

	    } else if (IM_NDIM(in) == 2) {
	        call projectim2( in, out, projaxis, average, cut, low, high)

	    } else {
		call projectimn( in, out, projaxis, average, cut, low, high)
	    }

	    # close files
	    call imunmap (out)
	    call imunmap (in)

	    if (verbose) {
		call printf(" projecting %s along axis %d ---> %s\n")
		    call pargstr( image1 )
		    call pargi( projaxis )
		    call pargstr( image2 )
		call flush( STDOUT)
	    }
	}

	call imtclose (inlist)
	call imtclose (outlist)
	
end

# proj_cre -- create output image
# This routine appends an image section string to the name of the input
# image (which therefore must not already include a section) and opens
# it read-only.  Then it opens the output image using the input as a
# template.  This is to set the output image size and to get the coordinate
# parameters set appropriately in the output image.  One pixel is written
# to the output, and then both input and output are closed.
#
# Phil Hodge, 27-May-1994  Subroutine created.
# Phil Hodge, 15-Apr-2011  Calling sequence modified; input image is no longer
#			   open when this function is called.

procedure proj_cre (image1, image2, projaxis, pixtype)

char	image1[ARB]	# i: name of input image
char	image2[ARB]	# i: name of output image
int	projaxis	# i: axis along which the projection is taken
int	pixtype		# i: data type for output image
#--
pointer sp
pointer input		# scratch for name of input with section
pointer in, out, x	# imhdr pointers, pointer to output data
pointer immap()
int	impnld()
long	v[IM_MAXDIM]	# for initializing output image
long	npix		# number of points along first output image axis
int	ndim		# dimension of input image
int	i
errchk	immap, impnld, imunmap

begin
	call smark (sp)
	call salloc (input, SZ_FNAME, TY_CHAR)

	# Open input image just to get the dimension.
	in = immap (image1, READ_ONLY, NULL)
	ndim = IM_NDIM(in)
	if (ndim < projaxis) {
	    call imunmap (in)
	    call error (1, "projaxis is too large")
	}
	call imunmap (in)

	do i = 1, IM_MAXDIM		# initialize for impnld
	    v[i] = 1

	# Append an image section to the input image name.

	call strcpy (image1, Memc[input], SZ_FNAME)

	call strcat ("[", Memc[input], SZ_FNAME)
	do i = 1, ndim {

	    if (i == projaxis)
		call strcat ("1:1", Memc[input], SZ_FNAME)
	    else
		call strcat ("*", Memc[input], SZ_FNAME)

	    if (i == ndim)
		call strcat ("]", Memc[input], SZ_FNAME)
	    else
		call strcat (",", Memc[input], SZ_FNAME)
	}

	# Create output image.
	in = immap (Memc[input], READ_ONLY, NULL)	# incl image section
	out = immap (image2, NEW_COPY, in)
	IM_PIXTYPE(out) = pixtype
	npix = IM_LEN(out, 1)
	while (impnld (out, x, v) != EOF) {
	    call amovkd (1.d0, Memd[x], npix)	# any value will do
	}
	call imunmap (out)
	call imunmap (in)

	call sfree (sp)
end

procedure projectim1( in, out, average, cut, low, high)

pointer	in		# i: input image
pointer	out		# i: output image
bool	average		# i: take average?
bool	cut		# i: do we really want to use the cuts?
double	high		# i: high cut
double	low		# i: low cut
#--
pointer	impl1d(), imgl1d()

begin
	call proj_avg( Memd[imgl1d(in)], IM_LEN(in, 1),
			average, cut, low, high, Memd[impl1d(out)] )
end

procedure projectim2( in, out, projaxis, average, cut, low, high)

pointer	in		# i: input image
pointer	out		# i: output image
int	projaxis	# i: axis to project along
bool	average		# i: take average?
bool	cut		# i: do we really want to use the cuts?
double	low		# i: low cut
double	high		# i: high cut
#--

long	npts, naxis1, i, j
pointer	data, outbuf, impl1d(), imgl2d(), nspt
double	temp

begin
	switch (projaxis) {
	case 1:
	    npts = IM_LEN (in, 1)
	    outbuf = impl1d (out)

	    do i = 1, IM_LEN(in, 2) {
		call proj_avg( Memd[imgl2d(in, i)], npts,
			average, cut, low, high, Memd[outbuf+i-1] )
	    }

	case 2:
	    npts = IM_LEN (in, 2)
	    naxis1 = IM_LEN( in, 1)
	    outbuf = impl1d (out)
	    if ( !cut) {
	        call aclrd( Memd[outbuf], naxis1)
	        do i = 1, npts
		    call aaddd( Memd[outbuf], Memd[imgl2d(in,i)], Memd[outbuf],
					naxis1)
		if ( average )
		    call adivkd(Memd[outbuf], double(npts), Memd[outbuf],
					naxis1)

	    } else {
		call salloc( nspt, naxis1, TY_LONG)
		call aclrl(Meml[nspt], naxis1)
		call aclrd( Memd[outbuf], naxis1)

		do i = 1, npts {
		    data = imgl2d( in, i)
		    do j = 0, naxis1-1 {
			temp = Memd[data + j]
			if ( !IS_INDEFD (temp) && temp <= high && temp >= low) {
			    Meml[nspt+j] = Meml[nspt+j] + 1
			    Memd[outbuf+j] = Memd[outbuf+j] + temp
			}
		    }
		}

		if ( average ) {
		    do i = 0, naxis1-1
		    	if ( Meml[nspt+i] == 0)
			    Memd[outbuf+i] = INDEF
			else
			    Memd[outbuf+i] = Memd[outbuf+i]/double(Meml[nspt+i])
		} else {
		    do i = 0, naxis1-1
		    	if ( Meml[nspt+i] == 0)
			    Memd[outbuf+i] = INDEF
		}
	    }
	}
end

procedure projectimn( in, out, projaxis, average, cut, low, high)

pointer	in		# i: input image
pointer	out		# i: output image
int	projaxis	# i: axis to project along
bool	average		# i: take average?
bool	cut		# i: do we really want to use the cuts?
double	high		# i: high cut
double	low		# i: low cut
#--

pointer	data, outbuf, imggsd(), impgsd()
long	vs[IM_MAXDIM],    ve[IM_MAXDIM]
long	vsout[IM_MAXDIM], veout[IM_MAXDIM]
long	one				# added by PEH, 1994 June 1
int	outndim, inaxref[IM_MAXDIM], inndim
long	ax, outax, inax, vsax
long	npts
bool	finished
errchk	imggsd, impgsd

begin
	one = 1
	call amovkl( one, vs, IM_MAXDIM)
	call amovkl( one, ve, IM_MAXDIM)
	call amovkl( one, vsout, IM_MAXDIM)
	call amovkl( one, veout, IM_MAXDIM)
	outndim = IM_NDIM( out)
	inndim  = IM_NDIM( in)
	npts    = IM_LEN(in, projaxis)

	vs[projaxis] = 1
	ve[projaxis] = IM_LEN(in, projaxis)

	# set output axes lengths	
	outax = 0
	do inax = 1, IM_NDIM(in) {
	    if ( inax != projaxis) {
		outax = outax + 1
		inaxref[outax] = inax
	        IM_LEN( out, outax) = IM_LEN( in, inaxref[outax])
	    }
	}

	# now read each column of input in turn until all read
	finished = false
	while ( !finished ) {

	    data   = imggsd( in, vs, ve, inndim)
	    outbuf = impgsd( out, vsout, veout, outndim)

	    # calculate the output pixel value
	    call proj_avg( Memd[data], npts,
			average, cut, low, high, Memd[outbuf] )

	    # increment vs, ve, vsout, veout, ensuring each component does
	    # not exceed corresponding axis length
	    for ( ax = 1; ax <= outndim ; ax = ax + 1) {
		vsout[ax] = vsout[ax] + 1
		if ( vsout[ax] <= IM_LEN(out, ax) ) {
		    veout[ax] = vsout[ax]
		    vsax      = inaxref[ax]
		    vs[vsax]  = vsout[ax]
		    ve[vsax]  = vs[vsax]
		    break
		} else {
		    if ( ax == outndim) {
			finished = true
			break
		    } else {
		        vsout[ax] = 1
		    	veout[ax] = 1
			vsax      = inaxref[ax]
		    	vs[vsax]  = 1
		    	ve[vsax]  = 1
		    }
		}
	    }
	    # The following is perhaps a more staightforward way of
	    # incrementing vs, ve, vsout, veout.
#bool	done
	#    ax = 1
	#    done = false
	#    while (!done) {
	#	if (ax > outndim) {
	#	    finished = true
	#	    break
	#	}
	#	vsax = inaxref[ax]
	#	if (veout[ax] >= IM_LEN(out, ax)) {
	#	    vsout[ax] = 1
	#	    veout[ax] = 1
	#	    vs[vsax]  = 1
	#	    ve[vsax]  = 1
	#	    ax = ax + 1
	#	} else {
	#	    vsout[ax] = vsout[ax] + 1
	#	    veout[ax] = vsout[ax]
	#	    vs[vsax]  = vsout[ax]
	#	    ve[vsax]  = vsout[ax]
	#	    done = true
	#	}
	#    }

	}
end

# PROJ_AVG -- calculate sum or average of an array, possibly with high
# and low cuts

procedure proj_avg( x, npts, average, cut, low, high, result)

double	x[ARB]		# i: input data
long	npts		# i: number of points
bool	average		# i: do we want average?
bool	cut		# i: do we want cut?
double	low		# i: low cut
double	high		# i: high cut
double	result		# o: resultant value
#--

long	nsum, i
double	temp, sum

begin
	nsum = 0
	sum  = 0.0

	if (cut) {
	    do i = 1, npts {
		temp = x[i]
		if ( !IS_INDEFD (temp) && temp <= high && temp >= low ) {
		    nsum = nsum + 1
		    sum  = sum + temp
		}
	    }
	} else {
	    do i=1, npts
		if ( !IS_INDEFD (x[i])) {
		    nsum = nsum + 1
		    sum = sum + x[i]
		}
	}

	if ( !average) {
	    if ( nsum == 0)
		result = INDEF
	    else
		result = sum
	} else {
	    if ( nsum == 0)
		result = INDEF
	    else
		result = sum/nsum
	}
end
