include	<imhdr.h>
define	REPLACE		1
define	ADD		2
define	MULTIPLY	3

# iminsert -- Insert a small image into a larger one in a number
# of different ways.
#
# Dave Giaretta		 Task created.
# Phil Hodge, 19-Sep-1989  Flush output after copying from first input image;
#			number of lines to copy is IM_LEN(im2,2).
# Phil Hodge, 13-Feb-1990  Change name from insertx to iminsert.
# Phil Hodge,  6-Aug-1993  After copying, read from output; also imflush.

procedure t_iminsert()

char	input1[SZ_FNAME]	# input large image into which input2 is 
				#   inserted
char	input2[SZ_FNAME]	# small 2D image to be inserted into input1
char	output[SZ_FNAME]	# output image
char	option[SZ_LINE]		# method of insertion - replace, add or multiply
char	coordfile[SZ_FNAME]	# file with coordinates for insertion
int	offset1, offset2	# offsets for coordinates
#--

pointer	in1, in2, out, immap()
pointer	imgs2r(), imps2r(), xpt, ypt, buff1, buff2, buffo
long	v_i[IM_MAXDIM], v_o[IM_MAXDIM]
long	xcoord_type()
int	imgnlr(), impnlr()
int	npos, kwindex, clgwrd(), clgeti()
int	line, line1, line2, pix11, pix12, pix21, pix22, delt, i
char	imtemp[SZ_FNAME]

begin

	call clgstr( "input1", input1, SZ_FNAME)	
	call clgstr( "input2", input2, SZ_FNAME)	
	call clgstr( "output", output, SZ_FNAME)	
	kwindex = clgwrd( "option", option, SZ_LINE, ",replace,add,multiply" )
	call clgstr( "coordfile", coordfile, SZ_FNAME)
	offset1 = clgeti("offset1")
	offset2 = clgeti("offset2")

	# watch out for output = input1
	call xt_mkimtemp ( input1, output, imtemp, SZ_FNAME)

	# open files
	in1 = immap( input1, READ_ONLY, 0)
	in2 = immap( input2, READ_ONLY, 0)
	out = immap( output, NEW_COPY,  in1)

	# read coordinates
	npos = xcoord_type( coordfile, "entry", "xcol", "ycol", xpt, ypt)

	# add offsets
	do i = 1, npos {
	    Meml[xpt + i-1] = Meml[xpt +i-1] + offset1
	    Meml[ypt + i-1] = Meml[ypt +i-1] + offset2
	}

	# copy input1 to output
	call amovkl( long(1), v_i, IM_MAXDIM)	
	call amovkl( long(1), v_o, IM_MAXDIM)	
	while( imgnlr( in1, buff1, v_i) != EOF && 
	       impnlr( out, buffo, v_o) != EOF    ) 
	    call amovr( Memr[buff1], Memr[buffo], IM_LEN( in1, 1) )

	call imflush (out)

	# Done with first input image.
	call imunmap( in1)

	# now insert sections for each position
	for ( i = 1; i<= npos; i=i+1) {

	    # first adjust starting position to use overlap area
	    line1 = Meml[ypt+i-1]
	    if ( line1 < 1 ) {
		line2 = 2 - line1
		line1 = 1
	    } else
		line2 = 1

	    pix11 = Meml[xpt+i-1] 
	    if ( pix11 < 1) {
		pix21 = 2 - pix11
		pix11 = 1
	    } else 
		pix21 = 1
	    
	    pix12 = min ( pix11 + IM_LEN( in2, 1) - 1, IM_LEN( out, 1) )
	    delt  = min ( pix12 - pix11 + 1, IM_LEN( in2, 1) - pix21 + 1)
	    pix12 = pix11 + delt - 1
	    pix22 = pix21 + delt - 1

	    # check if there is some overlap in pixels
	    if ( pix11 >= 1 && pix11 <= pix12 &&
		 pix21 >= 1 && pix21 <= pix22     ) {

		# This loop will only proceed if there are some overlap lines.
		# IM_LEN( in2, 1) changed to IM_LEN( in2, 2) by PEH 9/19/89
		# buff1 is read from out instead of from in1; PEH 8/6/93
	        for ( line = line2; 
		      line1 <= IM_LEN(out, 2) && line <= IM_LEN( in2, 2); 
		      line = line + 1) {
		    buff1 = imgs2r( out, pix11, pix12, line1, line1)
		    buff2 = imgs2r( in2, pix21, pix22, line , line )
		    buffo = imps2r( out, pix11, pix12, line1, line1)
		    switch (kwindex) {
		    case REPLACE:
			call amovr( Memr[buff2], Memr[buffo], pix12-pix11+1)
		    case ADD:
			call aaddr( Memr[buff1], Memr[buff2], Memr[buffo], 
								pix12-pix11+1 )
		    case MULTIPLY:
			call amulr( Memr[buff1], Memr[buff2], Memr[buffo], 
								pix12-pix11+1 )
		    }

		    # select next line
		    line1 = line1 + 1
		}
	    }

	    call imflush (out)			# added by PEH 8/6/93
	}		


	call imunmap( in2)
	call imunmap( out)

	call xt_delimtemp( output, imtemp)

end	    
