include <imhdr.h>
include <imset.h>
include <error.h>
include <mach.h>
include <ctype.h>
include <fset.h>        # FIO
include <tbset.h>       # TBTABLES

define  DEF_ERR		1.0		# Default error bar size.
define  SZ_BUF          512             # Initial pixel buffer size
define  LIST_OP         1
define  IMAGE_OP        2
define  TABLE_OP        3
define  RE_TRY          -99             # error flag value


# RC_RDCURVES -- Given a operand list as input, read in each of the 
# referenced tables, lists and/or image sections, in turn. Initialize 
# things the first time around.
#
#
#
# Original code: from original version of gfit1d by Chris Biemesderfer.
#
# Revision history
# ----------------
#
# 22-Feb-89   Takes into account physical units descriptor parameters 
#             in image headers (I. Busko).
# 25-Jun-90   Properly handles image sections (IB).
# 13-Jan-93   Fixed bug in section handling (IB).
# 24-May-94   Changed weighting from 0.5 to 1 for equal weights (J-C Hsu)
# 03-Jun-94   Fixed the bug of x array's index off by 1 (JCH)
# 03-Sep-97   Uses ksection to access FITS extensions (IB)
# 16-Sep-97   Accesses arrays in tables using selectors (IB)
# 19-Sep-97   Read from arrays in multiple rows (IB)
# 19-Sep-97   Fix memory leak on data vectors (IB)
# 23-Sep-97   Add support for row selectors with scalar columns (IB)
# 30-Jan-98   Uses clindex to access FITS extensions (IB)
# 28-Apr-00   Turn off log mode when projecting on other than disp axis (IB)


int procedure rc_rdcurves (oplist, image, ft, fty, nnames, word, x, y, size,
                           npix, axis, rdmarks, istat)

char    oplist[SZ_LINE]         # Operand list
char    image[SZ_FNAME]         # Image name, with section appended
pointer ft                      # file pointer
pointer fty                     # y data file pointer
int     nnames                  # number of names in input list
char    word[SZ_FNAME,4]        # Input line words
pointer x                       # Pointer to x vector
pointer y                       # Pointer to y vector
pointer size                    # Pointer to vector of marker sizes
int     npix                    # Number of values per vector
int     axis                    # Axis for projection
bool    rdmarks                 # Read mark sizes from list?
int     istat                   # OK?

#--
char    table[SZ_FNAME], column[SZ_COLNAME]
char    root[SZ_FNAME], rowselect[SZ_FNAME]
char    root2[SZ_FNAME], rowselect2[SZ_FNAME]
char    errcol[SZ_COLNAME], colselect[SZ_FNAME]
char    xtable[SZ_FNAME], xcolumn[SZ_COLNAME]
char    ytable[SZ_FNAME], ycolumn[SZ_COLNAME]
int     cdim1, cdim2, iw

int     rc_rdcurve(), imtopen(), imtgetim(), nscan()
int     fntopn(), fntgfn()
int     rc_rd1col(), rc_rd2col(), rc_rdxycol()
int     rc_rd1s(), rc_rdxys(), rc_rd2s(), rc_cdim()

errchk	rc_rdcurve, rc_rd1col, rc_rdxycol, rc_rd2col
errchk	rc_rd1s, rc_rdxys, rc_rd2s, rc_cdim

begin
        istat = NULL

	# Fix memory leak inherited from original code.
	if (x    != NULL) call mfree (x,    TY_REAL)
	if (y    != NULL) call mfree (y,    TY_REAL)
	if (size != NULL) call mfree (size, TY_REAL)

        if ( ft == NULL ) {
            # Parse the input
            call sscan (oplist)
            do iw = 1, 4
                call gargwrd (word[1,iw], SZ_FNAME)

            nnames = nscan()
        }

        switch (nnames ) {
            # How many words in input line?

        case 1:
            # Image or list file;  Read all curves into memory.

            if ( ft == NULL )
                ft = imtopen (word)
            if (imtgetim (ft, image, SZ_FNAME) == EOF) {
                call imtclose (ft)
                return ( EOF )
            }

            iferr {
                npix = rc_rdcurve (image,
                        x, y, size, axis, rdmarks, istat)
            } then {
                    call erract (EA_WARN)
                    istat = RE_TRY
                    return ( OK)
            }

        case 2:
            # Single table column

            call strcpy (word[1,2], column, SZ_COLNAME)
            if (rdmarks)
                call clgstr ("errcolumn", errcol, SZ_COLNAME)

	    # Must get rid of selectors, otherwise IRAF file name
	    # template routines bomb out.
	    call rdselect (word[1,1], root, rowselect, colselect, SZ_FNAME)

            if ( ft == NULL )
                ft = fntopn (root)
            if  (fntgfn (ft, table, SZ_FNAME) == EOF) {
                call fntcls (ft)
                return ( EOF )
            }
            call strcpy (word[1,1], image, SZ_FNAME)

	    # Get column dimensionality.
	    cdim1 = rc_cdim (root, column)

	    if (cdim1 == 1) {

	        # Scalar table.
                iferr {
                    npix = rc_rd1col (table, rowselect, column, rdmarks, 
                                      errcol, istat, x, y, size )
                } then {
	            call erract (EA_WARN)
                    istat = RE_TRY
                    return ( OK )
                }
	    } else {

	        # 3-D table.
	        iferr {
                    npix = rc_rd1s (table, rowselect, column, rdmarks, 
                                    errcol, istat, x, y, size)
                } then {
	            call erract (EA_WARN)
	            istat = RE_TRY
	            return ( OK )
                }
	    }

        case 3:
            # Single table;  X and Y column

            call strcpy (word[1,2], xcolumn, SZ_COLNAME)
            call strcpy (word[1,3], ycolumn, SZ_COLNAME)
            if (rdmarks)
                call clgstr ("errcolumn", errcol, SZ_COLNAME)

	    # Must get rid of selectors, otherwise IRAF file name
	    # template routines bomb out.
	    call rdselect (word[1,1], root, rowselect, colselect, SZ_FNAME)

            if ( ft == NULL )
                ft = fntopn (root)
            if (fntgfn (ft, table,  SZ_FNAME) == EOF) {
                call fntcls (ft)
                return ( EOF )
            }
            call strcpy (word[1,1], image, SZ_FNAME)

	    # Get column dimensionality.
	    cdim1 = rc_cdim (root, xcolumn)
	    cdim2 = rc_cdim (root, ycolumn)
	    if (cdim1 != cdim2)
	        call error (0, "Columns with incompatible dimensionality.")

	    if (cdim1 == 1) {

	        # Scalar table.
	        iferr {
                    npix = rc_rdxycol (table, rowselect, xcolumn, ycolumn, 
                                       rdmarks, errcol, istat, x, y, size )
                } then {
                    call erract (EA_WARN)
                    istat = RE_TRY
                    return ( OK )
                }
	    } else {

	        # 3-D table.
	        iferr {
                    npix = rc_rdxys (table, rowselect, xcolumn, ycolumn,
                                     rdmarks, errcol, istat, x, y, size)
                } then {
	            call erract (EA_WARN)
	            istat = RE_TRY
	            return ( OK )
                }
	    }

        case 4:
            # X table and column;  Y table and column

            call  strcpy (word[1,2], xcolumn, SZ_COLNAME)
            call  strcpy (word[1,4], ycolumn, SZ_COLNAME)
            if (rdmarks)
                call clgstr ("errcolumn", errcol, SZ_COLNAME)

	    # Must get rid of selectors, otherwise IRAF file name
	    # template routines bomb out.
	    call rdselect (word[1,1], root,  rowselect,  colselect, SZ_FNAME)
	    call rdselect (word[1,3], root2, rowselect2, colselect, SZ_FNAME)

            if ( ft == NULL )
                ft = fntopn (root)
            if ( fty == NULL )
                fty = fntopn (root2)
            if (fntgfn (ft , xtable, SZ_FNAME) == EOF ||
                fntgfn (fty, ytable, SZ_FNAME) == EOF   ) {
                call fntcls (ft )
                call fntcls (fty)
                return ( EOF)
            }
            call strcpy (word[1,1], image, SZ_FNAME)

	    # Get column dimensionality.
	    cdim1 = rc_cdim (root, xcolumn)
	    cdim2 = rc_cdim (root2, ycolumn)
	    if (cdim1 != cdim2)
	        call error (0, "Columns with incompatible dimensionality.")

	    if (cdim1 == 1) {

	        # Scalar tables.
                iferr {
                    npix = rc_rd2col (xtable, xcolumn, ytable, ycolumn,
                                     rowselect, rowselect2, rdmarks, 
                                     errcol, istat, x, y, size)
                } then {
                    call erract (EA_WARN)
                    istat = RE_TRY
                    return ( OK )
                }
	    } else {

	        # 3-D table.
	        iferr {
                    npix = rc_rd2s (xtable, xcolumn, ytable, ycolumn,
                                    rowselect, rowselect2, rdmarks, 
                                    errcol, istat, x, y, size)
                } then {
	            call erract (EA_WARN)
	            istat = RE_TRY
	            return ( OK )
                }
	    }
        }

        return ( npix)

end


# RC_RDCURVE -- Read a curve into memory.  The operand may specify either
# list or image input; we determine which and then call the appropriate
# input routine to access the data.

int procedure rc_rdcurve (operand, x, y, size, axis, rdmarks, istat)

char    operand[SZ_LINE]        # Operand to be fitted
pointer x, y, size              # Pointers to x, y and size arrays
int     axis                    # Axis of image projection
bool    rdmarks                 # Read marks from list?
int     istat                   # X or Y errors?

int     rc_rdlist2(), rc_rdimage2(), rc_optype(), npix
errchk  rc_rdlist2, rc_rdimage2, rc_optype

begin
        if (rc_optype (operand) == LIST_OP)
            npix = rc_rdlist2 (operand, x, y, size, rdmarks, istat)
        else
            npix = rc_rdimage2 (operand, x, y, size, axis)

        return ( npix)
end


# RC_OPTYPE -- Determine whether the operand argument is an image section
# or a list.  If the string is STDIN, it is a list; if a subscript is
# present, it is an image; otherwise we must test whether or not it is a
# binary file and make the decision based on that.

int procedure rc_optype (operand)

char    operand[SZ_LINE]        # Operand to be fitted
int     first, last, ip
int     access(), strncmp()
pointer im, immap()

errchk  immap

begin
        # Strip off any whitespace at the beginning or end of the string.
        for (ip=1;  IS_WHITE(operand[ip]);  ip=ip+1)
            ;
        first = ip
        for (last=ip;  operand[ip] != EOS;  ip=ip+1)
            if (!IS_WHITE(operand[ip]))
                last = ip

        if (first == last) {
            return (LIST_OP)
        } else if (strncmp (operand[first], "STDIN", 5) == 0) {
            return (LIST_OP)
        } else if (operand[last] == ']') {
            return (IMAGE_OP)
        } else if (access (operand, 0, TEXT_FILE) == YES) {
            iferr ( im = immap (operand, READ_ONLY, 0) ) {
                return (LIST_OP)
            } else {
                return (IMAGE_OP)
            }
        } else {
            return (IMAGE_OP)
        }
end


define  SZ_STRING       10

# RC_RDIMAGE2 -- Read an image section and compute the projection about
# one dimension, producing x and y vectors as output. Output vectors are
# returned in image units, taking into account the relevant parameters
# of the image header.

int procedure rc_rdimage2 (imsect, x, y, size, axis )

char    imsect[SZ_FNAME]        # Image section to be fitted
pointer x, y, size              # Pointer to x, y and size vector
int     axis                    # Axis about which the projection is to be taken

char    section[SZ_FNAME]
int     npix, dcflag, i, j, secaxis
long    x1, x2, stepx
real    crval, cdelt, crpix, bscale, bzero
pointer im

int     strlen(), getaxis()
pointer immap()

errchk  immap, rc_projection, malloc, rc_imparam

begin
        # Gets parameters from image header
        call rc_imparam (imsect, axis, npix, bscale, bzero, crval,
                         crpix, cdelt, dcflag)

        # Processing of image section specification.
        call imgsection (imsect, section, SZ_FNAME)
        if (strlen(section) > 0) {
            # Find projection axis specs.
            i = 1
            if (axis > 1) {
                do j = 2, axis {
                    while (section[i] != ',')
                        i = i + 1
                }
            }
            i = i + 1
            # Decode axis specs.
            while (IS_WHITE(section[i]))
                i = i + 1
            if (section[i] != '*') {
                call im_decode_subscript (section, i, x1, x2, stepx)
                npix = (x2 - x1) / stepx + 1
            }
            else {
                x1    = long(1)
                stepx = long(1)
            }
        } else {
            x1    = long(1)
            stepx = long(1)
        }

        # Generates y array from axis projection and constant error bar array.
        call malloc (x, npix, TY_REAL)
        call malloc (y, npix, TY_REAL)
        call malloc (size, npix, TY_REAL)

        # Get projection. First, compute relative section axis.
        secaxis = axis - getaxis (section) + 1
        im = immap (imsect, READ_ONLY, 0)
        call rc_projection (im, Memr[y], npix, secaxis)
        call imunmap (im)

        call amovkr (DEF_ERR, Memr[size], npix)

        # Generates x and y arrays in physical units.
        call rc_unit (Memr[x], Memr[y], npix, bscale, bzero, crval, crpix,
                      cdelt, dcflag, x1, stepx)

        return (npix)
end


# RC_IMPARAM -- Read image header parameters. Axis is the oiginal image
# axis, not a section axis.
#
# Installed: 22-Feb-89 - I.Busko

procedure rc_imparam (image, axis, npix, bscale, bzero, crval, crpix, cdelt,
                      dcflag)

char    image[ARB]              # i: Image name
int     axis                    # i: Image axis
int     npix                    # o: Number of pixels/line
real    bscale                  # o: Fits parameters
real    bzero
real    crval
real    crpix
real    cdelt
int     dcflag                  # o: Dispersion correction flag

char    name[SZ_LINE], ksection[SZ_LINE], section[SZ_LINE]
int     dispaxis, clindex, clsize
char    astring[SZ_STRING]
pointer im

int     imgeti()
real    imgetr()
pointer immap()

errchk  imgeti, imgetr, immap

begin
        # Get rid of image section spec.
#       call imgcluster (image, name, SZ_LINE)
	# imgcluster gets confused with FITS specs as "[SCI,1]" (9/3/97, IB)
	call imparse (image, name, SZ_LINE, ksection, SZ_LINE, section, 
                      SZ_LINE, clindex, clsize)
        if (clindex != -1) {
	    call sprintf (astring, SZ_STRING, "[%d]")
	    call pargi (clindex)
            call strcat (astring, name, SZ_LINE)
	}
        call strcat (ksection, name, SZ_LINE)

        im = immap (name, READ_ONLY, 0)

        if (axis < 1 || axis > IM_NDIM(im))
            call error (2, "Non existent axis")
        npix = IM_LEN (im, axis)

        # Determine if image is a spectrum which has been dispersion corrected.
        # If no DC-FLAG parameter is found, assume projection axis is in linear
        # wavelength units.
        iferr (dcflag = imgeti (im, "DC-FLAG"))
            dcflag = 0

        # Try to see if dispersion axis is marked.
        if (dcflag > -1) {
            iferr (dispaxis = imgeti (im, "DISPAXIS"))
                dispaxis = axis
            if (dispaxis != axis) {
                call eprintf (
                    "Warning: projection axis is not dispersion axis.")
	        # must turn off logarithmic mode
	        dcflag = 0
	    }
        }

        # Try to get axis parameters.
        crval = 1.
        crpix = 1.
        cdelt = 1.
        if (dcflag > -1) {
            call sprintf (astring, SZ_STRING, "CRPIX%1d")
                call pargi (axis)
            iferr (crpix = imgetr (im, astring))
                crpix = 1
            iferr (crval = imgetr (im, "W0")) {
                call sprintf (astring, SZ_STRING, "CRVAL%1d")
                    call pargi (axis)
                iferr (crval = imgetr (im, astring))
                    crval = 1.
            }
            iferr (cdelt = imgetr (im, "WPC")) {
                call sprintf (astring, SZ_STRING, "CD%1d_%1d")
                    call pargi (axis)
                    call pargi (axis)
                iferr (cdelt = imgetr (im, astring)){
                    call sprintf (astring, SZ_STRING, "CDELT%1d")
                        call pargi (axis)
                    iferr (cdelt = imgetr (im, astring))
                        cdelt = 1.
                }
            }
        } else
            call eprintf ("Warning: not dispersion corrected.")

        # Try to get pixel scale
        iferr (bscale = imgetr (im, "BSCALE"))
            bscale = 1.
        iferr (bzero = imgetr (im, "BZERO"))
            bzero = 0.

        call imunmap (im)
end


# RC_UNIT -- Converts raw pixel units to physical units. Y array must
# contain on entry raw pixel contents. X array is generated from FITS
# parameters and eventual image section spec.
#
# Installed: 22-Feb-89 - I.Busko

procedure rc_unit (x, y, npix, bscale, bzero, crval, crpix, cdelt, dcflag,
                   x1, stepx)

real    x[ARB]                  # o: Independent variable
real    y[ARB]                  # io: Dependent variable
int     npix                    # i: Number of pixels/line
real    bscale                  # i: Fits parameters
real    bzero
real    crval
real    crpix
real    cdelt
int     dcflag                  # i: Dispersion correction flag
long    x1, stepx               # i: Image section specs

int     i

begin
        # Converts y array.

        call amulkr (y, bscale, y, npix)
        call aaddkr (y, bzero, y, npix)

        # Generates x array.
        do i = 1, npix {
	    
	    # fix the index bug (from i to i-1), JC Hsu, 6/3/94
            x[i] = crval + cdelt * (((i-1) * int(stepx) + int(x1)) - crpix)

            if (dcflag == 1)                            # Log wavelength
                x[i] = 10. ** x[i]
        }
end


# RC_RDLIST2 -- Read a list of two dimensional data pairs into two type
# real arrays in memory.  Return pointers to the arrays and a count of the
# number of pixels.  If mark sizes are to be read from the input list,
# a third array of mark sizes is returned.  Mark sizes can only be given
# in two column (x,y) mode, with the mark size given as a third column.

int procedure rc_rdlist2 (fname, x, y, size, rdmarks, istat)

char    fname[SZ_FNAME]         # Name of list file
pointer x, y, size              # Pointers to x, y and size vectors
bool    rdmarks                 # Read markers from file?
int     istat                   # X or Y errors?

int     buflen, n, fd, ncols, lineno
pointer sp, lbuf, ip
real    xval, yval, szmark1, szmark2
int     getline(), nscan(), open()
errchk  open, sscan, getline, malloc

begin
        call smark (sp)
        call salloc (lbuf, SZ_LINE, TY_CHAR)

        fd = open (fname, READ_ONLY, TEXT_FILE)

        buflen = SZ_BUF
        iferr {
            call malloc (x, buflen, TY_REAL)
            call malloc (y, buflen, TY_REAL)
            call malloc (size, buflen, TY_REAL)
        } then
            call erract (EA_FATAL)

        n = 0
        ncols = 0
        lineno = 0
        szmark1 = 0.5
        szmark2 = 0.5

        while (getline (fd, Memc[lbuf]) != EOF) {
            # Skip comment lines and blank lines.
            lineno = lineno + 1
            if (Memc[lbuf] == '#')
                next
            for (ip=lbuf;  IS_WHITE(Memc[ip]);  ip=ip+1)
                ;
            if (Memc[ip] == '\n' || Memc[ip] == EOS)
                next

            # Decode the points to be fitted.
            call sscan (Memc[ip])
                call gargr (xval)
                call gargr (yval)
                if (rdmarks) {
                    # Find the marker or error bar size(s)
                    call gargr (szmark1)
                    call gargr (szmark2)
                }

            # The first line determines whether we have an x,y list or a
            # y-list.  It is an error if only one value can be decoded when
            # processing a two column list.

            if (ncols == 0 && nscan() > 0)
                ncols = nscan()

            switch (nscan()) {
            case 0:
                call eprintf ("no args; %s, line %d: %s\n")
                    call pargstr (fname)
                    call pargi (lineno)
                    call pargstr (Memc[lbuf])
                next
            case 1:
                if (ncols == 2) {
                    call eprintf ("only 1 arg; %s, line %d: %s\n")
                        call pargstr (fname)
                        call pargi (lineno)
                        call pargstr (Memc[lbuf])
                    next
                } else {
                    yval = xval
                    xval = n + 1.0
                }
            case 2:
                if (rdmarks) {
#                   call eprintf ("no mark size field; %s, line %d: %s\n")
#                       call pargstr (fname)
#                       call pargi (lineno)
#                       call pargstr (Memc[lbuf])
#		    change szmark1 from DEF_ERR/2. to DEF_ERR, JC Hsu 5/24/94
                    szmark1 = DEF_ERR
                }
            case 3:
                if (ncols == 4 && rdmarks) {
#                   call eprintf ("no upper error field; %s, line %d: %s\n")
#                       call pargstr (fname)
#                       call pargi (lineno)
#                       call pargstr (Memc[lbuf])
                    szmark2 = szmark1
                }
            }

           # only keep values if there are no INDEF values
            if ( !IS_INDEFR (xval) && !IS_INDEFR (yval) ) {
                n = n + 1
                if (n > buflen) {
                    buflen = buflen + SZ_BUF
                    call realloc (x, buflen, TY_REAL)
                    call realloc (y, buflen, TY_REAL)
                    call realloc (size, buflen, TY_REAL)
                }

                if (rdmarks) {
                    if (istat == 1) {
                        # Errors in X
                        Memr[size+n-1] = szmark1 + szmark2
                        xval = xval + (szmark2 - szmark1)/2.0
                    } else if (istat == 2) {
                        # Errors in Y
                        Memr[size+n-1] = szmark1 + szmark2
                        yval = yval + (szmark2 - szmark1)/2.0
                    } else
                        # Marker
                        Memr[size+n-1] = szmark1
                }
                Memr[x+n-1] = xval
                Memr[y+n-1] = yval
            }
        }

        call realloc (x, n, TY_REAL)
        call realloc (y, n, TY_REAL)
        call realloc (size, n, TY_REAL)

        call close (fd)
        call sfree (sp)
        return (n)
end


# RC_RD1COL -- Read fit data from a scalar column of a table.
#
#
# Modified: 23-Sep-97 - I.Busko  -  To process row selectors.

int procedure rc_rd1col (table, rowselect, column, rdmarks, 
                         errcol, istat, x, y, size)

char    table[SZ_FNAME]         # Table name
char    rowselect[SZ_FNAME]     # Row selector as parsed by rdselect
char    column[SZ_COLNAME]      # Data column name
bool    rdmarks                 # Read errors from file?
char    errcol[SZ_COLNAME]      # Error column name
int     istat                   # X or Y errors?
pointer x, y, size              # Pointers to x, y and size vectors

#--
pointer tdp                     # Pointers to table descriptor
pointer cdp, scdp               # Pointers to column descriptors
pointer	rcode			# Row selector structure
int     numrows                 # Maximum number of rows
int	nsrows			# Actual number of selected rows
int	irow, numcols, k
bool	null[1]

pointer tbtopn(), trsopen()
int     tbpsta(), selrows()
bool	trseval()

errchk	tbtopn, trsopen, malloc

begin
	numcols = 1
	tdp = tbtopn (table, READ_ONLY, 0)

	numrows = tbpsta (tdp, TBL_NROWS)

	# Find column.
	call tbcfnd (tdp, column, cdp, numcols)
	if (cdp <= 0)
	    call error (0, "Cannot find Y column.")
	if (rdmarks && (errcol[1] != EOS)) {
	    call tbcfnd (tdp, errcol, scdp, numcols)
	    if (scdp <= 0)
	        call error (0, "Cannot find error column.")
	}

	# Compute number of selected rows.
	nsrows = selrows (tdp, rowselect)
	if (nsrows == 0)
	    call error (0, "No rows selected.")

	# Allocate space for the fit arrays.
	iferr {
	    call malloc (x,    nsrows, TY_REAL)
	    call malloc (y,    nsrows, TY_REAL)
	    call malloc (size, nsrows, TY_REAL)
	    call amovkr (1.0, Memr[size], nsrows)
	} then
	    call erract (EA_FATAL)

	# Scan rows and read from selected ones.
	rcode = trsopen (tdp, rowselect)
	k = 0
	do irow = 1, numrows {

	    if (trseval (tdp, irow, rcode)) {

	        # Fill in the independent variable array.
	        Memr[x+k] = real (k+1)

	        # Read Y data.
	        call tbcgtr (tdp, cdp, Memr[y+k], null, irow, irow)

	        # If selected, read error. This only supports a single,
	        # assumed symmetric, error bar. This code is executed
	        # only when the row selector selects a subsample of the
	        # full row range in the table.
	        if (rdmarks && (errcol[1] != EOS) && (nsrows != numrows))
	            call tbcgtr (tdp, scdp, Memr[size+k], null, irow, irow)

	        # Bump index.
	        k = k + 1
	    }
	}

	# This is old code kept in place for backwards compatibility.
	# It is capable of reading asymmetric error bars. This code
	# executes only when the row selector assumes the default value, 
	# that is, it selects all rows in the table.
	if (rdmarks && (errcol[1] != EOS) && (nsrows == numrows))
	    call rc_ggrser (tdp, errcol, numrows, Memr[x], Memr[y], 
                            Memr[size], istat)

	call trsclose (rcode)
	call tbtclo (tdp)
	return (nsrows)
end


# RC_GGRSER -- Find the location and size of error bars from one or two
# columns in the table

procedure rc_ggrser (tdp, errcol, numrows, x, y, size, istat)

pointer tdp
char    errcol[SZ_FNAME]
int     numrows
real    x[numrows]
real    y[numrows]
real    size[numrows]
int     istat

char    lerrcol[SZ_COLNAME], herrcol[SZ_COLNAME]
pointer secp, lecp, hecp
pointer null                    # Pointer to null column values array
pointer lerr, herr, offs
char    errmsg[SZ_LINE]         # Error message

int     nscan()

begin
        # Read the errors column name(s)
        call malloc (null, numrows, TY_BOOL)
        call sscan (errcol)
            call gargwrd (lerrcol, SZ_COLNAME)
            call gargwrd (herrcol, SZ_COLNAME)

        switch (nscan ()) {
        case 1:
            # One column ==> symmetrical error bars
            call tbcfnd (tdp, lerrcol, secp, 1)
            if (secp <= 0) {
                call sprintf (errmsg, SZ_LINE,
                    "Cannot find errors column: %s\n")
                call pargstr (lerrcol)
                call error (0, errmsg)
            }
            call tbcgtr (tdp, secp, size, Memb[null], 1, numrows)
        case 2:
            # Two columns ==> asymmetrical errors (low high)
            call tbcfnd (tdp, lerrcol, lecp, 1)
            if (lecp <= 0) {
                call sprintf (errmsg, SZ_LINE,
                    "Cannot find low errors column: %s\n")
                call pargstr (lerrcol)
                call error (0, errmsg)
            }

            call tbcfnd (tdp, herrcol, hecp, 1)
            if (hecp <= 0) {
                call sprintf (errmsg, SZ_LINE,
                    "Cannot find high errors column: %s\n")
                call pargstr (herrcol)
                call error (0, errmsg)
            }

            # Total size of error bar
            call malloc (lerr, numrows, TY_REAL)
            call malloc (herr, numrows, TY_REAL)
            call tbcgtr (tdp, lecp, Memr[lerr], Memb[null], 1, numrows)
            call tbcgtr (tdp, hecp, Memr[herr], Memb[null], 1, numrows)
            call aaddr  (Memr[lerr], Memr[herr], size, numrows]

            # Position of centered marker
            call malloc (offs, numrows, TY_REAL)
            call asubr  (Memr[herr], Memr[lerr], Memr[offs], numrows)
            call adivkr (Memr[offs], 2.0, Memr[offs], numrows)

            if (istat == 1)
                # Errors in X;  adjust the horizontal position
                call aaddr (x, Memr[offs], x, numrows)
            else if (istat == 2)
                # Errors in Y;  adjust the vertical position
                call aaddr (y, Memr[offs], y, numrows)

            call mfree (lerr, TY_REAL)
            call mfree (herr, TY_REAL)
            call mfree (offs, TY_REAL)
        default:
            call amovkr(DEF_ERR, size, numrows)
        }
        call mfree (null, TY_BOOL)
end


# RC_RDXYCOL -- Read X and Y data from two scalar columns of the same table.
#
#
# Modified: 23-Sep-97 - I.Busko  -  To process row selectors.

int procedure rc_rdxycol (table, rowselect, xcolumn, ycolumn, rdmarks, 
                          errcol, istat, x, y, size)

char    table[SZ_FNAME]         # Table root name 
char    rowselect[SZ_FNAME]     # Row selector as parsed by rdselect
char    xcolumn[SZ_COLNAME]     # Column name
char    ycolumn[SZ_COLNAME]     # Column name
bool    rdmarks                 # Read errors from file?
char    errcol[SZ_COLNAME]      # Error column name
int     istat                   # X or Y errors?
pointer x, y, size              # Pointers to x, y and size vectors

pointer tdp                     # Pointer to table descriptor
pointer xcdp, ycdp, scdp        # Pointers to column descriptors
pointer	rcode			# Row selector structure
int     numrows                 # Maximum number of rows
int	nsrows			# Actual number of selected rows
int	irow, numcols, k
bool	null[1]

pointer tbtopn(), trsopen()
int     tbpsta(), selrows()
bool	trseval()

errchk	tbtopn, trsopen, malloc

begin
	numcols = 1
	tdp = tbtopn (table, READ_ONLY, 0)

	numrows = tbpsta (tdp, TBL_NROWS)

	# Find columns.
	call tbcfnd (tdp, xcolumn, xcdp, numcols)
	call tbcfnd (tdp, ycolumn, ycdp, numcols)
	if (xcdp <= 0)
	    call error (0, "Cannot find X column.")
	if (ycdp <= 0)
	    call error (0, "Cannot find Y column.")
	if (rdmarks && (errcol[1] != EOS)) {
	    call tbcfnd (tdp, errcol, scdp, numcols)
	    if (scdp <= 0)
	        call error (0, "Cannot find error column.")
	}

	# Compute number of selected rows.
	nsrows = selrows (tdp, rowselect)
	if (nsrows == 0)
	    call error (0, "No rows selected.")

	# Allocate space for the fit arrays.
	iferr {
	    call malloc (x,    nsrows, TY_REAL)
	    call malloc (y,    nsrows, TY_REAL)
	    call malloc (size, nsrows, TY_REAL)
	    call amovkr (1.0, Memr[size], nsrows)
	} then
	    call erract (EA_FATAL)

	# Scan rows and read from selected ones.
	rcode = trsopen (tdp, rowselect)
	k = 0
	do irow = 1, numrows {

	    if (trseval (tdp, irow, rcode)) {

	        # Read X and Y data.
	        call tbcgtr (tdp, xcdp, Memr[x+k], null, irow, irow)
	        call tbcgtr (tdp, ycdp, Memr[y+k], null, irow, irow)

	        # If selected, read error. This only supports a single,
	        # assumed symmetric, error bar. This code is executed
	        # only when the row selector selects a subsample of the
	        # full row range in the table.
	        if (rdmarks && (errcol[1] != EOS) && (nsrows != numrows))
	            call tbcgtr (tdp, scdp, Memr[size+k], null, irow, irow)

	        # Bump index.
	        k = k + 1
	    }
	}

	# This is old code kept in place for backwards compatibility.
	# It is capable of reading asymmetric error bars. This code
	# executes only when the row selector assumes the default value, 
	# that is, it selects all rows in the table.
	if (rdmarks && (errcol[1] != EOS) && (nsrows == numrows))
	    call rc_ggrser (tdp, errcol, numrows, Memr[x], Memr[y], 
                            Memr[size], istat)

	call trsclose (rcode)
	call tbtclo (tdp)
	return (nsrows)
end



# RC_RD2COL -- Read X and Y fit data from two scalar columns of 
#              separate tables.
#
#
# Modified: 23-Sep-97 - I.Busko  -  To process row selectors.

int procedure rc_rd2col (xtable, xcolumn, ytable, ycolumn, xrowselect,
                         yrowselect, rdmarks, errcol, istat, x, y, size)

char    xtable[SZ_FNAME]        # Table root names.
char	ytable[SZ_FNAME]
char    xcolumn[SZ_COLNAME]     # Column names.
char	ycolumn[SZ_COLNAME]
char    xrowselect[SZ_FNAME]	# Row selector in X table
char    yrowselect[SZ_FNAME]    # Row selector in Y table
bool    rdmarks                 # Read errors from file?
char    errcol[SZ_COLNAME]      # Error column name
int     istat                   # X or Y errors?
pointer x, y, size              # Pointers to x, y and size vectors

pointer xtdp, ytdp              # Pointers to table descriptors
pointer xcdp, ycdp, scdp        # Pointers to column descriptors
pointer	xrcode, yrcode		# Row selector structures
pointer	sp
pointer	xindex, yindex		# Arrays with row numbers.
int     xnumrows, ynumrows      # Maximum number of rows
int	nxsrows, nysrows	# Actual number of selected rows
int	irow, numcols, k
bool	null[1]

pointer tbtopn(), trsopen()
int     tbpsta(), selrows()
bool	trseval()

errchk	tbtopn, trsopen, malloc

begin
	numcols = 1
	xtdp = tbtopn (xtable, READ_ONLY, 0)
	ytdp = tbtopn (ytable, READ_ONLY, 0)

	xnumrows = tbpsta (xtdp, TBL_NROWS)
	ynumrows = tbpsta (ytdp, TBL_NROWS)

	# Find columns.
	call tbcfnd (xtdp, xcolumn, xcdp, numcols)
	call tbcfnd (ytdp, ycolumn, ycdp, numcols)
	if (xcdp <= 0)
	    call error (0, "Cannot find X column.")
	if (ycdp <= 0)
	    call error (0, "Cannot find Y column.")
	if (rdmarks && (errcol[1] != EOS)) {
	    # Errors are read from Y table.
	    call tbcfnd (ytdp, errcol, scdp, numcols)
	    if (scdp <= 0)
	        call error (0, "Cannot find error column.")
	}

	# Compute number of selected rows.
	nxsrows = selrows (xtdp, xrowselect)
	nysrows = selrows (ytdp, yrowselect)
	if (nxsrows == 0 || nysrows == 0)
	    call error (0, "No rows selected.")
	if (nxsrows != nysrows)
	    call error (0, "Selected different number of rows in each table.")

	# Allocate space for the fit arrays.
	iferr {
	    call malloc (x,    nxsrows, TY_REAL)
	    call malloc (y,    nxsrows, TY_REAL)
	    call malloc (size, nxsrows, TY_REAL)
	    call amovkr (1.0, Memr[size], nxsrows)
	} then
	    call erract (EA_FATAL)

	# Scan all rows in both tables and build arrays of row indices.
	# These are needed in order to one-to-one match selected rows 
	# from both tables.
	xrcode = trsopen (xtdp, xrowselect)
	yrcode = trsopen (ytdp, yrowselect)
	call smark (sp)
	call salloc (xindex, nxsrows, TY_INT)
	call salloc (yindex, nxsrows, TY_INT)
	k = 0
	do irow = 1, xnumrows {
	    if (trseval (xtdp, irow, xrcode)) {
	        Memi[xindex+k] = irow
	        k = k + 1
	    }
	}
	k = 0
	do irow = 1, ynumrows {
	    if (trseval (ytdp, irow, yrcode)) {
	        Memi[yindex+k] = irow
	        k = k + 1
	    }
	}
	call trsclose (xrcode)
	call trsclose (yrcode)

	# Read selected rows from both tables.
	do k = 0, nxsrows-1 {

	    # Read X and Y data.
	    call tbcgtr (xtdp, xcdp, Memr[x+k], null, Memi[xindex+k],
                                                      Memi[xindex+k])
	    call tbcgtr (ytdp, ycdp, Memr[y+k], null, Memi[yindex+k], 
                                                      Memi[yindex+k])

	    # If selected, read errors. This only supports a single,
	    # assumed symmetric, error bar. This code is executed
	    # only when the row selector selects a subsample of the
	    # full row range in the table.
	    if (rdmarks && (errcol[1] != EOS) && (nysrows != ynumrows))
	        call tbcgtr (ytdp, scdp, Memr[size+k], null, Memi[yindex+k], 
                                                             Memi[yindex+k])
	}

	# This is old code kept in place for backwards compatibility.
	# It is capable of reading asymmetric error bars. This code
	# executes only when the row selectors assume the default value, 
	# that is, they select all rows in the tables.
	if (rdmarks && (errcol[1] != EOS) && (nysrows == ynumrows))
	    call rc_ggrser (ytdp, errcol, ynumrows, Memr[x], Memr[y], 
                            Memr[size], istat)

	call tbtclo (xtdp)
	call tbtclo (ytdp)
	call sfree (sp)
	return (nxsrows)
end


# RC_PROJECTION -- Given an image section of arbitrary dimension, compute
# the projection along a single axis by taking the average over the other
# axes.  We do not know about bad pixels.

procedure rc_projection (im, pv, npix, axis)

pointer im                      # Pointer to image header
real    pv[npix]                # Receives the projection vector
int     npix                    # Length of projection vector
int     axis                    # The axis to be projected to (x=1)

int     i, lastv
long    v[IM_MAXDIM], nsum, totpix
pointer pix
real    asumr()
pointer imgnlr()
errchk  imgnlr

begin
        if (im == NULL)
            call error (1, "Image projection operator called with null im")
        if (axis < 1 || axis > IM_NDIM(im))
            call error (2, "Attempt to take projection over nonexistent axis")

        call aclrr (pv, npix)
        call amovkl (long(1), v, IM_MAXDIM)

        switch (axis) {
        case 1:
            # Since the image is read line by line, it is easy to compute the
            # projection along the x-axis (axis 1).  We merely sum all of the
            # image lines.
            while (imgnlr (im, pix, v) != EOF) 
                call aaddr (Memr[pix], pv, pv, npix)

        default:
            # Projecting along any other axis when reading the image line
            # by line is a bit difficult to understand.  Basically, the
            # element 'axis' of the V vector (position of the line in the
            # image) gives us the index into the appropriate element of
            # pv.  When computing the projection over multiple dimensions,
            # the same output element will be referenced repeatedly.  All
            # of the elmenents of the input line are summed and added into
            # this output element.

            for (lastv=v[axis];  imgnlr (im, pix, v) != EOF;  lastv=v[axis]) {
                i = lastv
                if (i <= npix)
                    pv[i] = pv[i] + asumr (Memr[pix], IM_LEN(im,1))
            }
        }

        # Now compute the number of pixels contributing to each element
        # of the output vector.  This is the number of pixels in the image
        # divided by the length of the projection.

        totpix = 1
        do i = 1, IM_NDIM(im)
            if (i == axis)
                totpix = totpix * min (npix, IM_LEN(im,i))
            else
                totpix = totpix * IM_LEN(im,i)
        nsum = totpix / min (npix, IM_LEN(im,axis))

        # Compute the average by dividing by the number if pixels summed at
        # each point.
        call adivkr (pv, real(nsum), pv, npix)
end


# GETAXIS  --  Returns the first non-singular axis number selected by
# a image section.

int procedure getaxis (section)

char    section[ARB]

int     i, j, len
int     commaindex[7], index

int     strlen()

begin
        index = 0
        len = strlen (section)
        if (len > 0) {
            j = 0
            do i = 1, len {
                if ((section[i] == ',') || (section[i] == ']')) {
                    j = j + 1
                    commaindex[j] = i
                }
            }
            do i = len, 1, -1 {
                if ((section[i] == '*') || (section[i] == ':'))
                    index = i
            }
            do i = 1, j {
                if (index < commaindex[i])
                    return (i)
            }
        }
        return (1)
end



# RC_RD1S -- Read fit data from a 1-D array in a row,column of a table.
#
# Installed: 16-Sep-97 - I.Busko

int procedure rc_rd1s (table, rowselect, column, rdmarks, errcol, istat, 
                       x, y, size)

char    table[SZ_FNAME]         # Table name
char    rowselect[SZ_FNAME]     # Row selector as parsed by rdselect
char    column[SZ_COLNAME]      # Data column name
bool    rdmarks                 # Read errors from file?
char    errcol[SZ_COLNAME]      # Error column name
int     istat                   # X or Y errors? (not used for now)
pointer x, y, size              # Pointers to x, y and size vectors

#--
char    colselect[SZ_FNAME]
char	newrow[SZ_FNAME], newtable[SZ_FNAME]
pointer	data[1], edata[1]
int	n, ncol, npix, npixe

errchk	rc_rdtb, malloc

begin
	# Cast column name in the form of a column selector.
	call sprintf (colselect, SZ_FNAME, "[c:%s]")
	    call pargstr (column)

	# Re-cast row selector in its original form between brackets,
	# since it was processed by rdselect. Skip if no row selector
	# is specified.
	if (rowselect[1] != EOS) {
	    call sprintf (newrow, SZ_FNAME, "[r:%s]")
	        call pargstr (rowselect)
	} else
	    newrow[1] = EOS

	# Build new table name with full selector construct in place.
	call strcpy (table, newtable, SZ_FNAME)
	call strcat (newrow, newtable, SZ_FNAME)
	call strcat (colselect, newtable, SZ_FNAME)

	# Read Y column.
	iferr {
	    call rc_rdtb (newtable, data, npix, ncol, 1)
	    y = data[1]
	} then
	    call erract (EA_WARN)
	if (ncol != 1)
	    call error (0, "Invalid Y column specification.")

	if (rdmarks && (errcol[1] != EOS)) {

	    # Read error array. Two-sided error bars are not supported yet.
	    call sprintf (colselect, SZ_FNAME, "[c:%s]")
	        call pargstr (errcol)
	    call strcpy (table, newtable, SZ_FNAME)
	    call strcat (newrow, newtable, SZ_FNAME)
	    call strcat (colselect, newtable, SZ_FNAME)
	    iferr {
	        call rc_rdtb (newtable, edata, npixe, ncol, 1)
	        size = edata[1]
	    } then
	        call erract (EA_WARN)
	    if (ncol != 1)
	        call error (0, "Invalid error column specification.")
	    if (npixe != npix)
	        call error (0, "Y and error arrays have different sizes.")
	} else {

	    # Build dummy error array.
	    iferr {
	        call malloc (size, npix, TY_REAL)
	    } then
	        call erract (EA_WARN)
	    do n = 1, npix
	        Memr[size+n-1] = 1.
	}

	# Generate X column (independent variable).
	iferr {
	    call malloc (x, npix, TY_REAL)
	} then
	    call erract (EA_WARN)
	do n = 1, npix
	    Memr[x+n-1] = n

	return (npix)
end



# RC_RD2S -- Read fit data from a 1-D array in row,column of 2 tables.
#
# Installed: 16-Sep-97 - I.Busko

int procedure rc_rd2s (xtable, xcolumn, ytable, ycolumn, 
                       rowselect, rowselect2, rdmarks, errcol, 
                       istat, x, y, size)

char    xtable[SZ_FNAME]        # X table name
char    ytable[SZ_FNAME]        # Y table name
char    xcolumn[SZ_COLNAME]     # X data column name
char    ycolumn[SZ_COLNAME]     # Y data column name
char    rowselect[SZ_FNAME]     # X table row selector as parsed by rdselect
char    rowselect2[SZ_FNAME]    # Y table row selector as parsed by rdselect
bool    rdmarks                 # Read errors from file?
char    errcol[SZ_COLNAME]      # Error column name
int     istat                   # X or Y errors? (not used for now)
pointer x, y, size              # Pointers to x, y and size vectors

#--
char    colselect[SZ_FNAME]
char	newrow[SZ_FNAME], newtable[SZ_FNAME]
pointer	data[1], edata[1]
int	n, ncol, npixx, npixy, npixe

errchk	rc_rdtb, malloc

begin
	# Cast column name in the form of a column selector.
	call sprintf (colselect, SZ_FNAME, "[c:%s]")
	    call pargstr (xcolumn)

	# Re-cast row selector in its original form between brackets,
	# since it was processed by rdselect. Skip if no row selector
	# is specified.
	if (rowselect[1] != EOS) {
	    call sprintf (newrow, SZ_FNAME, "[r:%s]")
	        call pargstr (rowselect)
	} else
	    newrow[1] = EOS

	# Build new table name with full selector construct in place.
	call strcpy (xtable, newtable, SZ_FNAME)
	call strcat (newrow, newtable, SZ_FNAME)
	call strcat (colselect, newtable, SZ_FNAME)

	# Read X column.
	iferr {
	    call rc_rdtb (newtable, data, npixx, ncol, 1)
	    x = data[1]
	} then
	    call erract (EA_WARN)
	if (ncol != 1)
	    call error (0, "Invalid X column specification.")

	# Now repeat the above for the Y data.
	call sprintf (colselect, SZ_FNAME, "[c:%s]")
	    call pargstr (ycolumn)
	call sprintf (newrow, SZ_FNAME, "[r:%s]")
	    call pargstr (rowselect2)
	call strcpy (ytable, newtable, SZ_FNAME)
	call strcat (newrow, newtable, SZ_FNAME)
	call strcat (colselect, newtable, SZ_FNAME)
	iferr {
	    call rc_rdtb (newtable, data, npixy, ncol, 1)
	    y = data[1]
	} then
	    call erract (EA_WARN)
	if (ncol != 1)
	    call error (0, "Invalid X column specification.")
	if (npixx != npixy)
	    call error (0, "X and Y arrays have different sizes.")

	if (rdmarks && (errcol[1] != EOS)) {

	    # Read error array. Two-sided error bars are not supported yet.
	    # The error is assumed to be the Y errror, so it is read from
	    # the Y table, using the Y row selector.
	    call sprintf (colselect, SZ_FNAME, "[c:%s]")
	        call pargstr (errcol)
	    call strcpy (ytable, newtable, SZ_FNAME)
	    call strcat (newrow, newtable, SZ_FNAME)
	    call strcat (colselect, newtable, SZ_FNAME)
	    iferr {
	        call rc_rdtb (newtable, edata, npixe, ncol, 1)
	        size = edata[1]
	    } then
	        call erract (EA_WARN)
	    if (ncol != 1)
	        call error (0, "Invalid error column specification.")
	    if (npixe != npixx)
	        call error (0, "Y and error arrays have different sizes.")
	} else {

	    # Build dummy error array.
	    iferr {
	        call malloc (size, npixx, TY_REAL)
	    } then
	        call erract (EA_WARN)
	    do n = 1, npixx
	        Memr[size+n-1] = 1.
	}

	return (npixx)
end



# RC_RDXYS -- Read fit data from a 1-D array in a row and two columns of 
#             a table.
#
# Installed: 16-Sep-97 - I.Busko

int procedure rc_rdxys (table, rowselect, xcolumn, ycolumn, rdmarks, 
                        errcol, istat, x, y, size)

char    table[SZ_FNAME]         # Table name
char    rowselect[SZ_FNAME]     # Row selector as parsed by rdselect
char    xcolumn[SZ_COLNAME]     # X column name
char    ycolumn[SZ_COLNAME]     # Y column name
bool    rdmarks                 # Read errors from file?
char    errcol[SZ_COLNAME]      # Error column name
int     istat                   # X or Y errors? (not used for now)
pointer x, y, size              # Pointers to x, y and size vectors

#--
char    colselect[SZ_FNAME]
char	newrow[SZ_FNAME], newtable[SZ_FNAME]
pointer	data[2], edata[1]
int	n, ncol, npix, npixe

errchk	rc_rdtb, malloc

begin
	# Cast column names in the form of a column selector.
	call sprintf (colselect, SZ_FNAME, "[c:%s,%s]")
	    call pargstr (xcolumn)
	    call pargstr (ycolumn)

	# Re-cast row selector in its original form between brackets,
	# since it was processed by rdselect. Skip if no row selector
	# is specified.
	if (rowselect[1] != EOS) {
	    call sprintf (newrow, SZ_FNAME, "[r:%s]")
	        call pargstr (rowselect)
	} else
	    newrow[1] = EOS

	# Build new table name with full selector construct in place.
	call strcpy (table, newtable, SZ_FNAME)
	call strcat (newrow, newtable, SZ_FNAME)
	call strcat (colselect, newtable, SZ_FNAME)

	# Read X and Y column.
	iferr {
	    call rc_rdtb (newtable, data, npix, ncol, 2)
	    x = data[1]
	    y = data[2]
	} then
	    call erract (EA_WARN)
	if (ncol != 2)
	    call error (0, "Invalid X,Y column specification.")

	if (rdmarks && (errcol[1] != EOS)) {

	    # Read error array. Two-sided error bars are not supported yet.
	    call sprintf (colselect, SZ_FNAME, "[c:%s]")
	        call pargstr (errcol)
	    call strcpy (table, newtable, SZ_FNAME)
	    call strcat (newrow, newtable, SZ_FNAME)
	    call strcat (colselect, newtable, SZ_FNAME)
	    iferr {
	        call rc_rdtb (newtable, edata, npixe, ncol, 1)
	        size = edata[1]
	    } then
	        call erract (EA_WARN)
	    if (ncol != 1)
	        call error (0, "Invalid error column specification.")
	    if (npixe != npix)
	        call error (0, "X,Y and error arrays have different sizes.")
	} else {

	    # Build dummy error array.
	    iferr {
	        call malloc (size, npix, TY_REAL)
	    } then
	        call erract (EA_WARN)
	    do n = 1, npix
	        Memr[size+n-1] = 1.
	}

	return (npix)
end



# RC_RDTB --  Read double or float arrays from columns and multiple rows.
#
# This routine replaces 'omnread' from the 'selectors' library. The important
# difference is that it can read data from more than one row in a single
# call, generating an output 1-D array which stores data from all selected
# rows. It does not pay attention to possible multi-dimensional arrays
# stored in the table; everything is handled as a 1-D vector. Only float 
# output is supported for now.
#
#
# Installed: 19-Sep-97 - I.Busko

procedure rc_rdtb (file, data, nelem, ncol, maxcol)

char	file[ARB]	# i: file name, including sections or selectors
pointer	data[ARB]	# o: pointers to columns of output data
int	nelem		# o: length of output columns
int	ncol		# o: number of output columns
int	maxcol		# i: maximum number of columns
#--
pointer	root, rowselect, colselect, sp, tp, col, rcode, length
int	irow, nrow, nsrow, icol, nscalar, size, blk

pointer	tbtopn(), trsopen()
int	selrows(), tcs_totsize(), tbpsta()
bool	trseval()

errchk	rdselect, tbtopn, tcs_open, trsopen, om_error

begin
	# Allocate temporary arrays.
	call smark (sp)
	call salloc (root,      SZ_FNAME, TY_CHAR)
	call salloc (rowselect, SZ_FNAME, TY_CHAR)
	call salloc (colselect, SZ_FNAME, TY_CHAR)
	call salloc (col,       maxcol,   TY_INT)
	call salloc (length,    maxcol,   TY_INT)

	# Break the table name into its various parts.
	call rdselect (file, Memc[root], Memc[rowselect], Memc[colselect], 
                       SZ_FNAME)

	# Open table.
	tp = tbtopn (Memc[root], READ_ONLY, NULL)

	# Get total array length of each selected column.
	# It is an error to select scalar columns.
	call tcs_open (tp, Memc[colselect], Memi[col], ncol, maxcol)
	if (ncol == 0)
	    call error (0, "No columns selected.")
	if (ncol != maxcol)
	    call error (0, "Error in column selector.")
	nscalar = 0
	do icol = 1, ncol {
	    Memi[length+icol-1] = tcs_totsize (Memi[col+icol-1])
	    if (Memi[length+icol-1] == 1)
		nscalar = nscalar + 1
	}
	if (nscalar > 0)
	    call error (0, "Scalar column selected.")

	# If more than one column selected, check array lengths 
	# to make sure they are equal.
	if (ncol > 1) {
	    do icol = 1, ncol-1 {
	        if (Memi[length+icol-1] != Memi[length+icol])
	            call error (0, 
                    "Selected columns have different array sizes.")
	    }
	}
	nelem = Memi[length]

	# Compute number of selected rows.
	nsrow = selrows (tp, Memc[rowselect])
	if (nsrow == 0)
	    call error (0, "No rows selected.")

	# Alloc output arrays. They must be sized to fit
	# all selected rows.
	do icol = 1, ncol
	    call malloc (data[icol], nelem * nsrow, TY_REAL)

	# Scan rows selected by row selector.
	rcode = trsopen (tp, Memc[rowselect])
	nrow = tbpsta (tp, TBL_NROWS)
	blk = 0
	do irow = 1, nrow {
	    if (trseval (tp, irow, rcode)) {

	        # If row matched, scan columns and read arrays. Store
	        # at appropriate offset in output buffer.
		do icol = 1, ncol {
	            call tcs_rdaryr (tp, Memi[col+icol-1], irow, nelem, size, 
                                     Memr[data[icol] + nelem * blk])
	        }

	        # Bump pointer in output arrays.
	        blk = blk + 1
	    }
	}

	# Total number of output pixels.
	nelem = nelem * blk

	call trsclose (rcode)
	call tbtclo (tp)
	call sfree (sp)
end




# RC_CDIM -- Get column dimensionality.
#
#
# Installed: 23-Sep-97 - I.Busko

int procedure rc_cdim (table, column)

char	table[ARB]	# i: root table name
char	column[ARB]	# i: column name 
#--
pointer	tp
int	col, ncol, length

pointer	tbtopn()
int	tcs_totsize()

errchk	rdselect, tbtopn, tcs_open, trsopen, 

begin
	# Open table.
	tp = tbtopn (table, READ_ONLY, NULL)

	# Select column.
	call tcs_open (tp, column, col, ncol, 1)
	if (ncol == 0)
	    call error (0, "No columns selected.")

	# Get array length.
	length = tcs_totsize (col)

	call tbtclo (tp)
	return (length)
end
