include	"icfit.h"

# IC_VSHOW -- Show fit parameters in verbose mode.

procedure ic_vshow (ic, file, nl, x, y, wts, npts, gt)

pointer	ic		# ICFIT pointer
char	file[ARB]	# Output file
pointer	nl		# NLFIT pointer
real	x[ARB]		# Ordinates
real	y[ARB]		# Abscissas
real	wts[ARB]	# Weights
int	npts		# Number of data points
pointer	gt		# Graphics tools pointer
#--

int	i, n, deleted, fd

int	open()

errchk	open()

begin
	# Do the standard ic_show option, then add on the verbose part.
	call ic_show (ic, file, gt, nl)

	if (npts == 0) {
	    call eprintf ("Incomplete output - no data points for fit\n")
	    return
	}

	# Open the output file.
	fd = open (file, APPEND, TEXT_FILE)

	# Count rejected and deleted points
	if (npts == IC_NFIT(ic)) 
	    n = npts
	else
	    n = IC_NFIT(ic)

	deleted = 0
	do i = 1, n {
	    if (wts[i] == 0.)
	        deleted = deleted + 1
	}

	# Print.
	call fprintf (fd, "total points  = %d\nsample points = %d\n")
	    call pargi (npts)
	    call pargi (n)
	call fprintf (fd, "rejected      = %d\ndeleted       = %d\n")
	    call pargi (IC_NREJECT(ic))
	    call pargi (deleted)

	# Print x,y pairs
	call ic_listxyr (fd, nl, x, y, npts)

	call close (fd)
end


# IC_LISTXY -- List data as x,y pairs on output.  Used for verbose
# show procedure.

procedure ic_listxyr (fd, nl, xvals, yvals, nvalues)

int	fd			# File descriptor of output file
pointer	nl			# Pointer to NLFIT structure
int	nvalues			# Number of data values
real	xvals[nvalues]		# Array of x data values
real	yvals[nvalues]		# Array of y data values

int	i
real	y, nl_zeval()

begin
	call fprintf (fd, "\n\t         X     \t   Yc   \t    Y    \t   Err \n")

	do i = 1, nvalues {
	    y = nl_zeval (nl, xvals[i], 0.)
	    call fprintf (fd, "\t%14.7e \t%14.7e \t%14.7e \t%14.7e\n")
		call pargr (xvals[i])
		call pargr (y)
		call pargr (yvals[i])
		call pargr (yvals[i] - y)
	}
end
