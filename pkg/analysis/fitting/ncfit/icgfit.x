include	<error.h>
include	<pkg/gtools.h>
include	"icfit.h"

# ICG_FIT -- Interactive curve fitting with graphics.  
# This is the main entry point for the interactive graphics part of the 
# icfit package. Returns ERR if 'a' key is returned by cursor.

int procedure icg_fit (ic, gp, cursor, gt, nl, x, y, wts, npts)

pointer	ic			# ICFIT pointer
pointer	gp			# GIO pointer
char	cursor[ARB]		# GIO cursor input
pointer	gt			# GTOOLS pointer
pointer	nl			# NLFIT pointer
real	x[npts]			# Ordinates
real	y[npts]			# Abscissas
real	wts[npts]		# Weights
int	npts			# Number of points
#--

real	wx, wy
int	wcs, key
char	cmd[SZ_LINE]
int	i, newgraph, newaxis, axes[2], xtype
real	x1, aux
real	rx1, rx2, ry1, ry2
pointer	sp, userwts

int	icg_cursor(), stridxs(), scan(), nscan()
int	icg_nearest(), nl_stati()
real	nl_zeval()

errchk	ic_fit, icg_graph

begin
	# Allocate memory for the fit and a copy of the weights.
	# The weights are copied because they are changed when points are
	# deleted.

	call smark (sp)
	call salloc (userwts, npts, TY_REAL)
	call amovr (wts, Memr[userwts], npts)

	# Initialize
	IC_OVERPLOT(ic) = NO
	IC_NEWX(ic) = YES
	IC_NEWY(ic) = YES
	IC_NEWWTS(ic) = YES
	IC_NEWFUNCTION(ic) = YES

	# Read cursor commands.

	key = '^'			# Don't fit, just plot it.
	newgraph = YES
	IC_FITERROR(ic) = NO
	axes[1] = IC_AXES(ic, IC_GKEY(ic), 1)
	axes[2] = IC_AXES(ic, IC_GKEY(ic), 2)
	xtype = 0

	repeat {
	    switch (key) {

	    case '^': # Draw graph for the first time
	        call ic_dosetup (ic, nl, x, wts, npts, YES, YES, YES)
		newgraph = YES

	    case '?': # Print help text.
		call gpagefile (gp, Memc[IC_HELP(ic)], IC_PROMPT)

	    case ':': # List or set parameters
		if (cmd[1] == '/')
	            call gt_colon (cmd, gp, gt, newgraph)
		else
		    call icg_colon (ic, cmd, newgraph, newaxis, gp,gt,nl, 
				    x,y,wts, npts)

	    case 'a': # Returns with error code
	        call sfree (sp)
	        return (ERR)

	    case 'b': # Set baseline for Gaussians
		rx1 = wx
	        rx2 = wx
		ry1 = wy
	        while (rx2 == rx1) {
	            call printf (" again:\n"); call flush(STDOUT)
		    i = icg_cursor(ic, gt, cursor, rx2, ry2, wcs, key, cmd, SZ_LINE)
	        }
	        call icg_fi2 (nl, rx1, ry1, rx2, ry2)
	        call printf ("\n"); call flush(STDOUT)
		newgraph = YES

	    case 'c': # Print data point nearest to cursor.
		i = icg_nearest (ic, gp, gt, nl, x, y, npts, wx, wy)

	    	if (i != 0) {
		    call printf ("x = %g  y = %g  fit = %g\n")
			call pargr (x[i])
	                switch (IC_YAXIS(ic)) {
	                case IC_LINEAR:
			    call pargr (y[i])
			    call pargr (nl_zeval (nl, x[i], 0.))
	                case IC_LOG:
	                    if (y[i] > 0.)
	                        aux = log10 (y[i])
	                    else
	                        aux = 0.
			    call pargr (aux)
	                    aux = nl_zeval (nl, x[i], 0.)
	                    if (aux > 0.)
	                        aux = log10 (aux)
	                    else
	                        aux = 0.
			    call pargr (aux)
	                case IC_MAG:
	                    if (y[i] > 0.)
	                        aux = IC_MAG0(ic) - 2.5 * log10 (y[i])
	                    else
	                        aux = 0.
			    call pargr (aux)
	                    aux = nl_zeval (nl, x[i], 0.)
	                    if (aux > 0.)
	                        aux = IC_MAG0(ic) - 2.5 * log10 (aux)
	                    else
	                        aux = 0.
			    call pargr (aux)
	                default:
	                }
		    call flush(STDOUT)
	        }

	    case 'd': # Delete data points.
		call icg_delete (ic, gp, gt, nl, x, y, wts, Memr[userwts],
		    npts, wx, wy)

	    case 'f': # Fit the function and reset the flags.
		iferr {
		    call ic_fit (ic, nl, x, y, wts, npts, IC_NEWX(ic),
			IC_NEWY(ic), IC_NEWWTS(ic), IC_NEWFUNCTION(ic))

		    IC_NEWX(ic) = NO
		    IC_NEWY(ic) = NO
		    IC_NEWWTS(ic) = NO
		    IC_NEWFUNCTION(ic) = NO
		    IC_FITERROR(ic) = NO
		    newgraph = YES
		} then {
		    IC_FITERROR(ic) = YES
		    call erract (EA_WARN)
		}

	    case 'g':	# Set graph axes types.
		call printf ("Graph key to be defined: ")
		call flush (STDOUT)
		if (scan() == EOF)
		    goto 10
		call gargc (cmd[1])

		switch (cmd[1]) {
		case '\n':
		case 'h', 'i', 'j', 'k', 'l':
		    switch (cmd[1]) {
		    case 'h':
		        key = 1
		    case 'i':
		        key = 2
		    case 'j':
		        key = 3
		    case 'k':
		        key = 4
		    case 'l':
			key = 5
		    }

		    call printf ("Set graph axes types (%c, %c): ")
		        call pargi (IC_AXES(ic, key, 1))
		        call pargi (IC_AXES(ic, key, 2))
		    call flush (STDOUT)
		    if (scan() == EOF)
		        goto 10
		    call gargc (cmd[1])

		    switch (cmd[1]) {
		    case '\n':
		    default:
		        call gargc (cmd[2])
		        call gargc (cmd[2])
		        if (cmd[2] != '\n') {
			    IC_AXES(ic, key, 1) = cmd[1]
			    IC_AXES(ic, key, 2) = cmd[2]
			    if (IC_GKEY(ic) == key)
				newgraph = YES
		        }
		    }
		default:
		    call printf ("Not a graph key\n"); call flush(STDOUT)
		}

	    case 'h':
		if (IC_GKEY(ic) != 1) {
		    IC_GKEY(ic) = 1
		    newgraph = YES
		}

	    case 'i':
		if (IC_GKEY(ic) != 2) {
		    IC_GKEY(ic) = 2
		    newgraph = YES
		}

	    case 'j':
		if (IC_GKEY(ic) != 3) {
		    IC_GKEY(ic) = 3
		    newgraph = YES
		}

	    case 'k':
		if (IC_GKEY(ic) != 4) {
		    IC_GKEY(ic) = 4
		    newgraph = YES
		}

	    case 'l':
		if (IC_GKEY(ic) != 5) {
		    IC_GKEY(ic) = 5
		    newgraph = YES
		}

	    case 't': # Initialize the sample string and erase from the graph.
		call icg_sample (ic, gp, gt, x, npts, 0)
		call sprintf (Memc[IC_SAMPLE(ic)], SZ_LINE, "*")
		IC_NEWX(ic) = YES

	    case 'o': # Set overplot flag
		IC_OVERPLOT(ic) = YES

	    case 'p': # Set Gaussians amplitude and center
	        if ( (nl_stati( nl, "fitfunc") == GAUSSIANS) ||
		     (nl_stati( nl, "fitfunc") == CGAUSS)) {
	            i = 1
	            call icg_fi1 (nl, wx, wy, i)
	            call printf ("\n Enter more Gaussians. Press 'q' to quit. ")
		    call flush(STDOUT)
		    while( icg_cursor(ic, gt, cursor, wx, wy, wcs, key, cmd, 
			   SZ_LINE) != EOF) {
	                if (key == 'q')
	                    break
	                i = i + 1
	                call icg_fi1 (nl, wx, wy, i)
	            }
	            newgraph = YES
	        }

	    case 'r': # Redraw the graph
		newgraph = YES

	    case 's': # Set sample regions with the cursor.
		if ((IC_AXES(ic,IC_GKEY(ic),1) == 'x') ||
		    (IC_AXES(ic,IC_GKEY(ic),2) == 'x')) {
		    if (stridxs ("*", Memc[IC_SAMPLE(ic)]) > 0)
		        Memc[IC_SAMPLE(ic)] = EOS

		    rx1 = wx
		    ry1 = wy
		    call printf (" again:\n"); call flush(STDOUT)
		    if (icg_cursor(ic, gt, cursor, wx, wy, wcs, key, cmd, 
			SZ_LINE) == EOF)
		        break
		    call printf ("\n"); call flush(STDOUT)
		    rx2 = wx
		    ry2 = wy

		    # Determine if the x vector is integer.
		    if (xtype == 0) {
			xtype = TY_INT
			do i = 1, npts
			    if (x[i] != int (x[i])) {
				xtype = TY_REAL
				break
			    }
		    }

		    if (IC_AXES(ic,IC_GKEY(ic),1) == 'x') {
		        if (xtype == TY_INT) {
		            call sprintf (cmd, SZ_LINE, " %d:%d")
		                call pargi (nint (rx1))
		                call pargi (nint (rx2))
			} else {
		            call sprintf (cmd, SZ_LINE, " %g:%g")
		                call pargr (rx1)
		                call pargr (rx2)
			}
		    } else {
		        if (xtype == TY_INT) {
		            call sprintf (cmd, SZ_LINE, " %d:%d")
		                call pargi (nint (ry1))
		                call pargi (nint (ry2))
			} else {
		            call sprintf (cmd, SZ_LINE, " %g:%g")
		                call pargr (ry1)
		                call pargr (ry2)
			}
		    }
		    call strcat (cmd, Memc[IC_SAMPLE(ic)], SZ_LINE)
		    call icg_sample (ic, gp, gt, x, npts, 1)
		    IC_NEWX(ic) = YES
		}

	    case 'u': # Undelete data points.
		call icg_undelete (ic, gp, gt, nl, x, y, wts, Memr[userwts],
		    npts, wx, wy)

	    case 'w':  # Window graph
		call gt_window (gt, gp, cursor, newgraph)

	    case 'x': # Reset the value of the x point.
		i = icg_nearest (ic, gp, gt, nl, x, y, npts, wx, wy)

	    	if (i != 0) {
		    call printf ("x = (%g) ")
			call pargr (x[i])
		    call flush (STDOUT)
		    if (scan() != EOF) {
		        call gargr (x1)
		        if (nscan() == 1) {
			    if (!IS_INDEF (x1)) {
			        x[i] = x1
			        IC_NEWX(ic) = YES
			    }
			}
		    }
		}

	    case 'y': # Reset the value of the y point.
		i = icg_nearest (ic, gp, gt, nl, x, y, npts, wx, wy)

	    	if (i != 0) {
		    call printf ("y = (%g) ")
			call pargr (y[i])
	            if (IC_YAXIS(ic) != LINEAR)
	                call printf ("  Beware of the y scale !")
		    call flush (STDOUT)
		    if (scan() != EOF) {
		        call gargr (x1)
		        if (nscan() == 1) {
			    if (!IS_INDEF (x1)) {
			        y[i] = x1
			        IC_NEWY(ic) = YES
			    }
			}
		    }
		}

	    case 'z': # Print cursor coordinates.
		call printf ("x = %g  y = %g \n")
		    call pargr (wx)
	            call pargr (wy)
	        call flush (STDOUT)

	    case 'I': # Interrupt
		call fatal (0, "Interrupt")

	    default: # Let the user decide on any other keys.
		call icg_user (ic, gp, gt, nl, wx, wy, wcs, key, cmd)
	    }

	    # Redraw the graph if necessary.
10	    if (newgraph == YES) {
		if ((IC_AXES(ic, IC_GKEY(ic), 1) != axes[1]) ||
		    (newaxis == YES)) {
		    axes[1] = IC_AXES(ic, IC_GKEY(ic), 1)
		    call gt_setr (gt, GTXMIN, INDEFR)
		    call gt_setr (gt, GTXMAX, INDEFR)
		}
		if ((IC_AXES(ic, IC_GKEY(ic), 2) != axes[2]) ||
		    (newaxis == YES)) {
		    axes[2] = IC_AXES(ic, IC_GKEY(ic), 2)
		    call gt_setr (gt, GTYMIN, INDEFR)
		    call gt_setr (gt, GTYMAX, INDEFR)
		}

	    	call icg_graph (ic, gp, gt, nl, x, y, wts, npts)
		newgraph = NO
	        newaxis  = NO
	    }
	} until (icg_cursor (ic, gt, cursor, wx, wy, wcs, key, cmd, 
		 SZ_LINE) == EOF)

	call sfree (sp)
	return (OK)
end


# ICG_FI1 -- Add one Gaussian to nlfit control structure.

procedure icg_fi1 (nl, wx, wy, n)

pointer	nl		# NLFIT pointer
real	wx, wy		# cursor position
int	n		# Gaussian number

pointer	coeff, cflags, sp
int	fitfunc, npts, npar, i
real	baseline

int	nl_stati()

begin
	# Get old coefficients.
	npar = 3 * n + 2
	npts    = nl_stati (nl, "npts")
	fitfunc = nl_stati (nl, "fitfunc")
	call smark (sp)
	call salloc (coeff,  npar, TY_REAL)
	call salloc (cflags, npar, TY_BOOL)
	call nl_gcoeff (nl, Memr[coeff], i)
	call nl_gflags (nl, Memb[cflags],i)

	# Compute baseline at given position
	baseline = Memr[coeff + NL_GA] + Memr[coeff + NL_GB] * wx

	# Store new Gaussian
	Memr[coeff  + NL_GAMPL(n)] = wy - baseline
	Memr[coeff  + NL_GCENT(n)] = wx
	Memr[coeff  + NL_GFWHM(n)] = Memr[coeff + NL_GFWHM(1)]
	Memb[cflags + NL_GAMPL(n)] = true
	Memb[cflags + NL_GCENT(n)] = true
	Memb[cflags + NL_GFWHM(n)] = false
	if ((fitfunc == CGAUSS) && (n != 1)) {
	    Memr[coeff + NL_GAMPL(n)] = Memr[coeff + NL_GAMPL(n)] / Memr[coeff + NL_GAMPL(1)]
	    Memr[coeff + NL_GCENT(n)] = Memr[coeff + NL_GCENT(n)] - Memr[coeff + NL_GCENT(1)]
	}

	call nl_init (nl, fitfunc, Memr[coeff], Memb[cflags], npar, npts)
	call sfree (sp)
end


# ICG_FI2 -- Modify baseline.

procedure icg_fi2 (nl, wx1, wy1, wx2, wy2)

pointer	nl		# NLFIT pointer
real	wx1, wy1	# cursor position 1
real	wx2, wy2	# cursor position 2

pointer	coeff, cflags, sp
int	fitfunc, npts, npar, i
real	zero, slope

int	nl_stati()

begin
	# Get old coefficients.
	fitfunc = nl_stati (nl, "fitfunc")
	npar    = nl_stati (nl, "npar")
	npts    = nl_stati (nl, "npts")
	call smark (sp)
	call salloc (coeff,  npar, TY_REAL)
	call salloc (cflags, npar, TY_BOOL)
	call nl_gcoeff (nl, Memr[coeff], i)
	call nl_gflags (nl, Memb[cflags],i)

	# Compute baseline coefficients
	slope = (wy2 - wy1) / (wx2 - wx1)
	zero  = ((wy1 - slope * wx1) + (wy2 - slope * wx2)) / 2.
	Memr[coeff + NL_GA] = zero
	Memr[coeff + NL_GB] = slope

	call nl_init (nl, fitfunc, Memr[coeff], Memb[cflags], npar, npts)
	call sfree (sp)
end

