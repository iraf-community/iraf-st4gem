include	<gset.h>
include	<pkg/gtools.h>
include	<mach.h>
include	"icfit.h"

define	NGRAPH		100		# Number of fit points to graph
define	SIZE_BOX	2.		# Box size
define	SIZE_REJECT	2.		# Mark size of rejected points.

# ICG_GRAPH -- Graph data and fit.

procedure icg_graph (ic, gp, gt, nl, x, y, wts, npts)

pointer	ic				# ICFIT pointer
pointer	gp				# GIO pointer
pointer	gt				# GTOOLS pointers
pointer	nl				# NLFIT pointer
real	x[npts]				# Independent variable
real	y[npts]				# Dependent variable
real	wts[npts]			# Weights
int	npts				# Number of points
#--

pointer sp, ytemp, ubars, lbars
pointer	xout, yout
int	i
real	size, aux1, aux2

int	nl_stati()

errchk	icg_axes, icg_param, icg_g1, icg_g2, icg_gf

begin
	# Prepare y data for plotting
	call smark (sp)
	call salloc (ytemp, npts, TY_REAL)
	call salloc (ubars, npts, TY_REAL)
	call salloc (lbars, npts, TY_REAL)
	call amovr (y, Memr[ytemp], npts)
	IC_MINY(ic) = MAX_REAL
	do i = 1, npts {
	    # Find minimum positive y.
	    if ((y[i] > 0.) && (y[i] < IC_MINY(ic)))
	        IC_MINY(ic) = y[i]
	    # Set error bar of each datum.
	    if (wts[i] > 0.) {
	        Memr[ubars+i-1] = 1.0/ wts[i]
	        Memr[lbars+i-1] = 1.0/ wts[i]
	    } else {
	        Memr[ubars+i-1] = 0.
	        Memr[lbars+i-1] = 0.
	    }
	}

	# Set axis type.
	call gseti (gp, G_YTRAN, GW_LINEAR)
	switch (IC_XAXIS(ic)) {
	case IC_LOG:
	    call gt_sets (gt, GTXTRAN, "logrithmic")
	default:
	    call gt_sets (gt, GTXTRAN, "linear")
	}
	switch (IC_YAXIS(ic)) {
	case IC_LOG:
	    call gt_sets (gt, GTYTRAN, "logrithmic")
	case IC_MAG:
	    call gt_sets (gt, GTYTRAN, "linear")
	    do i = 0, npts-1 {
	        if (Memr[ytemp+i] < IC_MINY(ic)) {
	            Memr[ytemp+i] = IC_MAG0(ic) - 2.5 * log10 (IC_MINY(ic))
	            Memr[ubars+i] = 0.
	            Memr[lbars+i] = 0.
	        } else {
	            aux1 = Memr[ytemp+i] + Memr[ubars+i]
	            aux2 = Memr[ytemp+i] - Memr[lbars+i]
	            aux1 = IC_MAG0(ic) - 2.5 * log10 (aux1)
	            if (aux2 > 0.) 
	                aux2 = IC_MAG0(ic) - 2.5 * log10 (aux2)
	            else
	                aux2 = INDEF
	            Memr[ytemp+i] = IC_MAG0(ic) - 2.5 * log10 (Memr[ytemp+i])
	            Memr[ubars+i]  = Memr[ytemp+i] - aux1
	            if (!IS_INDEFR (aux2))
	                Memr[lbars+i]  = aux2 - Memr[ytemp+i]
	            else
	                Memr[lbars+i]  = 0.
	        }
	    }
	default:
	    call gt_sets (gt, GTYTRAN, "linear")
	}
	call malloc (xout, npts, TY_REAL)
	call malloc (yout, npts, TY_REAL)
	call icg_axes (ic, gt, nl, 1, x, Memr[ytemp], Memr[xout], npts)
	call icg_axes (ic, gt, nl, 2, x, Memr[ytemp], Memr[yout], npts)
	call icg_params (ic, nl, x, Memr[ytemp], wts, npts, gt)
	call icg_g1 (ic, gp, gt, nl, Memr[xout], Memr[yout], wts,
		     Memr[ubars], Memr[lbars], npts)

	# Symbol size for averaged ranges
	size = abs (IC_NAVERAGE(ic) * (Memr[xout+npts-1] - Memr[xout]) / 
	       real (npts))

	if (npts != IC_NFIT(ic)) {
	    if ((abs (IC_NAVERAGE(ic)) > 1) || (IC_NREJECT(ic) > 0)) {
	        call realloc (xout, IC_NFIT(ic), TY_REAL)
		call realloc (yout, IC_NFIT(ic), TY_REAL)
		call icg_axes (ic, gt, nl, 1, Memr[IC_XFIT(ic)],
		    Memr[IC_YFIT(ic)], Memr[xout], IC_NFIT(ic))
		call icg_axes (ic, gt, nl, 2, Memr[IC_XFIT(ic)],
		    Memr[IC_YFIT(ic)], Memr[yout], IC_NFIT(ic))
		call icg_g2 (ic, gp, gt, Memr[xout], Memr[yout],
		    IC_NFIT(ic), size)
	    }

	} else if (IC_NREJECT(ic) > 0)
	    call icg_g2 (ic, gp, gt, Memr[xout], Memr[yout], npts, size)

	i = nl_stati (nl, "npts")
	call icg_gf (ic, gp, gt, nl, max (npts, NGRAPH))
	call nl_puti (nl, "npts", i)

	# Mark the the sample regions.

	call icg_sample (ic, gp, gt, x, npts, 1)

	call mfree (xout, TY_REAL)
	call mfree (yout, TY_REAL)
	call sfree (sp)
end


# Graph data points.

procedure icg_g1 (ic, gp, gt, nl, x, y, wts, ubars, lbars, npts)

pointer	ic				# ICFIT pointer
pointer	gp				# GIO pointer
pointer	gt				# GTOOLS pointer
pointer	nl				# NLFIT pointer
real	x[npts]				# Ordinates
real	y[npts]				# Abscissas
real	wts[npts]			# Weights
real	ubars[npts]			# Upper error bars
real	lbars[npts]			# Lower error bars
int	npts				# Number of points

pointer	sp, xr, yr, xr1, yr1, gt1, gt2
int	i
real	aux, x1, x2, y1, y2
bool	wflag

int	nl_stati()

begin
	call smark (sp)
	call salloc (xr, npts, TY_REAL)
	call salloc (yr, npts, TY_REAL)
	call salloc (xr1, 2, TY_REAL)
	call salloc (yr1, 2, TY_REAL)
	call achtrr (x, Memr[xr], npts)
	call achtrr (y, Memr[yr], npts)
	call gt_copy (gt, gt1)		# gt1 is for plotting deleted points
	call gt_copy (gt, gt2)		# gt2 is for plotting neg. mag. points
	call gt_sets (gt1, GTTYPE, "mark")
	call gt_sets (gt1, GTMARK, "cross")
	call gt_sets (gt2, GTTYPE, "mark")
	call gt_sets (gt2, GTMARK, "plus")
	if (IC_LTYPE(ic) == IC_BOX) {
	    call gt_sets (gt, GTTYPE, "mark")
	    call gt_sets (gt, GTMARK, "box")
	    call gt_setr (gt, GTXSIZE, SIZE_BOX)
	    call gt_setr (gt, GTYSIZE, SIZE_BOX)
	} else {
	    call gt_sets (gt, GTTYPE, "line")
	}

	if (IC_OVERPLOT(ic) == NO) {
	    # Start a new plot.

	    call gclear (gp)

	    # Set the graph scale and axis. 
	    aux = Memr[yr]
	    i = nl_stati(nl,"fitfunc")
	    if ((IC_COMP(ic) == YES) && (IC_YAXIS(ic) == IC_LINEAR) && 
	        ((i == COMPOSITE) || (i == TWOBBODY)))
	        Memr[yr] = 0.
	    call gascale (gp, Memr[xr], npts, 1)
	    call gascale (gp, Memr[yr], npts, 2)
	    if (IC_YAXIS(ic) == IC_MAG) {
	        call ggwind (gp, x1, x2, y1, y2)    # flip y axis
	        call gswind (gp, x1, x2, y2, y1)
	    }
	    Memr[yr] = aux
	    call gt_swind (gp, gt)
	    call gt_labax (gp, gt)
	}

	# Plot data points.
	if (IC_OVERPLOT(ic) == NO) {
	    Memr[xr1] = Memr[xr]
	    Memr[yr1] = Memr[yr]
	    aux = IC_MAG0(ic) - 2.5 * log10 (IC_MINY(ic))
	    do i = 1, npts {
	        # this plots deleted points
		if (wts[i] == 0.) {  
		    call gt_plot (gp, gt1, Memr[xr+i-1], Memr[yr+i-1], 1) 
	        # this plots negative points in mag scale
		} else if ((IC_LTYPE(ic) == IC_BOX) && 
			   (IC_YAXIS(ic) == IC_MAG) &&
			   (abs (Memr[yr+i-1] - aux) < EPSILON) &&
			   (abs (ubars[i]) < EPSILON) && 
			   (abs (lbars[i]) < EPSILON) ){
	            Memr[yr+i-1] = aux
		    call gt_plot (gp, gt2, Memr[xr+i-1], Memr[yr+i-1], 1) 
	        # this plots regular points
		} else {
		    Memr[xr1+1] = Memr[xr+i-1]
		    Memr[yr1+1] = Memr[yr+i-1]
		    call gt_plot (gp, gt, Memr[xr1], Memr[yr1], 2)
		    Memr[xr1] = Memr[xr1+1]
		    Memr[yr1] = Memr[yr1+1]
		}
	    }
	}

	# See if wheights are not all the same
	wflag = false
	if (IC_LTYPE(ic) == IC_BOX) {
	    do i = 1, npts {
	        if ((wts[i] != 0.) && (wts[i] != 1.))
	            aux = wts[i]
	    }
	    do i = 1, npts {
	        if ((wts[i] != 0.) && (wts[i] != 1.) && (wts[i] != aux))
	            wflag = true
	    }
	}

	# Plot error bars
	if (wflag) {
	    call gt_sets (gt, GTTYPE, "mark")
	    if (IC_GKEY(ic) == 2)
	        call gt_sets (gt, GTMARK, "hebar")
	    else
	        call gt_sets (gt, GTMARK, "vebar")
	    if (IC_OVERPLOT(ic) == NO) {
	        do i = 1, npts {
		    if (wts[i] != 0.) {
	                # upper
	                if (IC_GKEY(ic) == 2) {
	                    call gt_setr (gt, GTXSIZE, -ubars[i])
	                    call gt_plot (gp, gt, Memr[xr+i-1] - ubars[i]/2., 
					  Memr[yr+i-1], 1)
	                } else {
	                    call gt_setr (gt, GTYSIZE, -ubars[i])
	                    call gt_plot (gp, gt, Memr[xr+i-1], 
					  Memr[yr+i-1] - ubars[i]/2., 1)
	                }
	                # lower
	                if (IC_GKEY(ic) == 2) {
	                    call gt_setr (gt, GTXSIZE, -lbars[i])
	                    call gt_plot (gp, gt, Memr[xr+i-1] + lbars[i]/2., 
					  Memr[yr+i-1], 1)
	                } else {
	                    call gt_setr (gt, GTYSIZE, -lbars[i])
	                    call gt_plot (gp, gt, Memr[xr+i-1], 
					  Memr[yr+i-1] + lbars[i]/2., 1)
	               }
		    }
	        }
	    }
	}

	# Reset status flags.
	IC_OVERPLOT(ic) = NO

	call gt_free (gt1)
	call gt_free (gt2)
	call sfree (sp)
end


# Mark sample points.

procedure icg_g2 (ic, gp, gt, x, y, npts, size)

pointer	ic			# ICFIT pointer
pointer	gp			# GIO pointer
pointer	gt			# GTOOLS pointer
real	x[npts], y[npts]	# Data points
int	npts			# Number of symbols to plot
real	size			# Symbol size

int	i
pointer	sp, xr, yr, gt1

begin
	call smark (sp)
	call salloc (xr, npts, TY_REAL)
	call salloc (yr, npts, TY_REAL)
	call achtrr (x, Memr[xr], npts)
	call achtrr (y, Memr[yr], npts)

	call gt_copy (gt, gt1)
	call gt_sets (gt1, GTTYPE, "mark")

	# Mark the sample points.

	if (abs (IC_NAVERAGE(ic)) > 1) {
	    call gt_sets (gt1, GTMARK, "plus")
	    call gt_setr (gt1, GTXSIZE, -size)
	    call gt_setr (gt1, GTYSIZE, 1.)
	    call gt_plot (gp, gt1, Memr[xr], Memr[yr], npts)
	}

	# Mark the rejected points.

	if (IC_NREJECT(ic) > 0) {
	    call gt_sets (gt1, GTMARK, "diamond")
	    call gt_setr (gt1, GTXSIZE, SIZE_REJECT)
	    call gt_setr (gt1, GTYSIZE, SIZE_REJECT)
	    do i = 0, npts-1 {
		if (Memi[IC_REJPTS(ic)+i] == YES)
	            call gt_plot (gp, gt1, Memr[xr+i], Memr[yr+i], 1)
	    }
	}

	call gt_free (gt1)
	call sfree (sp)
end

# ICG_GF -- Plot fitted function and eventual components.

procedure icg_gf (ic, gp, gt, nl, npts)

pointer	ic			# ICFIT pointer
pointer	gp			# GIO pointer
pointer	gt			# GTOOL pointer
pointer	nl			# NLFIT pointer
int	npts			# Number of points to plot

pointer	sp, x, coef, cflag, errors, aux
int	i, ifunc, npar
real	dx

int	nl_stati()

errchk	nl_init, icg_gf1

begin
	if (IC_FITERROR(ic) == YES)
	    return

	call smark (sp)
	call salloc (x,  npts, TY_REAL)

	# Generate vector of independent variable values
	dx = (IC_XMAX(ic) - IC_XMIN(ic)) / (npts - 1)
	do i = 1, npts
	    Memr[x+i-1] = IC_XMIN(ic) + (i-1) * dx

	# Draw fitted function
	if (IC_LTYPE(ic) == IC_CONTINUOUS)
	    call icg_gf1 (ic, gp, gt, nl, Memr[x], npts, GL_DOTTED)
	else
	    call icg_gf1 (ic, gp, gt, nl, Memr[x], npts, GL_SOLID)

	# Draw individual function components.
	ifunc = nl_stati (nl, "fitfunc")
	if (((ifunc == COMPOSITE) || (ifunc == TWOBBODY) || 
	     (ifunc == GALPROF)   ||
	     (ifunc == GAUSSIANS) || (ifunc == CGAUSS)) && 
	     (IC_COMP(ic) == YES)) {

	    # Save coefficients and create auxiliary coefficient area
	    npar = nl_stati (nl, "npar")
	    call salloc (coef, npar, TY_REAL)
	    call salloc (errors, npar, TY_REAL)
	    call salloc (cflag, npar, TY_BOOL)
	    call nl_gcoeff (nl, Memr[coef], npar)
	    call nl_gerrors (nl, Memr[errors], npar)
	    call nl_gflags (nl, Memb[cflag], npar)
	    call salloc (aux, npar, TY_REAL)

	    if (ifunc == COMPOSITE) {

	        Memr[aux + NL_PINDEX] = Memr[coef + NL_CINDEX]
	        Memr[aux + NL_PAMPL]  = Memr[coef + NL_CPAMPL]
	        Memr[aux + NL_PREF]   = Memr[coef + NL_CPREF]
	        call nl_init (nl, POWERLAW, Memr[aux], Memb[cflag], 3, npts)
	        call icg_gf1 (ic, gp, gt, nl, Memr[x], npts, GL_DASHED)
	        Memr[aux + NL_BTEMP] = Memr[coef + NL_CTEMP]
	        Memr[aux + NL_BAMPL] = Memr[coef + NL_CBAMPL]
	        Memr[aux + NL_BREF]  = Memr[coef + NL_CBREF]
	        call nl_init (nl, BBODY, Memr[aux], Memb[cflag], 3, npts)
	        call icg_gf1 (ic, gp, gt, nl, Memr[x], npts, GL_DASHED)

	    } else if (ifunc == TWOBBODY) {

	        Memr[aux + NL_BTEMP] = Memr[coef + NL_TTEMP1]
	        Memr[aux + NL_BAMPL] = Memr[coef + NL_TAMPL1]
	        Memr[aux + NL_BREF]  = Memr[coef + NL_TREF1]
	        call nl_init (nl, BBODY, Memr[aux], Memb[cflag], 3, npts)
	        call icg_gf1 (ic, gp, gt, nl, Memr[x], npts, GL_DASHED)
	        Memr[aux + NL_BTEMP] = Memr[coef + NL_TTEMP2]
	        Memr[aux + NL_BAMPL] = Memr[coef + NL_TAMPL2]
	        Memr[aux + NL_BREF]  = Memr[coef + NL_TREF2]
	        call nl_init (nl, BBODY, Memr[aux], Memb[cflag], 3, npts)
	        call icg_gf1 (ic, gp, gt, nl, Memr[x], npts, GL_DASHED)

	    } else if (ifunc == GALPROF) {

	        call amovr (Memr[coef], Memr[aux], 6)
	        Memr[aux+NL_GBACKGR+1] = Memr[aux+NL_GBACKGR+1] / 2.
	        call nl_init (nl, GALBULGE, Memr[aux], Memb[cflag], 6, npts)
	        call icg_gf1 (ic, gp, gt, nl, Memr[x], npts, GL_DASHED)
	        call nl_init (nl, GALDISK, Memr[aux], Memb[cflag], 6, npts)
	        call icg_gf1 (ic, gp, gt, nl, Memr[x], npts, GL_DASHED)

	    } else if (ifunc == GAUSSIANS) {

	        call nl_init (nl, GAUSSIANS, Memr[coef], Memb[cflag], 5, npts)
	        Memr[aux + NL_GA] = Memr[coef + NL_GA]
	        Memr[aux + NL_GB] = Memr[coef + NL_GB]
	        do i = 1, (npar-2) / 3 {
	            Memr[aux + NL_GAMPL(1)] = Memr[coef + NL_GAMPL(i)]
	            Memr[aux + NL_GCENT(1)] = Memr[coef + NL_GCENT(i)]
	            Memr[aux + NL_GFWHM(1)] = Memr[coef + NL_GFWHM(i)]
	            call nl_scoeff (nl, Memr[aux], 5)
	            call icg_gf1 (ic, gp, gt, nl, Memr[x], npts, GL_DASHED)
	        }

	    } else if (ifunc == CGAUSS) {

	        call nl_init (nl, GAUSSIANS, Memr[coef], Memb[cflag], 5, npts)
	        call icg_gf1 (ic, gp, gt, nl, Memr[x], npts, GL_DASHED)
	        Memr[aux + NL_GA] = Memr[coef + NL_GA]
	        Memr[aux + NL_GB] = Memr[coef + NL_GB]
	        if (npar > 5) {
	            do i = 2, (npar-2) / 3 {
	                Memr[aux + NL_GAMPL(1)] = Memr[coef + NL_GAMPL(i)] * Memr[coef + NL_GAMPL(1)]
	                Memr[aux + NL_GCENT(1)] = Memr[coef + NL_GCENT(i)] + Memr[coef + NL_GCENT(1)]
	                Memr[aux + NL_GFWHM(1)] = Memr[coef + NL_GFWHM(i)]
	                call nl_scoeff (nl, Memr[aux], 5)
	                call icg_gf1 (ic, gp, gt, nl, Memr[x], npts, GL_DASHED)
	            }
	        }
	    }
	    call nl_init (nl, ifunc, Memr[coef], Memb[cflag], npar, npts)
	    call nl_serrors (nl, Memr[errors], npar)
	}

	call sfree (sp)
end


# ICG_GF1 -- Plot function or one individual component.

procedure icg_gf1 (ic, gp, gt, nl, x, npts, ltype)

pointer	nl, ic, gt, gp
int	ltype, npts, i
real	x[ARB]

pointer	sp, y, xr, yr, xo, yo, gt1, dummy

errchk	nl_vector

begin
	call smark (sp)
	call salloc (y,  npts, TY_REAL)
	call salloc (xr, npts, TY_REAL)
	call salloc (yr, npts, TY_REAL)
	call salloc (xo, npts, TY_REAL)
	call salloc (yo, npts, TY_REAL)
	call salloc (dummy, npts, TY_REAL)
	call amovkr (0., Memr[dummy], npts)

	# Calculate vector of function values. 
	call nl_vector (nl, x, Memr[dummy], Memr[y], npts)

	if (IC_YAXIS(ic) == IC_MAG) {
	    do i = 0, npts-1 {
	        if (Memr[y+i] >= IC_MINY(ic))
	            Memr[y+i] = IC_MAG0(ic) - 2.5 * log10 (Memr[y+i])
	        else
	            Memr[y+i] = IC_MAG0(ic) - 2.5 * log10 (IC_MINY(ic))
	    }
	}

	# Convert to user function or transpose axes.  Change type to reals
	# for plotting.
	call icg_axes (ic, gt, nl, 1, x, Memr[y], Memr[xo], npts)
	call icg_axes (ic, gt, nl, 2, x, Memr[y], Memr[yo], npts)
	call achtrr (Memr[xo], Memr[xr], npts)
	call achtrr (Memr[yo], Memr[yr], npts)

	# Plot
	call gt_copy (gt, gt1)
	call gt_sets (gt1, GTTYPE, "line")
	call gt_seti (gt1, GTLINE, ltype)
	call gt_plot (gp, gt1, Memr[xr], Memr[yr], npts)
	call gt_free (gt1)

	call sfree (sp)
end

