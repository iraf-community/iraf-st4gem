###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	procedure sfinteract(nfree, fpar, chisq)
#		int	nfree
#		real	fpar[ARB]
#		real	chisq
#
#  Description:	SFINTERACT uses IRAF plot routines to plot a spectrum
#		on STDGRAPH.  Fitting the data is driven by cursor commands.
#		
#  Arguments:	int	nfrees		- Number of free parameters
#		real	fpar[ARB]	- Array of free parameters
#
#  Returns:	real	chisq		- Last evaluated value of chi-square
#
#  Notes:	Adapted from splot in the IRAF onedspec package
#
#  History:	May	1989	Gerard Kriss
#		March	1991	G. Kriss - added d, -, +, and window options
#		Oct	1991	G. Kriss - added print, help,
#						 and change step options
#		June 	1994	Grimes - added Insert (I), List(G), Subtract 
#					(S), What to PLot(o) and modified 
#					Plot(p)
#		7/11/94		Grimes - Put in some error checking
#		8/15/94		Grimes - Changed defaults to old values when
#					entering a param value
#		8/25/94		Grimes - Added Numrecipe (n)
#		4/20/95		Grimes - Added clear screen in sfaddcomp
#		April 95	Grimes - Added new way to access key file
#
###########################################################################

define	RAWDATA	1
define	FITDATA	2
define	RESIDUALS	3
#define  KEY             "/export/home/student/grimes/specfit/specfit8/specfit.key"
define  PROMPT          "specfit options"

include	<gset.h>
include	<gio.h>
include	<pkg/gtools.h>
include	"specfit.h"

procedure sfinteract(nfree, fpar, chisq)
int	nfree
real	fpar[ARB]
real	chisq

bool	first
char	command[SZ_FNAME]
int	wc, key, newgraph, datatype
int	newcomp, newpar, nc, i, junk
real	wx, wy, dx, var
real	diff[MAXPTS]
pointer	gt, gfd, ranges, plot_ranges

int	gt_gcur(), clgeti(), clgstr(), strlen()
real	clgetr()
pointer	gopen(), gt_init()
int	sfcomponent_add()
char	key_file[SZ_FNAME]


include	"specfit.com"

begin
	junk = clgstr("key_file",key_file,SZ_FNAME)
	if (key_file[strlen(key_file)] == '/') {
		call strcat("specfit.key",key_file,SZ_FNAME)
	} else {
		call strcat("/specfit.key",key_file,SZ_FNAME)
	}

	#Set default plotting to include all components
	call clpstr("plot_range","*")

	# Open plotter and eliminate y-axis minor ticks.
	gfd = gopen ("stdgraph", NEW_FILE, STDGRAPH)
	call gseti (gfd, G_YNMINOR, 0)

	# Initialize graph format
	gt = gt_init()
	call gt_sets (gt, GTTITLE, specname)

	# Initial plot of data array

	datatype = RAWDATA
	call replot (gfd, gt, Memr[lambdas], Memr[spectrum], npts)

	# Evaluate current fit and overplot it
	call chispec(nfree, fpar, chisq)
	call printf("Initial chisq=%f\n")
	call pargr(chisq)
	call overplot(gfd, gt, Memr[lambdas], Memr[fitsp], npts)
	newgraph = NO
	first = true

	while (gt_gcur ("cursor", wx, wy, wc, key, command, SZ_FNAME) != EOF) {
	    switch (key) {
	    case 'I':		#interactive component addition
		
		if (YES==sfcomponent_add(nfree,fpar,gfd,gt)) {
			call replot(gfd,gt,Memr[lambdas],Memr[spectrum],npts)
			call overplot(gfd,gt,Memr[lambdas],Memr[fitsp],npts)
		}

	    case 'S':
		call sfcomponent_del(nfree,fpar)  #interactive component deletion


	    case 'G':
		call sfcomponent_list()		#list the component names

	    case 'a':
		newcomp = clgeti("comp_number")
		newpar = clgeti("param_number")
		call check_vals(newcomp,newpar)
		call clputr("param_tolerance",ptol[parptr[newpar,newcomp]])
		ptol[ parptr[ newpar, newcomp ] ] = clgetr("param_tolerance")
		call setpar(nfree, fpar)

	    case 'c':
		newcomp = clgeti("comp_number")
		newpar = clgeti("param_number")
		call check_vals(newcomp,newpar)
		call clputr("new_par_value",par0[parptr[newpar,newcomp]])
		par0[ parptr[ newpar, newcomp ] ] = clgetr("new_par_value")
		call setpar(nfree, fpar)

	    case 'e':

		call chispec(nfree, fpar, chisq)
		call printf("Current chisq = %g for %d parameters %d free and %d points\n")
		call pargr(chisq)
		call pargi(npar)
		call pargi(nfree)
		call pargi(nfitpts)

	    case 'f':
		call sfdosimplex(nfree, fpar, chisq)


	    case 'm':
		call sfdomarquadt(nfree, fpar, chisq)

	    case 'n':
		call sfdonumrec(nfree, fpar, chisq)
		
	    case 'i':
		newcomp = clgeti("comp_number")
		newpar = clgeti("param_number")
		call check_vals(newcomp,newpar)
		call printf("\n")
		call printf("%s%d %10.4g %10.4g %10.4g %10.4g %10.4g %2d %s\n")
		    call pargstr(compkeys[1,comtyp[newcomp]])
		    call pargi(newcomp)
		    call pargr(par0[ parptr[newpar, newcomp] ])
		    call pargr(blim[ parptr[newpar, newcomp] ])
		    call pargr(tlim[ parptr[newpar, newcomp] ])
		    call pargr(step[ parptr[newpar, newcomp] ])
		    call pargr(ptol[ parptr[newpar, newcomp] ])
		    call pargi(ifix[ parptr[newpar, newcomp] ])
		    call pargstr(comments[1,parptr[newpar,newcomp]])

	    case 'l':
		newcomp = clgeti("comp_number")
		newpar = clgeti("param_number")
		call check_vals(newcomp,newpar)
		call clputr("lower_limit",blim[parptr[newpar,newcomp]])
		blim[ parptr[ newpar, newcomp ] ] = clgetr("lower_limit")

	    case 'p':
		
		call chispec(nfree,fpar,chisq)				
		call overplot(gfd, gt, Memr[lambdas], Memr[fitsp], npts)

	    case 'o':

                call malloc(plot_ranges, 12 * MAXSAMPLE, TY_CHAR)
                nc = clgstr("plot_ranges",Memc[plot_ranges],12 * MAXSAMPLE - 1)
                call decode_plot ( Memc[plot_ranges] )
                call mfree(plot_ranges, TY_CHAR)


	    case 'r':
		datatype = RAWDATA
		newgraph = YES

	    case 's':		# Select new sample ranges
		call malloc(ranges, 12*MAXSAMPLE, TY_CHAR)
		nc = clgstr("new_sample", Memc[ranges], 12*MAXSAMPLE-1)
		call decode_ranges( Memc[ranges] )
		call mfree(ranges, TY_CHAR)
	
	    case 't':
		call clputi("change_iterations",itr)
		itr = clgeti("change_iterations")

	    case 'u':
		newcomp = clgeti("comp_number")
		newpar = clgeti("param_number")
		call check_vals(newcomp,newpar)
		call clputr("upper_limit",tlim[parptr[newpar,newcomp]])
		tlim[ parptr[ newpar, newcomp ] ] = clgetr("upper_limit")

	    case 'x':
		newcomp = clgeti("comp_number")
		newpar = clgeti("param_number")
		call check_vals(newcomp,newpar)
		call clputi("fix_or_free",ifix[parptr[newpar,newcomp]])
		ifix[ parptr[ newpar, newcomp ] ] = clgeti("fix_or_free")
		call setpar(nfree, fpar)

	    case 'z':
		newcomp = clgeti("comp_number")
		newpar = clgeti("param_number")
		call check_vals(newcomp,newpar)
		call clputr("step_size",step[parptr[newpar,newcomp]])
		step[ parptr[ newpar, newcomp ] ] = clgetr("step_size")

	    case 'd':	# Plot the distribution of Chi-square
		for ( i = 1; i <= npts; i = i + 1 ) {
		    if ( Memi[infit+i-1] == 0 ) {
		    	diff[i] = 0.0
		    } else {
			dx = (Memr[spectrum+i-1] - Memr[fitsp+i-1])
			if ( err_from_model ) {
				var = Memr[fitsp+i-1]
			} else {
				var = Memr[errors+i-1] * Memr[errors+i-1]
			}
			if ( var > 0. ) {
				diff[i] = dx * dx / var
			} else {
			 	diff[i] = 1.e10
			}
		    }
		}
		datatype = RESIDUALS
		call replot (gfd, gt, Memr[lambdas], diff, npts)

	    case '-':	# Plot the residual spectrum
		for ( i = 1; i <= npts; i = i + 1 ) {
			diff[i] = Memr[spectrum+i-1] - Memr[fitsp+i-1]
			if ( Memi[infit+i-1] == 0 ) diff[i] = 0.0
		}
		datatype = RESIDUALS
		call replot (gfd, gt, Memr[lambdas], diff, npts)

	    case '+':	# Plot the model alone
		datatype = FITDATA
		newgraph=YES

	    case 'w': # Window the graph
		call gt_window (gt, gfd, "cursor", newgraph)

	    case '?':
		call gpagefile(gfd, key_file, PROMPT)

	    case ':':
		if (command[1] == '/')
		    call gt_colon (command, gt, newgraph)
		else
		    call printf ("07")

	    case ' ':
		call printf ("x,y: %10.3f %10.4g")
		    call pargr (wx)
		    call pargr (wy)
		call flush (STDOUT)

	    default:
		# Default = ' '
		call printf ("x,y: %10.3f %10.4g")
		    call pargr (wx)
		    call pargr (wy)
		call flush (STDOUT)

	    }


	    if (newgraph == YES) {
	     switch ( datatype ) {
		case RAWDATA: 
		    call replot (gfd, gt, Memr[lambdas], Memr[spectrum], npts)
		case FITDATA:
		    call replot (gfd, gt, Memr[lambdas], Memr[fitsp], npts)
		case RESIDUALS:
		    call replot (gfd, gt, Memr[lambdas], diff, npts)
		default:
		    call replot (gfd, gt, Memr[lambdas], Memr[spectrum], npts)
		    datatype = RAWDATA
		}
		newgraph = NO
	    }
	}


# Close up graph window
	call gclear(gfd)
	call gclose(gfd)
	call gt_free(gt)
end



# This short procedure is just used to make sure that only legal values
# are entered for the component numbers and parameter numbers.  If one is
# out of range it is given the maximum possible value


procedure check_vals(comp,par)
int comp,par

include "specfit.com"

begin
if ( comp > ncomp )
	comp = ncomp
if ( par > ncpar[comtyp[comp]] )
	par = ncpar[comtyp[comp]]
end






















