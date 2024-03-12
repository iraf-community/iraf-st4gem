###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	call simplex(npar,par,step,maxerr,tol,maxiter,fx,chisq)
#
#  Description:	SIMPLEX is functional minimization routine which uses the
#		simplex algorithm. The user supplies the function to be
#		minimized.
#
#  Arguments:	int  npar	Number of parameters
#		real par[ARB]	Parameter initial values
#		real step[ARB]	Step size forthe search
#		real maxerr[ARB] Tolerance on individual parameters
#		real tol	Tolerance on Chi-square
#		int  maxiter	Maximum number of iterations
#		pointer fx	Pointer to the function to be minimized
#
#  Returns:	real chisq	Value of Chi-square after final iteration
#
#  Notes:	Information shared in common blocks defined in "specfit.com".
#		Most of the routines are written in SPP.
#
#  History:	1982		BYTE magazine publishes pascal version
#		1985		Gerard Kriss
#				Translated to "C"
#		May 1989	Gerard Kriss
#				Re-wrote in SPP for use in IRAF
#		10/12/89	G. Kriss
#				Change iteration logic
#		10/19/89	Change limit logic
#		7/9/94		J Grimes Changed initial guess algorithmn
#				to fix the tendency of the algorithmn to
#				think its found a false minimum
#
###########################################################################

include	"specfit.h"

define ROOT2 1.414213562
define SIN45 0.707106781
define REF -1.0
define CON 0.25
define EXP -4.0

real procedure simplex(npar,par,step,maxerr,tol,maxiter,fx,chisq)
int	npar, maxiter
real	par[ARB],step[ARB],maxerr[ARB],tol
extern	fx
real	chisq


int	high, low, expon
real	vertex[MAXFREE, MAXFREE+1], value[MAXFREE+1], nextv[MAXFREE+1]
real	trial, chisq2

int	i, j, niter, miter, limiter, done
real	p[MAXFREE], q[MAXFREE], qstep, pstep, trial2
real	error[MAXFREE], diff, val
real	center[MAXFREE], delta[MAXFREE]
real	sqrt()

bool	err_from_model
bool	interact
bool	debug
int	nlogfd
int	logfd[5]

common /users/ err_from_model, interact, debug, nlogfd, logfd
begin
#
# Outer loop iterations are controlled by the user.  The inner loop is limited
# to a number of iterations that are four times the number of parameters.
# This prevents the step sizes from getting too small when still far from a
# true minimum and decreasing the efficiency of the search.
#
    done = NO
    miter = 0
    limiter = 4 * npar
    expon = 0
    chisq = 0
    chisq2 = -100
    while (miter < maxiter && done == NO) {
	miter = miter + 1

# Before initialization, restrict parameters to the user-defined range
	call setlim(npar, par)

	qstep = (sqrt(npar + 1.) - 1.)/(npar * ROOT2)
	pstep = qstep + SIN45
	qstep = qstep / ( 10 ** expon )
	pstep = pstep / ( 10 ** expon )

	for( i = 1; i <= npar; i = i + 1) {  #Get offsets for starting vertices 
	    p[i] = step[i] * pstep
	    q[i] = step[i] * qstep
	}

	call fx(npar, par, val)    # Compute starting vertices
	value[1] = val
	for(i = 2; i <= npar+1; i = i + 1) {  
	    vertex[i-1, 1] = par[i-1]
	    for(j = 1; j <= npar; j = j + 1)
	    	vertex[j, i] = par[j] + q[j]
	    vertex[i-1, i] = par[i-1] + p[i-1]
	    call fx(npar, vertex[1, i], val)    # Compute starting vertices
	    value[i] = val
	}
	low = 1
	high = 1
	call extremes(npar, low, high, value)

# Setup is done.  Enter the main search section.
	niter = 0
	while (niter < limiter && done == NO) {
	    niter = niter + 1

	    for(i = 1; i <= npar; i = i + 1) {
	        center[i] = 0.0
	    }

#	    Compute centroid excluding worst vertex
	    for(i = 1; i <= npar+1; i = i + 1) {
	    	if( i != high) {
	    	    for(j = 1; j <= npar; j = j + 1)
	                center[j] = center[j] + vertex[j, i]
		}
	    }

#	    Next vertex is reflection of the worst
	    for(i = 1; i <= npar; i = i + 1) {
	    	center[i] = center[i] / npar
	    	delta[i] = vertex[i, high] - center[i]
	    	nextv[i] = center[i] + (REF * delta[i])
	    }
	    call fx(npar, nextv, val)
	    trial = val

	    if(trial <= value[low]) { # Better than current best?
	    	call new_vertex(trial, npar, high, nextv, value, vertex)  # Take it!
	    	for(i = 1; i <= npar; i = i + 1)   # Try an expansion
	    	    nextv[i] = center[i] + (EXP * delta[i])
	        call fx(npar, nextv, val)
	        trial2 = val
	    	if(trial2 <= trial) {  # If better yet, take it.
	    	    trial = trial2
	    	    call new_vertex(trial, npar, high, nextv, value, vertex)
	    	}
	    } else {			 # If not better than the best
	    	if(trial <= value[high]) { # accept if better than worst
	    	    call new_vertex(trial, npar, high, nextv, value, vertex)
	    	} else {			# If worse than worst, contract
	    	    for(i = 1; i <= npar; i = i + 1)
	    	    	nextv[i] = center[i] + CON * delta[i]
	            call fx(npar, nextv, val)
	            trial = val
	    	    if(trial <= value[high]) { # If contraction improves things
#						 accept it.
	    	    	call new_vertex(trial, npar, high, nextv, value, vertex)
	    	    } else {			# Things here are still bad.
#						   Shrink all vertices. 
	    	    	for(i = 1; i<=npar+1; i = i + 1) {
	    	    	    for(j = 1; j <= npar; j = j + 1) {
	    	    	    	vertex[j, i] = vertex[j, low] +
				     CON * (vertex[j, i] - vertex[j, low])
			    }
			    call fx(npar, vertex[1, i], val)
	    	    	    value[i] = val
	    	    	}
	    	    }
	    	}
	    }
	    call extremes(npar, low, high, value)

	    if ( debug ) {
		call printf("miter %4d, niter %4d. Chi = %10.2f\n")
		  call pargi(miter)
		  call pargi(niter)
		  call pargr(value[low])
	    }

	    done = YES
	    diff = (value[high] - value[low])/value[low]
	    if (diff > tol)
	    	done = NO
	    for(j = 1; j <= npar && done == YES; j = j + 1) {
	    	error[j] = (vertex[j, high] -
		            vertex[j, low]) / vertex[j, low]
	    	if(error[j] > maxerr[j])
	    	    done = NO
	    }
	}
	for(i = 1; i <= npar; i = i + 1)
	    par[i] = vertex[i, low]
	chisq2 = chisq
	chisq = value[low]
	call printf("%d %15.7g\n")	# Report on progress
		call pargi(miter)
		call pargr(chisq)

	if ( chisq == chisq2 )
		{ expon = expon + 1 } 
  }
end

# new_vertex: A routine to substitute vertices in simplex. 
procedure new_vertex(trial, npar, high, nextv, value, vertex)
real	trial
int	npar, high
real	nextv[ARB], value[ARB], vertex[MAXFREE, MAXFREE+1]

int	i

begin
	for(i = 1; i <= npar; i = i + 1) {
	    vertex[i, high] = nextv[i]
	}
	value[high] = trial
end

# extremes: A routine to find the highest and lowest values of each
#	    	parameter associated with a vertex in simplex. 
procedure extremes(npar, low, high, value)
int	npar, low, high
real	value[ARB]

int	i

begin
	for( i = 1; i <= npar + 1; i = i + 1) {
	    if (value[i] < value[low])
	    	low = i
	    if (value[i] > value[high])
	    	high = i
	}
end
