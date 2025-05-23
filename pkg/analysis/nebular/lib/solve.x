# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<error.h>
include <mach.h>
include	"../at.h"

define	DEBUG	false

#--------------------------------------------------------------------14 Dec 00--
.help solve.x Mar96 nebular/lib
.ih
NAME
               solve - Calculates level popl., line emissivities & N_crit 
coll_transition_prob - Compute collisional transition probs betw/atomic levels
  calc_critical_dens - Compute critical densities for each atomic level 
          pop_matrix - Populate the level population matrix
        solve_matrix - Invert a symmetric matrix and return the solution vector
        line_emissiv - Calculate line emissivities. 
.endhelp
#-------------------------------------------------------------------------------
#  SOLVE -	Calculate level populations and line emissivities for a 
#		specified N_e and T_e, as well as critical densities.  

procedure solve (at, n_e, t_e)

# Arguments:
pointer	at		# I: atomic data structure
real	n_e		# I: electron density
real	t_e		# I: electron temerature

# Declarations:
pointer	a1		# work array
pointer	dtl		# work arrays 
int	i		# generic
int	n_levels	# no. atomic energy levels for this ion
pointer	pop1		# storage array for level populations
pointer	sp		# top of stack memory
int	sz_2darray	# 
double	sum_pop		# normalization of level populations
real	t4		# scaled temp

# Memory management:
define	A1		Memd[a1]
define	Dtl		Memd[(dtl+$1-1)]
define	Pop1		Memd[(pop1+$1-1)]

define	SCALE_FACTOR	1.D+16

errchk	coll_pop, solve_matrix

begin
	if (AT_LOG_TE(at))
	    t4 = log10 (te)
	else
	    t4 = te / 1.e4

	n_levels = AT_NLVL(at)
	sz_2darray = n_levels*n_levels

	call smark (sp)
	call salloc (a1, sz_2darray, TY_DOUBLE)
	call salloc (dtl,  n_levels, TY_DOUBLE)
	call salloc (pop1, n_levels, TY_DOUBLE)

	# Evaluate collision strengths for specified T_e. 
	call coll_pop (CURVE(at,1), AT_NTRANS(at), COLL(at), n_levels, t4)

	# Calculate collisional transition probabilities. 
	call coll_transition_prob (COLL(at), WEIGHT(at), L_TRANS(at), 
					COLL_TR(at), n_levels, t_e)

	# Calculate critical densities. 
	call calc_critical_dens (RAD_TR(at), COLL_TR(at), N_CRIT(at,1), n_levels)

	if (DEBUG)
	    call at_debug (at)

	# Scale population matrix by arbitrary scale factor, & solve.
	call pop_matrix (COLL_TR(at), RAD_TR(at), A1, n_levels, n_e)

	call amulkd (A1, SCALE_FACTOR, A1, sz_2darray)
	call solve_matrix (A1, Dtl(1), n_levels)

	# Accumulate & normalize level populations. 
	Pop1(1) = 1.
	sum_pop = 1.
	do i = 2, n_levels {
	    Pop1(i) = Dtl(i) / Dtl(1)
	    sum_pop = sum_pop + Pop1(i)
	}
	call adivkd (Pop1(1), sum_pop, Pop1(1), n_levels)
	call amovd (Pop1(1), POP(at), n_levels)

	if (DEBUG)
	    call at_debug (at)

	# Calculate line emissivities. 
	call line_emissiv (POP(at), RAD_TR(at), L_TRANS(at), EMISS(at), n_levels)

	AT_TE(at) = t_e

	call sfree (sp)
end


#-------------------------------------------------------------------------------
#  COLL_TRANSITION_PROB - Compute collisional transition probabilities between 
#			each atomic level.  

procedure coll_transition_prob (coll, weight, energy_diff, prob, n_lvl, t_e)

# Arguments:
double	coll[n_lvl,n_lvl]		# I: collision strengths
int	weight[n_lvl]			# I: statistical weights
double	energy_diff[n_lvl,n_lvl]	# I: relative energy separation
double	prob[n_lvl,n_lvl]		# O: collision transition probabilities
int	n_lvl				# I: no. atomic energy levels
real	t_e				# I: electron tempertature

# Declarations:
int	i, j				# generic
double	kt 				# product of Boltzman's constant & T_e
real	sqrt_te				# sqrt of T_e

define	CQ	8.629d-6		# 

begin
	kt      = (1.38062d-16) * t_e
	sqrt_te = sqrt (t_e)

	# Compute collision transition probabilities. 
	do i = 1, n_lvl {
	    do j = 1, n_lvl {
		if (j > i) 
		    prob[i,j] = CQ * coll[j,i] / double (weight[i] * 
				sqrt_te) * exp (-energy_diff[j,i] / kt) 

		else if (j == i) 
		    prob[i,j] = 0.d0

		else if (j < i)
		    prob[i,j] = CQ * coll[i,j] / double (weight[i] * sqrt_te)
	    }
	}

end


#-------------------------------------------------------------------------------
#  CALC_CRITICAL_DENS - Compute critical densities for each atomic level.  

procedure calc_critical_dens (rad_prob, coll_prob, crit_density, n_lvl)

# Arguments:
double	rad_prob[n_lvl,n_lvl]		# I: radiative transition probabilities
double	coll_prob[n_lvl,n_lvl]		# I: collisional transition probabilities
double	crit_density[n_lvl]		# O: critical densities
int	n_lvl				# I: no. atomic energy levels

# Declarations:
double	c_sum, r_sum			# accumulators
int	i, j				# generic

begin
	do i = 2, n_lvl {
	    r_sum = 0.d0
	    c_sum = 0.d0
	    do j = 1, i-1 {
		r_sum = r_sum + rad_prob[i,j]
		c_sum = c_sum + coll_prob[i,j]
	    }

	    # Must loop over all j != i;	Fix 14-Dec-2000 by R.Rubin
	    for (j = i+1; j < n_lvl+1; j = j+1)
		c_sum = c_sum + coll_prob[i,j]

	    # Protect against divide by zero. 
	    if (c_sum > EPSILOND)
	    	crit_density[i] = r_sum / c_sum
	    else
	    	crit_density[i] = INDEFD
	}

end


#-------------------------------------------------------------------------------
#  POP_MATRIX -	Populate the level population matrix.  

procedure pop_matrix (coll_trans, rad_trans, pop, np, n_e)

# Arguments:
double	coll_trans[np,np]	# I: collisional transition probabilities
double	rad_trans[np,np]	# I: radiative transition probabilities
double	pop[np,np]		# O: populated matrix
int	np			# I: size of array dimensions
real	n_e			# I: electron density

# Declarations:
int	i, j			# generic

begin
	# Populate matrix to be solved.
	do i = 1, np {
	    do j = 1, np {
		if (j == i) 
		    pop[i,j] = 0.D0
		else if (j > i) 
		    pop[i,j] = n_e * coll_trans[i,j]
		else if (j < i)
		    pop[i,j] = n_e * coll_trans[i,j] + rad_trans[i,j]
	    }
	}
end


#-------------------------------------------------------------------------------
#  SOLVE_MATRIX - Invert a symmetric matrix and return the solution vector. 

procedure solve_matrix (a1, dtl, n_lvl)

# Arguments:
double	a1[n_lvl,n_lvl]		# I: matrix to be inverted
double	dtl[n_lvl]		# O: solution vector
int	n_lvl			# I: size of arrays

# Declarations:
pointer	a2			# work array
pointer	b1, b2			# work arrays 
int	i, j, k, l		# generic
pointer	sp			# top of stack memory
double	t1, t2			# accumulation variables

# Memory management
define	A2	Memd[(a2+($1-1)+(n_lvl*($2-1)))]	# Elements of a2 array
define	B1	Memd[(b1+$1-1)]		# Elements of b1 array
define	B2	Memd[(b2+$1-1)]		# Elements of b2 array

define	MIN_DOUBLE	1.d0 / MAX_DOUBLE

begin
	call smark (sp)
	call salloc (a2, n_lvl*n_lvl, TY_DOUBLE)
	call salloc (b1, n_lvl, TY_DOUBLE)
	call salloc (b2, n_lvl, TY_DOUBLE)

	do i = 1, n_lvl {
	    B1(i) = 0.d0
	    do j = 1, n_lvl
		B1(i) = B1(i) + a1[i,j]
	}

	do i = 1, n_lvl {
	    B2(i)   = -a1[1,i]
	    a1[i,i] = -B1(i)
	}

	do i = 1, n_lvl {
	    call amovd (a1, A2(1,1), n_lvl*n_lvl)

	    do j = 1, n_lvl
		A2(i,j) = B2(j)

	    dtl[i] = 1.
	    do j = 2, n_lvl-1 {
		t1 = A2(j,j)
		if (abs(t1) < MIN_DOUBLE) 
		    call error (2, "Floating underflow in routine SOLVE")

		dtl[i] = dtl[i] * t1
		do k = j, n_lvl
		    A2(j,k) = A2(j,k) / t1

		do k = j+1, n_lvl {
		    t2 = A2(k,j)
		    do l = j, n_lvl
			A2(k,l) = A2(k,l) - t2 * A2(j,l)
		}
	    }

	    dtl[i] = dtl[i] * A2(n_lvl,n_lvl)
	}

	call sfree (sp)
end


#-------------------------------------------------------------------------------
# LINE_EMISSIV - Calculate line emissivities. 

procedure line_emissiv (lvl_pop, radiative_tr, elevel_tr, emissiv, n_lvl)

# Arguments:
double	lvl_pop[n_lvl]		# I: level populations
double	radiative_tr[n_lvl,n_lvl] # I: radiative transition probabilities
double	elevel_tr[n_lvl,n_lvl]	# I: energy level transition matrix
double	emissiv[n_lvl,n_lvl]	# O: emissivities
int	n_lvl			# I: no. energy levels

# Declarations:
int	i, j			# generic

define	MIN_DOUBLE	1.D0 / MAX_DOUBLE

begin
	do i = 1, n_lvl {
	    do j = 1, n_lvl {
		emissiv[i,j] = lvl_pop[i] * radiative_tr[i,j] * elevel_tr[i,j]
		emissiv[i,j] = max (emissiv[i,j], MIN_DOUBLE)
	    }
	}
end


