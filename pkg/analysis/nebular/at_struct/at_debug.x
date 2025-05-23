# Copyright restrictions apply - see st4gem$copyright.stsdas

include	"../at.h"

#---------------------------------------------------------------------7 Sep 97--
.help at_debug.x Mar96 nebular/at_struct
.ih
NAME
.nf
 at_debug - Dereference elements of the AT structure for viewing
  at_view - Print elements of the AT structure to STDERR
.fi
.ih 
DESCRIPTION
These routines are intended to facilitate application debugging by 
dereferencing the elements of the AT structure and printing the 
contents to STDERR.  A good use might be to place statements like 
the following into an application: 
.nf
	if (DEBUG)
	    call at_debug (at)
.fi
If DEBUG is defined in the source file to be "true" then the contents 
of the AT structure, (and the first information for all elements of 
the first dimension of the various arrays) will be printed to STDERR.
.endhelp
#-------------------------------------------------------------------------------
#  AT_DEBUG - 	Dereference elements of the SF structure for viewing.  

procedure at_debug (at)

# Arguments:
pointer	at		# I: atomic data structure

begin
	call at_view (AT_ATOM(at), AT_ION(at), AT_NLVL(at), AT_NTRANS(at), 
		WEIGHT(at,1), E_TRANS(at), L_TRANS(at), RAD_TR(at), 
		AT_TE(at), AT_TE_MIN(at), AT_TE_MAX(at), AT_TY_FUNC(at), 
		CURVE(at,1), COLL(at), COLL_TR(at), N_CRIT(at,1), POP(at), 
		EMISS(at))
end


#-------------------------------------------------------------------------------
#  AT_VIEW - 	Print the elements of the SF structure to STDERR.  

procedure at_view (atom, ion, n_lvl, n_trans, weight, e_tr, l_tr, r_tr, 
			te, te_min, te_max, ty_func, cv, coll, coll_tr, 
			n_crit, pop, emiss)

# Arguments:
int	atom			# I: atomic number
int	ion			# I: ion (spectrum) number
int	n_lvl			# I: no. atomic levels
int	n_trans			# I: no. transitions between levels
int	weight[n_lvl]		# I: statistical weights @level
double	e_tr[n_lvl]		# I: energy @level
double	l_tr[n_lvl,n_lvl]	# I: 
double	r_tr[n_lvl,n_lvl]	# I: radiative transition probability
real	te			# I: electron temperature
real	te_min			# I: min allowed electron temperature
real	te_max			# I: max allowed electron temperature
int	ty_func			# I: type of fitting function
pointer	cv[ARB]			# I: curve fit structures
double	coll[n_lvl,n_lvl]	# I: collision strengths
double	coll_tr[n_lvl,n_lvl]	# I: collision transition probabilities
double	n_crit[n_lvl]		# I: critical densities
double	pop[n_lvl]		# I: level populations
double	emiss[n_lvl,n_lvl]	# I: line emissivities

# Local variables:
int	i			# generic

begin
	call eprintf (
	"\nAtom: %d; Ion: %d; N_lvls: %d; N_trans: %d; Ty_func: %d\n")
	    call pargi (atom)
	    call pargi (ion)
	    call pargi (n_lvl)
	    call pargi (n_trans)
	    call pargi (ty_func)

	call eprintf ("T_e: %8.2f; T_e_min: %8.2f; T_e_max: %8.2f\n")
	    call pargr (te)
	    call pargr (te_min)
	    call pargr (te_max)

	call eprintf ("\ni  Weight  Energy   N_crit      Pop  \n")
	for (i=1; i<=n_lvl; i=i+1) {
	    call eprintf ("%d     %d  %8g %8g %8g\n")
		call pargi (i)
		call pargi (weight[i])
		call pargd (e_tr[i])
		call pargd (n_crit[i])
		call pargd (pop[i])
	}

	call eprintf ("\ni  Coll_Str  Coll_Prb   Rad_Prb  E_trans   Emiss\n")
	for (i=1; i<=n_lvl; i=i+1) {
	    call eprintf ("%d  %8g  %8g  %8g %8g %8g\n")
		call pargi (i)
		call pargd (coll[i,1])
		call pargd (coll_tr[i,1])
		call pargd (r_tr[i,1])
		call pargd (l_tr[i,1])
		call pargd (emiss[i,1])
	}

	call flush (STDERR)
end


