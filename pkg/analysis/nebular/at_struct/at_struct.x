# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<error.h>
include "../at.h"
include "../neberr.h"

#--------------------------------------------------------------------06 Feb 98--
.help at_struct.x May96 nebular/at_struct
.ih
NAME
.nf
    at_alloc - Allocate an atomic data structure
     at_free - Free an atomic data structure
at_get_nlvls - Determine no. energy levels for atomic data structure
.fi
.endhelp
#-------------------------------------------------------------------------------
#  AT_ALLOC -- 	Allocate the atomic data structure.  If the structure exists, 
#		free it and then re-alloc for the new atom/ion.  

procedure at_alloc (atom, ion, at)

#  Arguments:
int	atom			# I: atomic number
int	ion			# I: ionization stage
pointer	at			# O: atomic data object

#  Declarations:
int	at_get_func_typ()	# get type of fitting function
int	at_get_nlvls()		# Determine no. energy levels 
int	get_atom_str()		# get atom/ion name
int	i, j			# generic
char	ion_str[SZ_LINE]	# name of atom/ion
int	n_lvl			# no. atomic energy levels for this ion
int	n_trans			# no. possible downward energy tranistions
pointer	open_atom_file()	# open atomic datafile
int	strdic()		# search a string dictionary
char	sx[SZ_LINE]		# generic
char	sx2[SZ_LINE]		# generic
int	sz_arr			# size of 2-D arrays
pointer	tpa			# table descriptor for generic atomic data
pointer	tpc			# table descriptor for collision strengths

errchk	read_atomic

begin
	if (at != NULL)
	   call at_free (at)

	# Allocate atomic data structure. 
	call calloc (at, LEN_AT, TY_STRUCT)
	AT_ATOM(at) = atom
	AT_ION(at)  = ion
	i = get_atom_str (atom, ion, ion_str, SZ_LINE)

	# Open the datafiles & determine no. atomic energy levels to expect. 
	# Failure to open a datafile means the atom/ion is not supported. 
	iferr {
	    tpa = open_atom_file (atom, ion, AT_DATA)
	    tpc = open_atom_file (atom, ion, AT_COLLISION)

	} then {
	    call mfree (at, TY_STRUCT)
	    i = get_atom_str (atom, ion, sx2, SZ_LINE)
	    call sprintf (sx, SZ_LINE, "Cannot locate atomic data for %s")
		call pargstr (sx2)
	    call error (INVLD_ATOM_ION, sx)
	}

	n_lvl = at_get_nlvls (tpa, tpc)
	if (n_lvl > 0)
	    AT_NLVL(at) = n_lvl

	else {
	    call mfree (at, TY_STRUCT)
	    call tbtclo (tpa)
	    call tbtclo (tpc)
	    call sprintf (sx, SZ_LINE, "Corrupted data for atom %s")
		call pargstr (ion_str)
	    call error (BAD_REF_DATA, sx)
	}

	# Determine the ground-state electron configuration.
	call tbhgtt (tpa, "GSCONFIG", sx2, SZ_LINE)
	call strupr (sx2)
	i = strdic (sx2, sx2, SZ_LINE, GS_CONFIG)
	if (i > 0)
	    AT_GSCONFIG(at) = i

	else {
	    call mfree (at, TY_STRUCT)
	    call sprintf (sx, SZ_LINE, "Corrupted data for atom %s")
		call pargstr (ion_str)
	    call error (BAD_REF_DATA, sx)
	}

	# Compute number of possible radiative transitions. 
	n_trans = 0
	do i = AT_NLVL(at), 1, -1 {
	    do j = 1, i-1
	    	n_trans = n_trans + 1
	}
	AT_NTRANS(at) = n_trans

	# Allocate arrays.
	sz_arr = n_lvl * n_lvl 
	call calloc (AT_WT(at),     n_lvl,  TY_INT)
	call calloc (AT_E_TR(at),   n_lvl,  TY_DOUBLE)
	call calloc (AT_L_TR(at),   sz_arr, TY_DOUBLE)
	call calloc (AT_RAD_TR(at), sz_arr, TY_DOUBLE)

	# Read T_e-independent atomic data from file. 
	call read_atomic (tpa, at)
	call tbtclo (tpa)

	# Fetch the T_e dependent fitting function. 
	call tbhgtt (tpc, "FIT_FUNC", sx, SZ_LINE)
	AT_TY_FUNC(at) = at_get_func_typ (sx, SZ_LINE)

	# Allocate arrays for T_e dependent atomic data. 
	AT_TE(at)     = INDEFR
	AT_TE_MIN(at) = INDEFR
	AT_TE_MAX(at) = INDEFR

	call calloc (AT_CV(at),      n_trans, TY_POINTER)
	call calloc (AT_COLL(at),    sz_arr,  TY_DOUBLE)
	call calloc (AT_COLL_TR(at), sz_arr,  TY_DOUBLE)
	call calloc (AT_NCRIT(at),   n_lvl,   TY_DOUBLE)
	call calloc (AT_POP(at),     n_lvl,   TY_DOUBLE)
	call calloc (AT_EMISS(at),   sz_arr,  TY_DOUBLE)

	# Read T_e dependent atomic data from file. 
	iferr ( call fit_collisions (tpc, at) ) {
	    call erract (EA_WARN)
	    call tbtclo (tpc)
	    call at_free (at)
	    i = get_atom_str (atom, ion, sx2, SZ_LINE)
	    call sprintf (sx, SZ_LINE, 
		"Error fitting collision strengths for atom %s")
		call pargstr (sx2)
	    call error (BAD_REF_DATA, sx)
	}
	call tbtclo (tpc)

	# Compute the transition energy differences. 
	call calc_energy_diff (E_TRANS(at), L_TRANS(at), AT_NLVL(at))
end


#-------------------------------------------------------------------------------
#  AT_FREE -- Free an atomic data structure.

procedure at_free (at)

#  Calling arguments:
pointer	at			# I: atomic data object

#  Declarations:
int	i			# generic

begin
	# Free T_e dependent data. 
	do i = 1, AT_NTRANS(at) {
	    if (CURVE(at,i) != NULL)
	    	call cvfree (CURVE(at,i))
	}

	call mfree (AT_CV(at),      TY_POINTER)
	call mfree (AT_COLL(at),    TY_DOUBLE)
	call mfree (AT_COLL_TR(at), TY_DOUBLE)
	call mfree (AT_NCRIT(at),   TY_DOUBLE)
	call mfree (AT_POP(at),     TY_DOUBLE)
	call mfree (AT_EMISS(at),   TY_DOUBLE)

	# Free T_e independent data. 
	call mfree (AT_WT(at),     TY_INT)
	call mfree (AT_E_TR(at),   TY_DOUBLE)
	call mfree (AT_L_TR(at),   TY_DOUBLE)
	call mfree (AT_RAD_TR(at), TY_DOUBLE)

	call mfree (AT_CV(at), TY_POINTER)
	call mfree (at, TY_STRUCT)
end


#-------------------------------------------------------------------------------
#  AT_GET_NLVLS - Determine no. energy levels for atomic data structure from 
#		two input tables. 

int procedure at_get_nlvls (tpa, tpc)

#  Arguments:
pointer	tpa		# table descriptor for generic atomic data
pointer	tpc		# table descriptor for collision strengths

#  Declarations:
int	nlvl_a		# no. atomic levels available from generic atomic data
int	nlvl_c		# no. atomic levels available from collision strengths
int	tbhgti()	# fetch table header parameter, TY_INT

begin
	nlvl_a = tbhgti (tpa, "N_LEVELS")
	nlvl_c = tbhgti (tpc, "N_LEVELS")

	return (min( nlvl_a, nlvl_c ))
end


