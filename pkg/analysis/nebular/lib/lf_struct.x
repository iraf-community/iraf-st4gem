# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<error.h>
include "../flux.h"
include "../neberr.h"

#--------------------------------------------------------------------06 Feb 98--
.help lf_struct.x Feb97 nebular/lib
.ih
NAME
.nf
    lf_alloc - Allocate a flux data structure
     lf_free - Free a flux object data structure
 lf_add_line - Add observed emission line strength to LF object. 
    lf_debug - Dereference values from the flux object for viewing
.fi
.endhelp
#-------------------------------------------------------------------------------
#  LF_ALLOC -- 	Allocate the flux data structure.  If the structure exists, 
#		free it and then re-alloc for the new atom/ion.  

procedure lf_alloc (atom, ion, o)

#  Arguments:
int	atom		# I: atomic number
int	ion		# I: ionization stage
pointer	o		# O: flux data object

#  Declarations:
int	sz_arr		# size of observed flux arrays

begin
	if (o != NULL)
	   call lf_free (o)

	# Allocate atomic data structure. 
	call calloc (o, LEN_LF, TY_STRUCT)
	if (o == NULL)
	    call error (1, "Memory allocation failed for LF")
	LF_ATOM(o)  = atom
	LF_ION(o)   = ion
	LF_ZONE(o)  = INDEFI
	LF_ABUND(o) = INDEFR
	LF_N(o)     = 0

	# Allocate arrays.
	sz_arr       = LF_GROW 
	LF_SZ_ARR(o) = sz_arr
	call calloc (LF_WAVE_PTR(o),   sz_arr, TY_REAL)
	call calloc (LF_FLUX_PTR(o),   sz_arr, TY_REAL)
	call calloc (LF_WT_PTR(o),     sz_arr, TY_REAL)
	call calloc (LF_NLINES_PTR(o), sz_arr, TY_INT)
	call calloc (LF_WIDTH_PTR(o),  sz_arr, TY_REAL)

end


#-------------------------------------------------------------------------------
#  LF_FREE -- Free flux object data structure.

procedure lf_free (o)

#  Calling arguments:
pointer	o		# I: flux data object

begin
	if (o != NULL) {
	    call mfree (LF_WAVE_PTR(o),   TY_REAL)
	    call mfree (LF_FLUX_PTR(o),   TY_REAL)
	    call mfree (LF_WT_PTR(o),     TY_REAL)
	    call mfree (LF_NLINES_PTR(o), TY_INT)
	    call mfree (LF_WIDTH_PTR(o),  TY_REAL)

	    call mfree (o, TY_STRUCT)
	}
end


#-------------------------------------------------------------------------------
#  LF_SET_ZONE -- Set zone value in flux object data structure.

procedure lf_set_zone (o, zone)

#  Calling arguments:
pointer	o		# I: flux data object
int	zone		# I: zone number

begin
	if (o != NULL) {
	    if (zone > 0)
	    	LF_ZONE(o) = zone
	}
end


#-------------------------------------------------------------------------------
#  LF_ADD_LINE - Add observed emission line strength to LF array. Returns 
#		"yes" if line was added.  

int procedure lf_add_line (o, wave, flux, weight, nlines, width)

#  Arguments:
pointer	o		# I: flux data object
real	wave		# I: nominal wavelength of observed multiplet
real	flux		# I: observed flux in multiplet
real	weight		# I: relative weight of flux for abundance determination
int	nlines		# I: no. lines in multiplet
real	width		# I: relative weight of flux for abundance determination

#  Declarations:
int	n		# generic

begin
	# Add line only if flux is not INDEFR.
	if (IS_INDEFR (flux))
	    return (NO)

	n = LF_N(o) + 1
	LF_N(o) = n

	if (n > LF_SZ_ARR(o)) {
	    LF_SZ_ARR(o) = LF_SZ_ARR(o) + LF_GROW
	    call realloc (LF_WAVE_PTR(o),   LF_SZ_ARR(o), TY_REAL)
	    call realloc (LF_FLUX_PTR(o),   LF_SZ_ARR(o), TY_REAL)
	    call realloc (LF_WT_PTR(o),     LF_SZ_ARR(o), TY_REAL)
	    call realloc (LF_NLINES_PTR(o), LF_SZ_ARR(o), TY_INT)
	    call realloc (LF_WIDTH_PTR(o),  LF_SZ_ARR(o), TY_REAL)
	}

	LF_WAVE(o,n)   = wave
	LF_FLUX(o,n)   = flux
	LF_WT(o,n)     = weight
	LF_NLINES(o,n) = nlines
	LF_WIDTH(o,n)  = width

	return (YES)
end


#-------------------------------------------------------------------------------
#  LF_DEBUG -	Dereference values from the flux object data structure 
#		for viewing.

procedure lf_debug (atom, ion, zone, abund, npts, wave, flux, wt, nlines, width)

#  Calling arguments:
int	atom		# I: atomic number
int	ion		# I: ion
int	zone		# I: ionization zone
real	abund		# I: abundance of ion
int	npts		# I: no. line fluxes
real	wave[ARB]	# I: wavelengths of lines
real	flux[ARB]	# I: fluxes of lines
real	wt[ARB]		# I: summation weights of lines
int	nlines[ARB]	# I: number of lines in multiplet
real	width[ARB]	# I: wavelength interval of multiplet

#  Declarations:
int	i		# generic

begin
	call printf ("Atomic no.: %d;  Ion: %2d; Zone: %d;  Abundance: %8.4g;  No. fluxes: %d\n")
	    call pargi (atom)
	    call pargi (ion)
	    call pargi (zone)
	    call pargr (abund)
	    call pargi (npts)

	if (npts == 0)
	    return

	do i = 1, npts {
	    call printf ("   Wave: %7.1f; Flux: %f; Wt: %5.1f; N_lines: %d; Width: %4.1f\n")
	    call pargr (wave[i])
	    call pargr (flux[i])
	    call pargr (wt[i])
	    call pargi (nlines[i])
	    call pargr (width[i])
	}

	call flush (STDOUT)
end


