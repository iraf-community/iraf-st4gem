# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<mach.h>
include	<tbset.h>
include	"../fivel.h"
include	"../neberr.h" 

define	DEBUG		false

#---------------------------------------------------------------------7 Sep 97--
.help f_ratios.x Jul97 nebular/abund
.ih
NAME
 f_ratios - Compute unreddened line fluxes and normalize to H-beta 
line_flux - Fetch input fluxes from input table & correct I.S. extinction
.endhelp
#-------------------------------------------------------------------------------
#  F_RATIOS -	Fetch line intensities from input table and calculate 
#		emission line flux ratios for subsequent use in ionic 
#		abundance programs.  Note: fluxes are normalized to the scale 
#		I(4861) = 1.0.	

procedure f_ratios2 (itp, row, lfl)

#  Calling arguments
pointer	itp			# I: input table descriptor
int	row			# I: input, output table row indexes
pointer	lfl			# I: pointer to flux object list

#  Declarations
pointer	atl_get_at()		# get LF for atom, ion
int	atom			# atomic number
pointer	atp			# abundance table descriptor
real	c_ext			# nebular extinction constant
char	colname[SZ_COLNAME]	# table column names
pointer	colptr[8]		# table column pointers
int	extn			# index for extinction law
real	flux			# scaled, unreddened emission line flux
pointer	lf			# flux data object
int	lf_add_line()		# add line to flux data object
real	get_flux()		# fetch flux values from input table
real	hbeta			# H-beta (4861 Ang) flux
int	i			# generic
int	ion			# ion (spectrum) number
real	line_flux()		# fetch line flux & correct for extinction
real	norm			# H-beta normalization factor
int	nlines			# no. lines in multiplet
int	n_rows			# no. abundance table rows
real	S_opt_uv		# optical/UV flux scale factor
real	wave			# nominal wavelength of observed multiplet
real	width			# wavelength interval
real	wt			# relative weight of flux 
int	status			# return status
int	zone			# ionization zone

begin

	# Fetch the optical/UV flux scale factor; default value is 1.0. 
	S_opt_uv = get_flux (itp, row, "opt2uv_col")
	if ( IS_INDEF(S_opt_uv) ) 
	    S_opt_uv = 1.0

	# Fetch the extinction law & constant. 
	call get_extinction (itp, row, extn, c_ext)

	# Insist that the H-beta flux be present and non-zero.
	hbeta = get_flux (itp, row, "h4861_col")
	if ( IS_INDEFR (hbeta) ) 
	    call error (NO_HBETA, "Missing H-beta flux")
	else
	    hbeta = 1. / hbeta 

	# Open abundance reference file.
	call open_abund (atp, colptr, n_rows)

	# Fetch values from the abundance table, renormalize to H-beta=100, and 
	# correct for interstellar extinction. Note that the normalization 
	# must include the optical/UV scale factor for wavelengths blueward 
	# of 3000 Ang. 

	do i = 1, n_rows {
	    call get_abund_line (atp, colptr, i, atom, ion, zone, wave, wt, 
				nlines, width, colname)
	    if (wave < 3000.)
		norm = hbeta * S_opt_uv
	    else
		norm = hbeta

	    flux = line_flux (itp, row, colname, wave, extn, c_ext, norm)
	    if (!IS_INDEFR (flux)) {

	    	# Fetch the LF for this atom, ion; or add one if not available.
	    	lf = atl_get_at (lfl, atom, ion)
	    	if (lf == NULL) {
	    	    call lf_alloc (atom, ion, lf)
	    	    call lf_set_zone (lf, zone)
		    call atl_add (lfl, lf)
	    	}

	    	status = lf_add_line (lf, wave, flux, wt, nlines, width)
	    }
	}

	if (DEBUG)
	    call lfl_debug (lfl)

	call tbtclo (atp)
end


#-------------------------------------------------------------------------------
#  LINE_FLUX - Fetch input fluxes from input table & correct I.S. extinction. 

real procedure line_flux (itp, row, c_name, wave, ext_fn, c_ext, norm)

#  Arguments:
pointer	itp		# I: input table descriptor
int	row		# I: input table row
char	c_name[ARB]	# I: user column name
real	wave		# I: wavelength of flux
int	ext_fn		# I: index of I.S. extinction function
real	c_ext		# I: extinction constant
real	norm		# I: flux normalization
real	val		# O: returned value

#  Declarations:
real	flux		# generic
real	get_flux()	# fetch flux values from input table

begin
	# Fetch flux from the input table.
	flux = get_flux (itp, row, c_name)

	# Correct flux for interstellar extinction. 
	if (IS_INDEFR (flux))
	    val = INDEFR

	else if (abs (c_ext) < EPSILONR)
	    val = flux * norm

	else {
	    call dered (flux, wave, val, 1, ext_fn, c_ext) 
	    val = val * norm
	}

	return (val)
end


