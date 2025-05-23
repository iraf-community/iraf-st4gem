# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<error.h>
include	<mach.h>
include	"../atom.h"
include	"../fivel.h"

define	MIN_REAL	1./MAX_REAL

#---------------------------------------------------------------------7 Sep 97--
.help l_ratios.x Mar97 nebular/zones
.ih
NAME
.nf
  l_ratios - Read line intensities & calculate diagnostic line ratios.
  te_ratio - Calculate T_e diagnostic ratio.  
get_dratio - Fetch a diagnostic ratio given the atom, ion, and type.
.fi
.ih
DESCRIPTION
.endhelp
#-------------------------------------------------------------------------------
#  L_RATIOS -	Fetch line intensities from input table and calculate 
#		diagnostic line ratios for subsequent use in zone programs. 

procedure l_ratios (itp, ilr, row)

#  Calling arguments
pointer	itp			# input table descriptor
pointer	ilr			# pointer to line ratio structure
int	row			# input, output table row indexes

#  Local variables
real	c_ext			# nebular extinction constant
int	extn			# index for extinction law
real	S_opt_uv		# optical/UV flux scale factor
real	val[4]			# working variable

#  Functions used:
real	get_flux()		# fetch flux values from input table
real	line_flux()		# fetch line flux & correct for extinction
real	te_ratio()		# Calculate T_e diagnostic ratio


begin
	# Fetch the optical/UV flux scale factor; default value is 1.0. 
	S_opt_uv = get_flux (itp, row, "opt2uv_col")
	if ( IS_INDEF(S_opt_uv) ) 
	    S_opt_uv = 1.0

	# Fetch the extinction law & constant. 
	call get_extinction (itp, row, extn, c_ext)


	## DENSITY RATIOS

	# Fetch ratios for: 	C ii], Ciii], [Ni], N iii], [Oii], O iv], 
	# 			SI iii], [Sii], [CLiii], & [NE iv].
	# If the line ratios are present but zero, set to INDEF. 

	RN_C2(ilr)  = get_flux (itp, row, "r_c2_col")
	RN_C3(ilr)  = get_flux (itp, row, "r_c3_col")
	RN_N1(ilr)  = get_flux (itp, row, "r_n1_col")
	RN_N3(ilr)  = get_flux (itp, row, "r_n3_col")
	RN_O2(ilr)  = get_flux (itp, row, "r_o2_col")
	RN_O4(ilr)  = get_flux (itp, row, "r_o4_col")
	RN_SI3(ilr) = get_flux (itp, row, "r_si3_col")
	RN_S2(ilr)  = get_flux (itp, row, "r_s2_col")
	RN_CL3(ilr) = get_flux (itp, row, "r_cl3_col")
	RN_NE4(ilr) = get_flux (itp, row, "r_ne4_col")

	# Calculate ratio for Ne_[ARiv].
	# Fetch dereddened line fluxes.

	val[1] = line_flux (itp, row, "ar4711_col", 4711., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "ar4740_col", 4740., extn, c_ext, 1.)

	if (!IS_INDEF(val[1]) && !IS_INDEF(val[2])) 
	    RN_AR4(ilr) = val[1] / val[2]
	else 
	    RN_AR4(ilr) = INDEFR


	## TEMPERATURE RATIOS

	# Fetch [N ii] line intensities. Auroral line is val[3]. 
	val[1] = line_flux (itp, row, "n6583_col", 6583., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "n6548_col", 6548., extn, c_ext, 1.)
	val[3] = line_flux (itp, row, "n5755_col", 5755., extn, c_ext, 1.)

	# Compute Te_[N ii].
	RT_N2(ilr) = te_ratio (val[1], val[2], val[3], 2.944)

	# Fetch [O i] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "o6300_col", 6300., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "o6363_col", 6363., extn, c_ext, 1.)
	val[3] = line_flux (itp, row, "o5577_col", 5577., extn, c_ext, 1.)

	# Compute Te_[O i].
	RT_O1(ilr) = te_ratio (val[1], val[2], val[3], 3.035)

	# Fetch [O ii] line intensities. 
	val[1] = line_flux (itp, row, "o3727_col", 3727., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "o7325_col", 7325., extn, c_ext, 1.)

	# Compute Te_[O ii].
	if ( IS_INDEF(val[1]) || IS_INDEF(val[2]) )
	    RT_O2(ilr) = INDEFR
	else 
	    RT_O2(ilr) = val[1] / val[2]

	# Fetch [O iii] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "o5007_col", 5007., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "o4959_col", 4959., extn, c_ext, 1.)
	val[3] = line_flux (itp, row, "o4363_col", 4363., extn, c_ext, 1.)

	# Compute Te_[O iii].
	RT_O3(ilr) = te_ratio (val[1], val[2], val[3], 2.880)

	# Fetch [NE iii] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "ne3869_col", 3869., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "ne3969_col", 3969., extn, c_ext, 1.)
	val[3] = line_flux (itp, row, "ne3342_col", 3342., extn, c_ext, 1.)

	# Compute Te_[NE iii].
	RT_NE3(ilr) = te_ratio (val[1], val[2], val[3], 2.395)

	# Fetch [NE iv] line intensities. 
	val[1] = line_flux (itp, row, "ne2424_col", 2424., extn, c_ext, S_opt_uv)
	val[2] = line_flux (itp, row, "ne1601_col", 1601., extn, c_ext, S_opt_uv)

	# Compute Te_[NE iv].
	if ( IS_INDEF(val[1]) || IS_INDEF(val[2]) )
	    RT_NE4(ilr) = INDEFR
	else 
	    RT_NE4(ilr) = val[1] / val[2]

	# Fetch [NE v] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "ne3426_col", 3426., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "ne3346_col", 3346., extn, c_ext, 1.)
	val[3] = line_flux (itp, row, "ne2975_col", 2975., extn, c_ext, S_opt_uv)

	# Compute Te_[NE v].
	RT_NE5(ilr) = te_ratio (val[1], val[2], val[3], 2.722)

	# Fetch [NA iv] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "na3242_col", 3242., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "na3362_col", 3362., extn, c_ext, 1.)
	val[3] = line_flux (itp, row, "na2805_col", 2805., extn, c_ext, S_opt_uv)

	# Compute Te_[NA iv].
	RT_NA4(ilr) = te_ratio (val[1], val[2], val[3], 2.937)

	# Fetch [NA vi] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "na2871_col", 2871., extn, c_ext, S_opt_uv)
	val[2] = line_flux (itp, row, "na2970_col", 2970., extn, c_ext, S_opt_uv)
	val[3] = line_flux (itp, row, "na2569_col", 2569., extn, c_ext, S_opt_uv)

	# Compute Te_[NA vi].
	RT_NA6(ilr) = te_ratio (val[1], val[2], val[3], 0.33075)

	# Fetch [MG v] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "mg2783_col", 2783., extn, c_ext, S_opt_uv)
	val[2] = line_flux (itp, row, "mg2928_col", 2928., extn, c_ext, S_opt_uv)
	val[3] = line_flux (itp, row, "mg2418_col", 2418., extn, c_ext, S_opt_uv)

	# Compute Te_[MG v].
	RT_MG5(ilr) = te_ratio (val[1], val[2], val[3], 3.346)

	# Fetch [MG vii] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "mg2506_col", 2506., extn, c_ext, S_opt_uv)
	val[2] = line_flux (itp, row, "mg2626_col", 2626., extn, c_ext, S_opt_uv)
	val[3] = line_flux (itp, row, "mg2262_col", 2262., extn, c_ext, S_opt_uv)

	# Compute Te_[MG vii].
	RT_MG7(ilr) = te_ratio (val[1], val[2], val[3], 0.39924)

	# Fetch [AL ii] line intensities. Auroral line is val[2].
	val[1] = line_flux (itp, row, "al2666_col", 2666., extn, c_ext, S_opt_uv)
	val[2] = line_flux (itp, row, "al1671_col", 1671., extn, c_ext, S_opt_uv)

	# Compute Te_[AL ii].
	if ( IS_INDEF(val[1]) || IS_INDEF(val[2]) )
	    RT_AL2(ilr) = INDEFR
	else 
	    RT_AL2(ilr) = val[1] / val[2]

	# Fetch SI iii] line intensities. Auroral line is val[2].
	val[1] = line_flux (itp, row, "si1887_col", 1887., extn, c_ext, S_opt_uv)
	val[2] = line_flux (itp, row, "si1206_col", 1206., extn, c_ext, S_opt_uv)

	# Compute Te_SI iii].
	if ( IS_INDEF(val[1]) || IS_INDEF(val[2]) )
	    RT_SI3(ilr) = INDEFR
	else 
	    RT_SI3(ilr) = val[1] / val[2]

	# Fetch [S ii] line intensities. 
	val[1] = line_flux (itp, row, "s6723_col", 6723., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "s4072_col", 4072., extn, c_ext, 1.)

	# Compute Te_[S ii].
	if ( IS_INDEF(val[1]) || IS_INDEF(val[2]) )
	    RT_S2(ilr) = INDEFR
	else 
	    RT_S2(ilr) = val[1] / val[2]

	# Fetch [S iii] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "s9532_col", 9532., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "s9069_col", 9069., extn, c_ext, 1.)
	val[3] = line_flux (itp, row, "s6312_col", 6312., extn, c_ext, 1.)

	# Compute Te_[S iii].
	RT_S3(ilr) = te_ratio (val[1], val[2], val[3], 2.480)

	# Fetch [CL iv] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "cl8045_col", 8045., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "cl7531_col", 7531., extn, c_ext, 1.)
	val[3] = line_flux (itp, row, "cl5323_col", 5323., extn, c_ext, 1.)

	# Compute Te_[CL iv].
	RT_CL4(ilr) = te_ratio (val[1], val[2], val[3], 2.318)

	# Fetch [AR iii] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "ar7136_col", 7136., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "ar7751_col", 7751., extn, c_ext, 1.)
	val[3] = line_flux (itp, row, "ar5192_col", 5192., extn, c_ext, 1.)

	# Compute Te_[AR iii].
	RT_AR3(ilr) = te_ratio (val[1], val[2], val[3], 4.100)

	# Fetch [AR iv] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "ar4711_col", 4711., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "ar4740_col", 4740., extn, c_ext, 1.)
	val[3] = line_flux (itp, row, "ar2861_col", 2861., extn, c_ext, S_opt_uv)

	# Compute Te_[AR iv]. 
	RT_AR4(ilr) = te_ratio (val[1], val[2], val[3], 1.292)

	# Fetch [AR v] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "ar7006_col", 7006., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "ar6435_col", 6435., extn, c_ext, 1.)
	val[3] = line_flux (itp, row, "ar4626_col", 4626., extn, c_ext, 1.)

	# Compute Te_[AR v].
	RT_AR5(ilr) = te_ratio (val[1], val[2], val[3], 2.143)

	# Fetch [K iv] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "k6102_col", 6102., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "k6796_col", 6796., extn, c_ext, 1.)
	val[3] = line_flux (itp, row, "k4511_col", 4511., extn, c_ext, 1.)

	# Compute Te_[K iv].
	RT_K4(ilr) = te_ratio (val[1], val[2], val[3], 4.597)

	# Fetch [K v] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "k2495_col", 2495., extn, c_ext, S_opt_uv)
	val[2] = line_flux (itp, row, "k2515_col", 2515., extn, c_ext, S_opt_uv)
	val[3] = line_flux (itp, row, "k4123_col", 4123., extn, c_ext, 1.)
	val[4] = line_flux (itp, row, "k4163_col", 4163., extn, c_ext, 1.)

	# Compute Te_[K v].
	if ( IS_INDEF(val[1]) || IS_INDEF(val[2]) )
	    RT_K5(ilr) = INDEFR
	else 
	    RT_K5(ilr) = (val[1] + val[2]) / (val[3] + val[4])

	# Fetch [Ca v] line intensities. Auroral line is val[3].
	val[1] = line_flux (itp, row, "ca5309_col", 5309., extn, c_ext, 1.)
	val[2] = line_flux (itp, row, "ca6087_col", 6087., extn, c_ext, 1.)
	val[3] = line_flux (itp, row, "ca3996_col", 3996., extn, c_ext, 1.)

	# Compute Te_[CA v].
	RT_CA5(ilr) = te_ratio (val[1], val[2], val[3], 5.138)
end


#-------------------------------------------------------------------------------
#  TE_RATIO -	Calculate T_e diagnostic ratio. Return INDEF if auroral or 
#		both nebular lines are INDEF. 

real procedure te_ratio (neb_hi, neb_lo, auroral, neb_ratio)

#  Arguments:
real	neb_hi, neb_lo		# I: strong, weak nebular flux values
real	auroral			# I: auroral flux value
real	neb_ratio		# I: ratio of strong/weak nebular lines
real	diag_ratio		# O: temperature diagnostic

begin

	# If auroral line is INDEF, diagnostic ratio is also INDEF.
	if ( IS_INDEF(auroral) || auroral < EPSILONR)
	    return (INDEFR)

	# If both nebular lines are INDEF, diagnostic ratio is also INDEF.
	if ( (IS_INDEF(neb_hi) && IS_INDEF(neb_lo)) ) 
	    return (INDEFR)

	if ( !IS_INDEF(neb_ratio) ) {
	    if ( IS_INDEF(neb_lo) )
	    	neb_lo = neb_hi / neb_ratio 

	    else if ( IS_INDEF(neb_hi) )
	    	neb_hi = neb_lo * neb_ratio 
	}

	diag_ratio = (neb_hi + neb_lo) / auroral

	return (diag_ratio)
end


#-------------------------------------------------------------------------------
#  GET_DRATIO -	Fetch a diagnostic ratio given the atom, ion, and type.

real procedure get_dratio (ilr, atom, ion, diag_type)

#  Calling arguments:
pointer	ilr		# I: diagnostic line ratio structure
int	atom		# I: atomic number
int	ion		# I: spectrum no.
int	diag_type	# I: diagnostic type
real	ratio		# O: diagnostic line ratio

begin
	ratio = INDEFR

	# Set diagnostic line ratio
	switch (atom) {
	case CARBON:
	    if (ion == 1 && diag_type == DENSITY) 
		ratio = RN_C2(ilr)

	    else if (ion == 2 && diag_type == DENSITY)
		ratio = RN_C3(ilr)

	case NITROGEN:
	    if (ion == 0 && diag_type == DENSITY) 
		ratio = RN_N1(ilr)

	    else if (ion == 1 && diag_type == TEMPERATURE)
		ratio = RT_N2(ilr)

	    else if (ion == 2 && diag_type == DENSITY)
		ratio = RN_N3(ilr)

	case OXYGEN:
	    if (ion == 0 && diag_type == TEMPERATURE) 
		ratio = RT_O1(ilr)

	    else if (ion == 1 && diag_type == DENSITY)
		ratio = RN_O2(ilr)

	    else if (ion == 1 && diag_type == TEMPERATURE)
		ratio = RT_O2(ilr)

	    else if (ion == 2 && diag_type == TEMPERATURE)
		ratio = RT_O3(ilr)

	    else if (ion == 3 && diag_type == DENSITY)
		ratio = RN_O4(ilr)

	case NEON:
	    if (ion == 2 && diag_type == TEMPERATURE) 
		ratio = RT_NE3(ilr)

	    else if (ion == 3 && diag_type == DENSITY)
		ratio = RN_NE4(ilr)

	    else if (ion == 3 && diag_type == TEMPERATURE)
		ratio = RT_NE4(ilr)

	    else if (ion == 4 && diag_type == TEMPERATURE)
		ratio = RT_NE5(ilr)

	case SODIUM:
	    if (ion == 3 && diag_type == TEMPERATURE) 
		ratio = RT_NA4(ilr)

	    else if (ion == 5 && diag_type == TEMPERATURE)
		ratio = RT_NA6(ilr)

	case MAGNESIUM:
	    if (ion == 4 && diag_type == TEMPERATURE) 
		ratio = RT_MG5(ilr)

	    else if (ion == 6 && diag_type == TEMPERATURE)
		ratio = RT_MG7(ilr)

	case ALUMINUM:
	    if (ion == 1 && diag_type == DENSITY) 
		ratio = RN_AL2(ilr)

	    else if (ion == 1 && diag_type == TEMPERATURE)
		ratio = RT_AL2(ilr)

	case SILICON:
	    if (ion == 2 && diag_type == DENSITY) 
		ratio = RN_SI3(ilr)

	    else if (ion == 2 && diag_type == TEMPERATURE)
		ratio = RT_SI3(ilr)

	case SULFUR:
	    if (ion == 1 && diag_type == DENSITY) 
		ratio = RN_S2(ilr)

	    else if (ion == 1 && diag_type == TEMPERATURE)
		ratio = RT_S2(ilr)

	    else if (ion == 2 && diag_type == TEMPERATURE)
		ratio = RT_S3(ilr)

	case CHLORINE:
	    if (ion == 2 && diag_type == DENSITY) 
		ratio = RN_CL3(ilr)

	    else if (ion == 3 && diag_type == TEMPERATURE)
		ratio = RT_CL4(ilr)

	case ARGON:
	    if (ion == 2 && diag_type == TEMPERATURE) 
		ratio = RT_AR3(ilr)

	    else if (ion == 3 && diag_type == DENSITY)
		ratio = RN_AR4(ilr)

	    else if (ion == 3 && diag_type == TEMPERATURE)
		ratio = RT_AR4(ilr)

	    else if (ion == 4 && diag_type == TEMPERATURE)
		ratio = RT_AR5(ilr)

	case POTASSIUM:
	    if (ion == 3 && diag_type == TEMPERATURE) 
		ratio = RT_K4(ilr)

	    else if (ion == 4 && diag_type == DENSITY)
		ratio = RN_K5(ilr)

	case CALCIUM:
	    if (ion == 4 && diag_type == TEMPERATURE) 
		ratio = RT_CA5(ilr)

	default:
	    ratio = INDEFR
	}

	return (ratio)
end


