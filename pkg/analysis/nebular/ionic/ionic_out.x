# Copyright restrictions apply - see st4gem$copyright.stsdas

include <mach.h>
include	"../at.h"
include	"../fivel.h"

define	DEBUG	false

#--------------------------------------------------------------------13 Feb 98--
.help ionic_out.x Feb97 nebular/fivel
.ih
NAME
.nf
      ion_out - Write IONIC output to STDOUT. 
format_labels -	Format label strings for IONIC output.
.fi
.endhelp
#-------------------------------------------------------------------------------
#  ION_OUT -	Write IONIC output to STDOUT. 

procedure ion_out (at, flux, wave, dens, temp, result, width, verb)

#  Calling arguments:
pointer	at		# I: atomic data structure
real	flux		# I: diagnostic line ratio
real	wave		# I: wavelength of emission line
real	dens		# I: electron density
real	temp		# I: electron temperature 
real	result		# I: result of abundance calculation
real	width		# I: wave interval for input flux
bool	verb		# I: print level populations & critical densities?

#  Declarations:
double	atg_emiss()	# compute the emissivity for an atomic transition
double	atg_wave()	# compute the wavelength for an atomic transition
real	emisv()		# compute emissivity from one or more lines
real	hb_emiss()	# compute emissivity for H-beta
int	i, j		# loop indexes
real	jhb		# vol emissivity of H-beta
char	symbol[SZ_FNAME] # matched keyword
double	min_double	# Min useful REAL value
int	n_lvls		# no. atomic levels
pointer	s_emis		# emissivity output
pointer	s_label		# label output
pointer	s_wave		# wavelength output
pointer	sp		# top of stack memory
pointer	xl		# calculated wavelengths
pointer	xi		# calculated intensities

#  Memory management
define	LVL_POP		Memd[AT_POP(at)+($1-1)]
define	S_emis		Memc[s_emis]
define	S_label		Memc[s_label]
define	S_wave		Memc[s_wave]
define	Xl		Memd[xl+$1-1]
define	Xi		Memd[xi+$1-1]

begin
	min_double = log10 (1. / MAX_DOUBLE)
	n_lvls = AT_NLVL(at)

	call smark (sp)
	call salloc (s_emis,  SZ_LINE, TY_CHAR)
	call salloc (s_label, SZ_LINE, TY_CHAR)
	call salloc (s_wave,  SZ_LINE, TY_CHAR)
	call salloc (xl, n_lvls, TY_DOUBLE)
	call salloc (xi, n_lvls, TY_DOUBLE)

	# Fetch atomic symbol. 
	call get_atomic_symbol (AT_ATOM(at), symbol, SZ_FNAME)

	if (DEBUG)
	    call at_debug (at)

	# Write out banner.
	call printf ("\n# Volume Emissivities for: %s^%d+ \n")
	    call pargstr (symbol)
	    call pargi (AT_ION(at))
	call printf ("   T_e: %8.1f; N_e: %8.4e\n")
	    call pargr (temp)
	    call pargr (dens)

	# Write level populations and critical densities. 
	if (verb) {
	    call printf ("\n# Level Populations - Critical Densities\n")
	    call printf ("\n Level 1: %8.2e \n")
	    	call pargd (LVL_POP(1))

	    do i = 2, n_lvls {
	    	call printf (" Level %d: %8.4e    %8.4e\n")
		    call pargi (i)
	    	    call pargd (LVL_POP(i))
	    	    call pargd (N_CRIT(at,i))
	    }
	}

	do i = 2, n_lvls {

	    # Derive wavelengths & intensities. 
	    do j = 1, i-1 {
		Xl(j) = atg_wave (at, i, j)
		Xi(j) = atg_emiss (at, i, j, dens)
	    }

	    # Write emissivities.
	    # Build up the format strings.
	    call format_labels (S_emis, S_label, S_wave, SZ_LINE, i-1)

	    if (i == 2) {
		call strcat ("# Wavelength\n", S_wave, SZ_LINE)
		call strcat ("# Upper->Lower Level\n", S_label, SZ_LINE)
		call strcat ("# Volume Emissivity\n", S_emis, SZ_LINE)
	    } else {
		call strcat ("\n", S_wave, SZ_LINE)
		call strcat ("\n", S_label, SZ_LINE)
		call strcat ("\n", S_emis, SZ_LINE)
	    }

	    call printf (S_wave)
	    do j = 1, i-1 
		call pargd (Xl(j))

	    call printf (S_label)
	    do j = 1, i-1 {
		call pargi (i)
		call pargi (j)
	    }

	    call printf (S_emis)
	    do j = 1, i-1 
		call pargd (Xi(j))
	}

	# Write H-beta emissivity & log10 of "X"
	jhb = hb_emiss (temp)
	call printf (
	    "\n# H-beta Volume Emissivity:\n %9.4e N(H+) * N(e-)  (erg/s)\n")
		call pargr (jhb)
	call printf ("\n Log10(x) = %9.4e\n")
	    call pargr ( log10 (dens/sqrt(temp)) )

	# Derive relative ionic abundance if "wave" and "flux" specified. 
	result = emisv (flux, wave, width, dens, E_TRANS(at), EMISS(at), n_lvls)
	if (IS_INDEFR (result)) 
	    flux = INDEFR

	else if (result > 0.) 
	    result = flux * jhb / result

	else
	    flux = 0.

	# Print abundance determination.
	if (!IS_INDEF (flux) && !IS_INDEF (wave)) {
	    call clputr ("result", result)
	    call printf ("\n Ionic Abundance: N(%s^%d+) / N(H+) = %9.4e \n")
	    	call pargstr (symbol)
	    	call pargi (AT_ION(at))
	    	call pargr (result)
	}

	call sfree (sp)
end


#-------------------------------------------------------------------------------
#  FORMAT_LABELS -	Format label strings for IONIC output. 

procedure format_labels (s_emis, s_label, s_wave, sz_label, n_report)

#  Arguments:
char	s_emis[ARB]	# I: emissivity label
char	s_label[ARB]	# I: transition label
char	s_wave[ARB]	# I: wavelength label
int	sz_label	# I: size of input label strings
int	n_report	# I: no. levels to report

#  Local variables:
int	i		# generic

begin
	call strcpy ("\n", s_wave,  sz_label)
	call strcpy ("", s_label, sz_label)
	call strcpy ("", s_emis,  sz_label)

	do i = 1, n_report {
	    call strcat ("%9.2f   ", s_wave, sz_label)
	    call strcat ("  (%d-->%d)   ", s_label, sz_label)
	    call strcat ("%9.4e   ", s_emis, sz_label)
	}

end


