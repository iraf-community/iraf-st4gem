include <iraf.h>
#  flux.h	Flux object definition.				9-Jul-97

#-------------------------------------------------------------------------------
#  Data structure of emission line fluxes.
define	LF_ATOM		Memi[$1]	# atomic no. for this species.
define	LF_ION		Memi[$1+1]	# ion stage for this species.
define	LF_ZONE		Memi[$1+2]	# ionization zone 
define	LF_ABUND	Memr[P2R($1+3)]	# ionic abundance
define	LF_N		Memi[$1+4]	# no. lines stored in structure
define	LF_SZ_ARR	Memi[$1+5]	# size of obs. line arrays
define	LF_GROW		4		# incremental array growth size

define	LF_WAVE_PTR	Memi[$1+6]	# array of nominal wavelengths
define	LF_FLUX_PTR	Memi[$1+7]	# array of fluxes @wave
define	LF_WT_PTR	Memi[$1+8]	# array of weights in abundance avg
define	LF_NLINES_PTR	Memi[$1+9]	# no. lines in each multiplet
define	LF_WIDTH_PTR	Memi[$1+10]	# wavelength interval of lines in mult.

define	LEN_LF		12		# size of structure

# Memory management
define	LF_WAVE		Memr[LF_WAVE_PTR($1)+$2-1]
define	LF_FLUX		Memr[LF_FLUX_PTR($1)+$2-1]
define	LF_WT		Memr[LF_WT_PTR($1)+$2-1]
define	LF_NLINES	Memi[LF_NLINES_PTR($1)+$2-1]
define	LF_WIDTH	Memr[LF_WIDTH_PTR($1)+$2-1]

#-------------------------------------------------------------------------------
# Structure to define a list of flux objects.  The lfl_* routines operate
# on lists of flux objects.
define	LFL_A_PTR	Memi[$1]	# Pointer to array of flux data objects
define	LFL_A_SZ	Memi[$1+1]	# Size of the array
define	LFL_N		Memi[$1+2]	# Number of flux data objects in array
define	LEN_LFL		4

# Access the array.
define	LFL_A		Memi[LFL_A_PTR($1)+$2-1]

# How much the array should grow when adding flux data objects.
define	LFL_GROW	8
#-------------------------------------------------------------------------------
