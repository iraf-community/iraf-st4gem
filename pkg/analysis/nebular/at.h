include <iraf.h>
#	at.h	Atomic data object definition.			22-May-97

#  Ground-state electron configurations:
define	GS_CONFIG	"|S1|S2|P1|P2|P3|P4|P5|P6|"
define	GS_S1		1
define	GS_S2		2
define	GS_P1		3
define	GS_P2		4
define	GS_P3		5
define	GS_P4		6
define	GS_P5		7
define	GS_P6		8

#  Diagnostic types:
define	TEMDEN		"|density|temperature|emissivity|intensity|"
define	DENSITY		1
define	TEMPERATURE	2
define	EMISSIVITY	3
define	INTENSITY	4

define	AT_DATA		1
define	AT_COLLISION	2

#-------------------------------------------------------------------------------
#  Data structure for T_e independent atomic parameters.
define	AT_ATOM		Memi[$1]	  # atomic no. for this species.
define	AT_ION		Memi[$1+1]	  # ion stage for this species.
define	AT_NLVL		Memi[$1+2]	  # no. energy levels for this ion
define	AT_NTRANS	Memi[$1+3]	  # no. transitions for this ion
define	AT_FLUX_LIST	Memi[$1+4]	  # ptr to list of line fluxes
define	AT_WT		Memi[$1+5]	  # statistical weights
define	AT_E_TR		Memi[$1+6]	  # transition energies
define	AT_L_TR		Memi[$1+7]	  # level transition matrix for this ion
define	AT_RAD_TR	Memi[$1+8]	  # radiative transition probabilities
define	AT_GSCONFIG	Memi[$1+9]	  # ground-state configuration

#  T_e dependent atomic parameters.
define	AT_TE		Memr[P2R($1+10)]  # electron temperature
define	AT_TE_MIN	Memr[P2R($1+11)]  # minimum valid T_e
define	AT_TE_MAX	Memr[P2R($1+12)]  # maximum valid T_e
define	AT_LOG_TE	Memb[$1+13]	  # flag for T_e tabulated as log
define	AT_TY_FUNC	Memi[$1+14]	  # type of functional fit for cross sections
define	AT_CV		Memi[$1+15]	  # curve fit data structure
define	AT_COLL		Memi[$1+16]	  # collision strengths
define	AT_COLL_TR	Memi[$1+17]	  # collisional transition probabilities
define	AT_NCRIT	Memi[$1+18]	  # critical densities
define	AT_POP		Memi[$1+19]	  # level populations
define	AT_EMISS	Memi[$1+20]	  # line emissivities
define	LEN_AT		24		  # size of structure

#  Memory management.
define	COLL		Memd[AT_COLL($1)]
define	COLL_TR		Memd[AT_COLL_TR($1)]
define	CURVE		Memi[AT_CV($1)+$2-1]
define	EMISS		Memd[AT_EMISS($1)]
define	E_TRANS		Memd[AT_E_TR($1)]
define	L_TRANS		Memd[AT_L_TR($1)]
define	N_CRIT		Memd[AT_NCRIT($1)+$2-1]
define	POP		Memd[AT_POP($1)]
define	RAD_TR		Memd[AT_RAD_TR($1)]
define	WEIGHT		Memi[AT_WT($1)]

#-------------------------------------------------------------------------------
# Structure to define a list of atomic data objects.  The atl_* routines operate
# on lists of atomic data objects.
define	ATL_A_PTR	Memi[P2I($1+0)]	# Pointer to array of atomic data objects
define	ATL_A_SZ	Memi[P2I($1+1)]	# Size of the array
define	ATL_N		Memi[P2I($1+2)]	# Number of atomic data objects in array
define	ATL_SZ		8

# Access the array.
define	ATL_A			Memi[ATL_A_PTR($1)+$2-1]

# How much the array should grow when adding atomic data objects.
define	ATL_GROW	8
#-------------------------------------------------------------------------------
