include <iraf.h>
#  fivel.h --	Include file for line ratio structure		25-Aug-2008

define	N_LINES		 4

#  Interstellar reddening correction function
define	EXTN_MODEL	"|gal|ccm|lmc|smc|jbk|"
define	GAL		1
define	CCM		2
define	LMC		3
define	SMC		4
define	JBK		5

#  Structure for line ratios
define	RN_C2		Memr[P2R($1+0)]	        # Density ratio for   C ii]
define	RN_C3		Memr[P2R($1+1)]	        # Density ratio for   C iii]
define	RN_N1		Memr[P2R($1+2)]	        # Density ratio for  [N i]
define	RN_N3		Memr[P2R($1+3)]	        # Density ratio for   N iii]
define	RN_O2		Memr[P2R($1+4)]	        # Density ratio for  [O ii]
define	RN_O4		Memr[P2R($1+5)]	        # Density ratio for   O iv]
define	RN_NE2		Memr[P2R($1+6)]	        # Density ratio for [Ne ii]
define	RN_NE4		Memr[P2R($1+7)]	        # Density ratio for [Ne iv]
define	RN_AL2		Memr[P2R($1+8)]	        # Density ratio for [Al ii]
define	RN_SI3		Memr[P2R($1+9)]	        # Density ratio for  Si iii]
define	RN_S2		Memr[P2R($1+10)]	# Density ratio for  [S ii]
define	RN_S4		Memr[P2R($1+11)]	# Density ratio for  [S iv]
define	RN_CL3		Memr[P2R($1+12)]	# Density ratio for [Cl iii]
define	RN_AR2		Memr[P2R($1+13)]	# Density ratio for [Ar ii]
define	RN_AR4		Memr[P2R($1+14)]	# Density ratio for [Ar iv]
define	RN_K5		Memr[P2R($1+15)]	# Density ratio for  [K v]

define	RT_N2		Memr[P2R($1+21)]	# Temperature ratio for  [N ii]
define	RT_O1		Memr[P2R($1+22)]	# Temperature ratio for  [O i]
define	RT_O2		Memr[P2R($1+23)]	# Temperature ratio for  [O ii]
define	RT_O3		Memr[P2R($1+24)]	# Temperature ratio for  [O iii]
define	RT_NE3		Memr[P2R($1+25)]	# Temperature ratio for [Ne iii]
define	RT_NE4		Memr[P2R($1+26)]	# Temperature ratio for [Ne iv]
define	RT_NE5		Memr[P2R($1+27)]	# Temperature ratio for [Ne v]
define	RT_NA4		Memr[P2R($1+28)]	# Temperature ratio for [Na iv]
define	RT_NA6		Memr[P2R($1+29)]	# Temperature ratio for [Na vi]
define	RT_MG5		Memr[P2R($1+30)]	# Temperature ratio for [Mg v]
define	RT_MG7		Memr[P2R($1+31)]	# Temperature ratio for [Mg vii]
define	RT_AL2		Memr[P2R($1+32)]	# Temperature ratio for [Al ii]
define	RT_SI3		Memr[P2R($1+33)]	# Temperature ratio for  Si iii]
define	RT_S2		Memr[P2R($1+34)]	# Temperature ratio for  [S ii]
define	RT_S3		Memr[P2R($1+35)]	# Temperature ratio for  [S iii]
define	RT_CL4		Memr[P2R($1+36)]	# Temperature ratio for [Cl iv]
define	RT_AR3		Memr[P2R($1+37)]	# Temperature ratio for [Ar iii]
define	RT_AR4		Memr[P2R($1+38)]	# Temperature ratio for [Ar iv]
define	RT_AR5		Memr[P2R($1+39)]	# Temperature ratio for [Ar v]
define	RT_K4		Memr[P2R($1+40)]	# Temperature ratio for  [K iv]
define	RT_K5		Memr[P2R($1+41)]	# Temperature ratio for  [K v]
define	RT_CA5		Memr[P2R($1+42)]	# Temperature ratio for [Ca v]

define	LEN_ILR		44		# No. input line ratios

