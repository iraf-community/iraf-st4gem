# flerr.h	--	Include file for error messages		17-Jan-97

#  Status codes for temperature/density calculation:
#define	OK			0
define	NOT_CONVERGENT		1
define	ARITHMETIC_ERROR	2
define	LR_OUT_OF_BOUNDS	3
define	LR_NEAR_LIMIT		4
define	ION_INVALID		5
define	DIAG_INVALID		6
define	DIAG_FAILED		7
define	ABUND_FAILED		8

#  File i/o error codes: 
define	REF_TABLE_ACCESS	11	# Unable to access reference data
define	TABLE_ACCESS		12	# Unable to access user data
define	NO_HBETA		13	# Missing H-beta flux

#  Input errors:
define	INVLD_ATOM_ION		21	# Unsupported atom/ion combination
define	BAD_REF_DATA		22	# Bad reference data for atom/ion
define	USER_FUNCTION		23	# Error in user function
define	BAD_OBJ_ID		24	# Missing or bad object ID/Region
define	BAD_NT			25	# Error in N_e or T_e input
define	BAD_N_E			26	# Missing N_e input
define	BAD_T_E			27	# Missing T_e input
define	BAD_RATIO		28	# Missing/invalid input line ratio
define	BAD_EXTN_LAW		29	# Invalid I.S. extinction law
define	BAD_EXTN_WAVE		30	# Invalid wavelength extinction law

define	BAD_PLOT_RESOL		31	# Bad plot resolution
define	INVLD_PLOT		32	# No valid points to plot

define	ATL_NOCURRENT		41	# no current AT object
define	ATL_NOSUCHID		42	# no AT object with that ID
