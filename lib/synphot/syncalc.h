# SYNCALC -- Constants use by synphot expression calculator

# Operation tokens, copied from syncompile.x

define	Y_LPAR		257
define	Y_RPAR		258
define	Y_COMMA		259
define	Y_WRONG		260
define	Y_DONE		261
define	Y_STR		262
define	Y_FILE		263
define	Y_NUM		264
define	Y_VAR		265
define	Y_FUNC		266
define	Y_ADD		267
define	Y_SUB		268
define	Y_MUL		269
define	Y_DIV		270
define	Y_NEG		271

# Function tokens stored in pseudocode

define	FN_BAND		1 
define	FN_BB 		2
define	FN_BOX 		3
define	FN_CAT		4
define	FN_EBMV 	5
define	FN_EBMVX 	6
define  FN_EM		7
define	FN_GAUSS 	8
define	FN_GRID 	9
define	FN_HI 		10
define	FN_ICAT		11
define	FN_LGAUSS 	12
define	FN_PL 		13
define	FN_POLY 	14
define	FN_RN 		15
define	FN_SPEC		16
define	FN_THRU		17
define	FN_TILT 	18
define	FN_UNIT		19
define	FN_Z		20

define	MAX_FUNC	FN_Z

define	CAT_TEMPLATE	"crrefer$grid/*/catalog"

# Function names

define	FUNCNAMES	"band,bb,box,cat,ebmv,ebmvx,em,gauss,grid,hi,icat,\
lgauss,pl,poly,rn,spec,thru,tilt,unit,z"
