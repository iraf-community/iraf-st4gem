# FOC.H
# FOC specific definitions

# Image specific defines:
# Include for FOC coordinates, etc

define	FOC_COORD_OK		1
define	FOC_NO_OVERLAP		-1
define	FOC_PART_OVERLAP	2
define	FOC_OK			0
define	FOC_ERROR		99999

define	SZ_FOC_COORD		20

# FOC_SAMPOFF and FOC_LINEOFF in this struct are the values of sampoff
# and lineoff as gotten directly from the image header.  The value of
# FOC_SAMPBEG, however, is "flipped" so it refers to the mirror image
# of the photocathode relative to FOC_SAMPOFF.  FOC_LINEBEG is just
# FOC_LINEOFF + 1.

define	FOC_IMPT		Memi[$1]
define	FOC_SOFF1		Memi[$1+1]
define	FOC_SOFF2		Memi[$1+2]
define	FOC_NPIX1		Memi[$1+3]
define	FOC_NPIX2		Memi[$1+4]
define	FOC_INC1		Memi[$1+5]
define	FOC_INC2		Memi[$1+6]
define	FOC_LINEBEG		Memi[$1+7]
define	FOC_SAMPBEG		Memi[$1+8]
define	FOC_ACCMODE		Memi[$1+9]
define	FOC_COORD_INIT		Memb[$1+10]
define	FOC_PC_COORDS		Memb[$1+11]
define	FOC_SAMPOFF		Memr[P2R($1+12)]
define	FOC_LINEOFF		Memr[P2R($1+13)]

# end of foc_coordinate structure

# definitions for reading text files

define	INIT_COORD_SIZE	2000
define	INC_COORD_SIZE	1000
define	MAX_MODELS	50
define	MAX_MODEL_DIM	50
define	SZ_MODEL_INFO	4
define	MODEL_NUM	Meml[$1]
define	MODEL_IN_RES	Memb[$1+SZ_MODEL_INFO*$2-3]
define	MODEL_PTR	Meml[$1+SZ_MODEL_INFO*$2-2]
define	MODEL_NAXIS1	Meml[$1+SZ_MODEL_INFO*$2-1]
define	MODEL_NAXIS2	Meml[$1+SZ_MODEL_INFO*$2  ]

# definition for adding group parameter

define	FOC_GROUP_PAR	1

# Reseau file specific defines:
#
# length of character columns

define	SZ_R_PATTERN		63	#Length of pattern for entry name
define	SZ_R_ENTRY		15	#Length of character columns
define	SZ_R_TRACKING		15
define	SZ_R_ENTNAME		63
define	SZ_R_INFOCOLS		15
define	SZ_R_CMAXSIZE		23
define	SZ_RES_MESS		131
define	SZ_R_TIME		22

define	RES_P_ENTRY		1
define	RES_P_TRACKING		2
define	RES_P_DATE		3
define	RES_P_NROWS		4
define	RES_P_NCOLS		5
define	RES_P_SFIELD1		6
define	RES_P_SFIELD2		7
define	RES_P_SFIELD3		8
define	RES_P_SFIELD4		9

define	NUM_R_RINFOCOLS		4
define	NUM_R_INFOCOLS		9
define	NUM_R_EXTRAROWS		10
define	NUM_R_EXTRACOLS		10

define	RES_ER_OVEND		88100
define	RES_ER_NOENTRIES	88200
define	RES_ER_NOTRESFILE	88300
define	RES_ER_WRONGSIZE	88400
define	RES_F_USEDWANTED	 1
define	RES_F_FIRST		-1
define	RES_F_LAST		-2
define	RES_F_USEDANOTHER	-3
define	RES_F_NONEXTENTRY	-1000

# The reseau files structure:

define	LEN_RES_INFO		5+NUM_R_INFOCOLS
define	RES_FILE		Meml[($1)]
define	RES_NROWS		Meml[($1)+1]
define	RES_NCOLS		Meml[($1)+2]
define	RES_ALL_NAMES		Meml[($1)+3]
define	RES_EXTENT		Meml[($1)+4]
define	RES_COLPTR		Meml[($1)+5+$2-1]

define	LEN_RES_TINDEX		1		# additional length for index

# The reseau entry info columns structure

define	LEN_ENTRY_INFO		80

define	RES_ENT_NROWS		Memi[($1)]
define	RES_ENT_NCOLS		Memi[($1) + 1]
define	RES_ENT_DATE		Meml[P2L($1 + 2)]
define	RES_ENT_SFIELD1		Memr[P2R(($1) + 3)]
define	RES_ENT_SFIELD2		Memr[P2R(($1) + 4)]
define	RES_ENT_SFIELD3		Memr[P2R(($1) + 5)]
define	RES_ENT_SFIELD4		Memr[P2R(($1) + 6)]
define	RES_ENT_ENTRY		Memc[P2C(($1)+7)]
define	RES_ENT_TRACKING	Memc[P2C(($1)+15)]
define	RES_XPT			Meml[$1 + 23]
define	RES_YPT			Meml[$1 + 24]
define	RES_FPT			Meml[$1 + 25]
define	RES_DATA_TYPE		Memi[$1 + 26]

define	REF_DATA	Meml[$1]
define	RES_DATA	Meml[$1+1]
define	RES_TEMPL_COUNT	Memi[$1]		# keeps pointer for template use
# Type for reseau data type

define	RES_TYPE		real
define	RES_TYPE_SZ		SZ_REAL
define	RES_TYPE_TYPE		TY_REAL
define	RES_TYPE_MEM		Memr
define	RES_ROWSPERENTRY	1		# 1 row per grid
define	RES_COLSPERRES		2		# x and y col each row
define	RES_USERPARS		10
define	RES_INDEF		INDEFR
