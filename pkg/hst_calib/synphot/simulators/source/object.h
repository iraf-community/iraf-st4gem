# OBJECT.H -- Object list structure

define	SZ_OBJSTRUCT	7

define	OBJ_NUMBER	Memi[$1]	# number of objects
define	OBJ_NUMPSF	Memi[$1+1]	# number of psf weights per object
define	OBJ_XPOSARY	Memi[$1+2]	# x position of object center
define	OBJ_YPOSARY	Memi[$1+3]	# y position of object center
define	OBJ_FLUXARY	Memi[$1+4]	# object integrated flux
define	OBJ_PWTARY	Memi[$1+5]	# psf weight array
define	OBJ_SHPARY	Memi[$1+6]	# object shape

define	OBJ_XPOS	Memr[OBJ_XPOSARY($1)+($2)-1]
define	OBJ_YPOS	Memr[OBJ_YPOSARY($1)+($2)-1]
define	OBJ_FLUX	Memr[OBJ_FLUXARY($1)+($2)-1]
define	OBJ_PWEIGHT	Memr[OBJ_PWTARY($1)+OBJ_NUMPSF($1)*(($2)-1)+($3)-1]
define	OBJ_SHAPE	Memi[OBJ_SHPARY($1)+($2)-1]
