# GRFLIST.H -- Definition of the graph list structure

define	LEN_GRFSTRUCT	7

define	GRF_SIZE	Memi[$1]
define	GRF_CMPLEN	Memi[$1+1]
define	GRF_KEYLEN	Memi[$1+2]
define	GRF_INARRAY	Memi[$1+3]
define	GRF_OUTARRAY	Memi[$1+4]
define	GRF_CMPARRAY	Memi[$1+5]
define	GRF_KEYARRAY	Memi[$1+6]
define	GRF_INNODE	Memi[GRF_INARRAY($1)+($2)-1]
define	GRF_OUTNODE	Memi[GRF_OUTARRAY($1)+($2)-1]
define	GRF_COMPNAME	Memc[GRF_CMPARRAY($1)+GRF_CMPLEN($1)*(($2)-1)]
define	GRF_KEYNAME	Memc[GRF_KEYARRAY($1)+GRF_KEYLEN($1)*(($2)-1)]
