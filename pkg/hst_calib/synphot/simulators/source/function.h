# FUNCTION.H -- Structure containing the function parameters

define	LEN_FUNSTRUCT	3

define	FUN_NPAR	Memi[$1]	# number of parameters to function
define	FUN_TYPARY	Memi[$1+1]	# array of parameter types
define	FUN_PARARY	Memi[$1+2]	# array of ptrs to parameters

define	FUN_TYPE	Memi[FUN_TYPARY($1)+$2-1]
define	FUN_PARAM	Memi[FUN_PARARY($1)+$2-1]

define	FUN_STR		Memc[FUN_PARAM($1,$2)]
define	FUN_NUM		Memr[FUN_PARAM($1,$2)]

