#ADJLIST.H -- Definition of the instrument adjacency list structure

define	LEN_ADJSTRUCT	5

define	ADJ_SIZE	Memi[$1]
define	ADJ_NODPTR	Memi[$1+1]
define	ADJ_NUMPTR	Memi[$1+2]
define	ADJ_NAMPTR	Memi[$1+3]
define	ADJ_BACKPTR	Memi[$1+4]
define	ADJ_NODARY	Memi[ADJ_NODPTR($1)+($2)-1]
define	ADJ_NAMARY	Memi[ADJ_NAMPTR($1)+($2)-1]
define	ADJ_BACKARY	Memi[ADJ_BACKPTR($1)+($2)-1]
define	ADJ_NODE	Memi[ADJ_NODARY($1,$2)+($3)-1]
define	ADJ_NUMBER	Memi[ADJ_NUMPTR($1)+($2)-1]
define	ADJ_NAME	Memc[ADJ_NAMARY($1,$2)]
define	ADJ_BACK	Memi[ADJ_BACKARY($1,$2)+($3)-1]
