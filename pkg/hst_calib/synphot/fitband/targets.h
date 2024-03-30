# Common block used to pass target ids to bandfunk
include	"../plspec/plspec.h"
int	ntarget
char	targetid[SZ_FNAME,MAXPHOT]
char	specfile[SZ_FNAME,MAXPHOT]
common	/targetcom/ ntarget, targetid, specfile
