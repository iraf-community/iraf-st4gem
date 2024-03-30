include	"../grflist.h"

define	MAXLIST		100		# Maximum number of keywords in list
define	MAXPARAM	3		# Max number of parameters / keyword

# STARTKEY -- Get starting keyword name from observation mode

procedure startkey (gp, mode, snode, keyname)

pointer	gp		# i: graph table descriptor
char	mode[ARB]	# i: instrument mode string
int	snode		# o: starting node number
char	keyname[ARB]	# o: starting keyword name
#--
int	nmode, complen, keylen
pointer	sp, modenum, nparam, paramlist, modelist

begin
	keylen = GRF_KEYLEN(gp)
	complen = GRF_CMPLEN(gp)

	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (modenum, MAXLIST, TY_INT)
	call salloc (nparam, MAXLIST, TY_INT)
	call salloc (paramlist, MAXLIST*MAXPARAM, TY_REAL)
	call salloc (modelist, MAXLIST*(keylen+1), TY_CHAR)

	# Break mode string into list of keywords and parameters

	call breakmode (mode, MAXLIST, MAXPARAM, keylen, nmode, 
			Memi[nparam], Memr[paramlist], Memc[modelist])

	# Use keywords to trace path thru instrument graph

	call startpath (gp, complen, keylen, nmode, Memc[modelist], 
			snode, keyname)

	call sfree (sp)
end
