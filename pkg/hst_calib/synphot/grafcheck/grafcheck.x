include	"../grflist.h"

# GRAFCHECK -- Check the validity of a graph table 
#
# B.Simon	22-Jul-88	First Code

procedure grafcheck ()

pointer	grftable	# i: Graph table name
#--
pointer	sp, gp

pointer	grf_open()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (grftable, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgstr ("grftable", Memc[grftable], SZ_FNAME)

	# Read the instrument graph into a graph list structure

	gp = grf_open (Memc[grftable])

	# Check the validity of the graph list structure

	call nullchk (gp)
	call backchk (gp)
	call ambigchk (gp)
	call cnectchk (gp)

	# Free dynamic memory

	call grf_close (gp)
	call sfree (sp)

end
