include	<ctype.h>
include	"../adjlist.h"

# GETNODE -- Get the node associated with a component name

# This procedure returns the first node which matches an input component 
# name. The component name string may optionally contain a input node number,
# which allows the user to specify nodes with non-unique names. If the
# the component name string is blank, the procedure returns zero, which
# indicates that the program should take its default action.
#
# B.Simon	21-Jul-88	First Code
# B.Simon	26-Sep-91	Fixed upper limit on loop over nodes

int procedure getnode (gp, compname)

pointer	gp		# i: Pointer to adjacency list structure
char	compname[ARB]	# i: Component name
#--
bool	found
int	ic, list, number
pointer	sp, name

string	badnum	  "Input node number is invalid"
string	notfound  "Component name not found in graph"

bool	streq()
int	ctowrd(), ctoi()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (name, SZ_LINE, TY_CHAR)

	# The first word in the compname string is the component name.
	# If the name is not found, assume the user wants to take
	# default action.

	ic = 1
	if (ctowrd (compname, ic, Memc[name], SZ_LINE) == 0)
	    list = 0

	else {

	    # Convert the component name to lower case
	    # for a case insensitive match

	    call strlwr (Memc[name])

	    # Optionally, the compname string may be further qualified
	    # by a input node number. This allows the user to distinguish
	    # between nodes with the same component name.

	    if (ctoi (compname, ic, number) == 0)
		number = INDEFI

	    if (! (IS_WHITE(compname[ic]) || compname[ic] == EOS))
		call printerr_str (badnum, compname)

	    # Find the first node which matches the component name and
	    # optionally, the input node number

	    list = 0
	    found = false
	    while (list < ADJ_SIZE(gp) && ! found) {
		list = list + 1
		found = streq (ADJ_NAME(gp,list), Memc[name]) &&
			(ADJ_NUMBER(gp,list) == number || IS_INDEFI (number))
	    }

	    if (! found)
		call printerr_str (notfound, compname)
	}

	call sfree (sp)
	return (list)
end
