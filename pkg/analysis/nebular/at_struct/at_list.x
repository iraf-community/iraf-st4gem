# Copyright restrictions apply - see st4gem$copyright.stsdas

include	"../at.h"
include	"../neberr.h"

# Memory management.
define	Nid		Memi[nid+($1)-1]
define	Oid		Memi[oid+($1)-1]

#---------------------------------------------------------------------6 Feb 97--
.help at_list Oct96 nebular/fivel
.ih
NAME
at_list -- Atomic data object list manipulation routines
.ih
DESCRIPTION
Along with 'at.h', this set of routines define the semantics of 
lists of atomic data objects.  Below is described the atomic 
data (AT) list structure along with the routines to manipulate 
a list.

An AT list is implemented as an array, each element of which
contains a pointer to a single AT.  The routines provide a
method of creating and destroying the list, adding and removing
elements, and of finding an element.

Also implemented is the concept of the "current data object".  
An AT object is current if it has just been added or 
queried for.
.ih
MEMORY STRUCTURE
This section describes the memory structure used to store the
AT list and any macros associated with access to the memory
structure.  All variables/macros begin with 'ATL_'.
.ls ATL_A_PTR(o)
This contains the pointer to the array of elements.  See also 
ATL_A below.
.le
.ls ATL_A(o,i)
This is the AT whose index in the array is 'i'.
.le
.ls ATL_A_SZ(o)
This contains the size of the array.  Since ATs can be
added/deleted dynamically, the size of the array is not necessarily
the same as the number of ATs stored in the array.  
.le
.ls ATL_N(o)
This holds the number of ATs currently stored in the list.  
.le
.ls ATL_SZ
Size of the memory structure.
.le
.ls ATL_GROW
The amount the array will grow (if necessary) when ATs 
are added.
.le
.ih
ROUTINE OVERVIEW
The following routines are available for manipulating AT lists.
Information contained in an AT list may also be accessed
directly from the memory structure.  However, for adding/deleting
elements from the list, it is strongly suggested that the routines
described below be utilized.

.nf
atl_alloc   -- Create an AT list
atl_free    -- Destroy an AT list
atl_add     -- Add an AT to the list
atl_del     -- Delete an AT from the list
atl_find    -- Find the index of an AT in the list, given its id
atl_get_cmp -- Get a pointer to an AT from the list, given its id
.fi
.ih
ROUTINE DETAILS
Below, each routine is described in detail.
.ls pointer atl_alloc (size)
Create an AT list
.ls ARGUMENTS
.ls size (input: integer)
Initial number of ATs the list will hold.  If unknown, specify
INDEFI.
.le
.le
.ls RETURNS
A pointer to the AT list structure.  This pointer should be
used in all macros and routines which use the AT list.
.le
.le
.ls atl_free (o, free)
Destroy an AT list.
.ls ARGUMENTS
.ls o (input/output: pointer)
The pointer to the AT list structure to be freed.  On return,
the value will be NULL.
.le
.ls free (input: bool)
If 'true', each AT contained in the list will be destroyed 
using the 'at_free' routine.  If false, the ATs will not be destroyed.  
This is useful if there are other structures which point to the 
same ATs. 
.le
.le
.le
.ls atl_add (o, c, newid)
Add an AT to the list.
.ls ARGUMENTS
.ls o (input: pointer)
The AT list pointer.
.le
.ls c (input: pointer)
The pointer to the AT to add to the list.
.le
.ls newid (input: bool)
If 'true', the AT will be assinged a new identifier.  If false,
the AT will retain its identifier.
.le
.le
.le
.ls atl_del (o, id, free)
Remove an AT from the list based on the AT id.
.ls ARGUMENTS
.ls o (input: pointer)
The AT list pointer.
.le
.ls id (input: int)
The id of the AT to remove from the list.  If INDEFI, the
current AT will be removed.  If there is no current AT,
an error will be generated.
.le
.ls free (input: bool)
If 'true', the AT will be destroyed using the 'at_free' routine.
Else, only the list reference to the AT will be destroyed,
leaving the AT itself intact.  Useful if other memory
structures refer to the same AT.
.le
.le
.le
.ls pointer atl_copy (o, copy)
Create a new AT list that is a copy of the given one.
.ls ARGUMENTS
.ls o (input: pointer)
The AT list to be copied.
.le
.ls int atl_find (o, id)
Get the index into the array of an AT with the specified id.
.ls ARGUMENTS
.ls o (input: pointer)
The AT list to search.
.le
.ls id (input: int)
The id of the AT to find.  If INDEFI, the current AT will be 
returned.  If there is no current AT, an error will be generated.
.le
.le
.ls RETURNS
An index into the array pointing to the specified AT.
.le
.le
.ls pointer atl_get_at (o, atom, ion)
Get the AT from the list with the specified id.
.ls ARGUMENTS
.ls o (input: pointer)
The AT list to search.
.le
.ls atom (input: int)
The atom of the AT to retrieve.  
.le
.ls ion (input: int)
The ion of the AT to retrieve.  
.le
.le
.ls RETURNS
The pointer to the desired AT.
.le
.le
.ls int aindxi (v, a, n)
Find a value in an array and return its index position in the array.
If the value is not found, an error is generated.
.ls ARGUMENTS
.ls v (input: int)
The value to find in the array.
.le
.ls a (input: int[n])
The array to look for the value in.
.le
.ls n (input: int)
The number of elements in the array
.le
.le
.ls RETURNS (int)
The index into the array where the value is located.  
.le
.le
.endhelp
#---------------------------------------------------------------------------
pointer procedure atl_alloc (size)

int	size			# I: Initial number of elements to allocate.

# Declarations
pointer	o			# AT list object.

errchk	malloc

begin
	call malloc (o, ATL_SZ, TY_STRUCT)
	if (IS_INDEFI(size))
	    ATL_A_SZ(o) = ATL_GROW
	else
	    ATL_A_SZ(o) = size

	call malloc (ATL_A_PTR(o), ATL_A_SZ(o), TY_POINTER)
	ATL_N(o) = 0

	return (o)
end
#---------------------------------------------------------------------------
# End of atl_alloc
#---------------------------------------------------------------------------
procedure atl_free (o, free)

pointer	o			# IO: AT List Object, NULL on exit.
bool	free			# I:  'true' to free the individual elements.

# Declarations
int	i			# Generic.

errchk	at_free, mfree

begin
	# Make sure there is a list to free.
	if (o == NULL)
	    return
	
	if (free)
	    do i = 1, ATL_N(o)
		call at_free (ATL_A(o, i))

	call mfree (ATL_A_PTR(o), TY_POINTER)
	call mfree (o, TY_STRUCT)
end
#---------------------------------------------------------------------------
# End of atl_free
#---------------------------------------------------------------------------
procedure atl_add (o, c)

pointer	o			# I:  AT list object.
pointer	c			# I:  AT to add.

errchk	realloc

begin
	ATL_N(o) = ATL_N(o) + 1
	if (ATL_N(o) > ATL_A_SZ(o)) {
	    ATL_A_SZ(o) = ATL_A_SZ(o) + ATL_GROW
	    call realloc (ATL_A_PTR(o), ATL_A_SZ(o), TY_POINTER)
	}

	ATL_A(o, ATL_N(o)) = c
end
#---------------------------------------------------------------------------
# End of atl_add
#---------------------------------------------------------------------------
procedure atl_del (o, atom, ion, free)

pointer	o			# I:  AT list object.
int	atom			# I:  atom of AT.
int	ion			# I:  ion of AT.
bool	free			# I:  'true' to free AT.

# Declarations
int	atl_find()		# Get index of specified AT.
int	i			# Index of compoent to remove

errchk	at_free, atl_find

begin
	i = atl_find (o, atom, ion)

	if (free)
	    call at_free (ATL_A(o,i))
	
	call amovi (ATL_A(o, i+1), ATL_A(o,i), ATL_N(o) - i)
	ATL_N(o) = ATL_N(o) - 1
end
#---------------------------------------------------------------------------
# End of atl_del
#---------------------------------------------------------------------------
int procedure atl_find (o, atom, ion)

pointer	o			# I: AT list object.
int	atom			# I: atom of AT to find.
int	ion			# I: ion of AT to find.

# Declarations.
int	i			# Generic.
char	sx[SZ_LINE]		# Generic.

begin
	
	# Look for the AT.
	do i = 1, ATL_N(o) {
	    if ( atom == AT_ATOM(ATL_A(o, i)) ) {
	    	if ( ion == AT_ION(ATL_A(o, i))) 
		    return (i)
	    }
	}

	call sprintf (sx, SZ_LINE, 
	    "atl_find: no ATs for atom %d, ion %d")
	    call pargi (atom)
	    call pargi (ion)
	call error (ATL_NOSUCHID, sx)
end
#---------------------------------------------------------------------------
# End of atl_find
#---------------------------------------------------------------------------
pointer	procedure atl_get_at (o, atom, ion)

pointer	o			# I: AT list object.
int	atom			# I: atom of AT to find.
int	ion			# I: ion of AT to find.

# Declarations
int	atl_find()		# Find index of AT in array.
int	i			# Generic.

errchk	atl_find

begin
	iferr (i = atl_find (o, atom, ion) )
	    return (NULL)
	else
	    return (ATL_A(o, i))
end
#---------------------------------------------------------------------------
# End of atl_get_at
#---------------------------------------------------------------------------
int procedure aindex (v, a, n)

int	v			# I: The value to find in the array.
int	a[n]			# I: The array of values.
int	n			# I: Number of values.

# Declarations
int	i			# Generic.

begin
	do i = 1, n
	    if (a[i] == v)
		return (i)

	call error (1, "aindex: value not found in array")
end
#---------------------------------------------------------------------------
# End of aindex
#---------------------------------------------------------------------------
