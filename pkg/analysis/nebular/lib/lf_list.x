# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<error.h>
include "../flux.h"

#---------------------------------------------------------------------9 Jul 97--
.help lf_list.x Jul97 nebular/lib
.ih
NAME
.nf
   lfl_alloc - Allocate a flux list data structure
    lfl_free - Free a flux object list data structure
   lfl_debug - Dereference values from the flux object list for viewing
.fi
.endhelp
#-------------------------------------------------------------------------------
#  LFL_ALLOC --	Allocate the flux list data structure. 

pointer procedure lfl_alloc (size)

#  Arguments:
int	size			# I: initial no. elements to allocate
pointer	o			# O: flux list data object

begin
	call calloc (o, LEN_LFL, TY_STRUCT)

	if (IS_INDEFI(size) || size < 1)
	    LFL_A_SZ(o) = LFL_GROW
	else
	    LFL_A_SZ(o) = size

	# Allocate array of flux object pointers.
	call calloc (LFL_A_PTR(o), LFL_A_SZ(o), TY_POINTER)
	LFL_N(o) = 0

	return (o)
end


#-------------------------------------------------------------------------------
#  LFL_FREE -- Free flux object list data structure.

procedure lfl_free (o)

#  Calling arguments:
pointer	o		# I: flux data object list

#  Declarations:
int	i		# generic

begin
	if (o == NULL)
	    return

	do i = LFL_N(o), 1, -1 
	    call lf_free ( LFL_A(o, i) )

	call mfree (LFL_A_PTR(o), TY_POINTER)
	call mfree (o, TY_STRUCT)
end


#-------------------------------------------------------------------------------
#  LFL_DEBUG -- Dereference values from the flux object list data structure 
#		for viewing.

procedure lfl_debug (o)

#  Calling arguments:
pointer	o		# I: flux data object list

#  Declarations:
pointer	lf		# single flux object
int	i		# generic

begin
	if (o == NULL) {
	    call printf ("Flux object is NULL\n")
	    return
	}

	call printf ("No. flux objects in list: %d\n")
	    call pargi (LFL_N(o))

	do i = 1, LFL_N(o) {
	    lf = LFL_A(o, i)
	    call lf_debug ( LF_ATOM(lf), LF_ION(lf), LF_ZONE(lf), LF_ABUND(lf), 
		LF_N(lf), LF_WAVE(lf,1), LF_FLUX(lf,1), LF_WT(lf,1), 
		LF_NLINES(lf,1), LF_WIDTH(lf,1))
	}
end


