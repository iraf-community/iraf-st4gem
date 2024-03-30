include "function.h"

# SHAPEFUNCS -- Return pointers to shape functions

procedure shapefuncs (func, newfunc, delfunc, valfunc, sizfunc)

pointer	func		# i: the function descriptor
pointer	newfunc		# o: pointer to initialization function
pointer	delfunc		# o: pointer to delete function
pointer	valfunc		# o: pointer to evaluation function
pointer	sizfunc		# o: pointer to size function
#--
int	index

string	funclist "gauss,moffat,exp,devauc,prof,img"
string  badfunc  "Unrecognized shape function"

int	word_match(), locpr()

extern	newgauss, delgauss, valgauss, sizgauss
extern	newmoffat, delmoffat, valmoffat, sizmoffat
extern	newexp, delexp, valexp, sizexp
extern	newdevauc, deldevauc, valdevauc, sizdevauc
extern	newlist, dellist, vallist, sizlist
extern	newimage, delimage, valimage, sizimage

begin
	newfunc = NULL
	delfunc = NULL
	valfunc = NULL
	sizfunc = NULL

	if (func == NULL)
	    return

	# Select proper functions

	call strfix (FUN_STR(func,1))
	index = word_match (FUN_STR(func,1), funclist)

	switch (index) {
	case 1: # gaussian shape
	    newfunc = locpr (newgauss)
	    delfunc = locpr (delgauss)
	    valfunc = locpr (valgauss)
	    sizfunc = locpr (sizgauss)

	case 2: # moffat shape
	    newfunc = locpr (newmoffat)
	    delfunc = locpr (delmoffat)
	    valfunc = locpr (valmoffat)
	    sizfunc = locpr (sizmoffat)

	case 3: # exponential shape
	    newfunc = locpr (newexp)
	    delfunc = locpr (delexp)
	    valfunc = locpr (valexp)
	    sizfunc = locpr (sizexp)

	case 4: # Devaucalour shape
	    newfunc = locpr (newdevauc)
	    delfunc = locpr (deldevauc)
	    valfunc = locpr (valdevauc)
	    sizfunc = locpr (sizdevauc)

	case 5: # Read shape from file
	    newfunc = locpr (newlist)
	    delfunc = locpr (dellist)
	    valfunc = locpr (vallist)
	    sizfunc = locpr (sizlist)

	case 6: # Read shape from image
	    newfunc = locpr (newimage)
	    delfunc = locpr (delimage)
	    valfunc = locpr (valimage)
	    sizfunc = locpr (sizimage)

	default: # OOPS !!
	    call printerr_str (badfunc, FUN_STR(func,1))
	}


end
