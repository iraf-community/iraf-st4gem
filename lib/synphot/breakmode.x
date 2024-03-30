include "libsynphot.h"

# BREAKMODE -- Break instrument mode string into separate keywords

procedure breakmode (mode, mxlist, mxparam, keylen, nmode, 
		     nparam, paramlist, modelist)

char	mode[ARB]			# i: instrument mode string
int	mxlist				# i: Maximum number of keywords in list
int	mxparam				# i: Max number of params per keyword
int	keylen				# i: Max length keyword string
int	nmode				# o: number of keywords in mode string
int	nparam[ARB]			# o: number of parameters per keyword
real	paramlist[mxparam,ARB]		# o: parameters of each keyword
char	modelist[keylen,ARB]		# o: list of instrument mode keywords
#--
int	ic

string	overlist  "Too many keywords in list. Last keyword read"

int	word_fetch()
errchk	synphoterr, breakparam

begin
	# Break instrument mode string into keywords

	ic = 1
	for (nmode = 1; nmode <= mxlist; nmode = nmode + 1) {
	    if (word_fetch (mode, ic, modelist[1,nmode], keylen) == 0) {
		nmode = nmode - 1
		break
	    }

	    # Extract any parameters from keyword string

	    call strlwr (modelist[1,nmode])

	    call breakparam (modelist[1,nmode], mxparam, nparam[nmode], 
			     paramlist[1,nmode])
	}

	if (nmode > mxlist)
	    call synphoterr (overlist, modelist[1,mxlist])

end
