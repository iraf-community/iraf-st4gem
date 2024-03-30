include "libsynphot.h"

# ADJPARLIS -- Adjust the parameter list so it matches the component list

procedure adjparlis (mxparam, nmode, ncomp, nparam, 
		      modelist, modenum, paramlist)

int	mxparam				# i: Max number of params per keyword
int	nmode				# i: number of keywords in mode string
int	ncomp				# i: number of instrument components
int	nparam[ARB]			# u: number of parameters per keyword
char	modelist[SZ_KEYWRD,ARB]		# i: list of instrument mode keywords
int	modenum[ARB]			# u: component number of mode keyword
real	paramlist[mxparam,ARB]		# u: parameters of each keyword
#--
int	ic, imode, icomp, ipar
pointer	sp, tparnum, tparlist, badmodes, tparam

int	gstrcpy()

errchk	synphoterr

begin
	# Allocate memory for temporary arrays and strings

	call smark (sp)
	call salloc (tparnum, ncomp, TY_INT)
	call salloc (tparlist, mxparam*ncomp, TY_REAL)
	call salloc (badmodes, SZ_LINE, TY_CHAR)

	# Clear temporary arrays

	call amovki (0, Memi[tparnum], ncomp)
	call amovkr (0.0, Memr[tparlist], mxparam*ncomp)

	# Copy the parameters into the temporary arrays from the
	# locations corresponding to the mode list to the locations
	# corresponding to the component list

	# Create a list of mode keywords that do not have 
	# corresponding components

	ic = 0
	do imode = 1, nmode {
	    icomp = modenum[imode]

	    if (icomp > 0 && nparam[imode] > 0) {
		Memi[tparnum+icomp-1] = nparam[imode]
		tparam = tparlist + mxparam * (icomp - 1)

		do ipar = 1, nparam[imode] {
		    Memr[tparam] = paramlist[ipar,imode]
		    tparam = tparam + 1
		}

	    } else if (icomp == 0) {
		ic = ic + gstrcpy (modelist[1,imode], Memc[badmodes+ic], 
				   SZ_LINE-(ic+1))
		Memc[badmodes+ic] = ','
		ic = ic + 1
	    }
	}

	# Copy the temporary arrays back to the original arrays

	call amovi (Memi[tparnum], nparam, ncomp)
	call amovr (Memr[tparlist], paramlist, mxparam*ncomp)

	# Write a warning message if any mode keyword did not have 
	# a corresponding component

	if (ic > 0) {
	    Memc[badmodes+ic-1] = EOS
	    call modewarn (Memc[badmodes])
	}

	call sfree (sp)
	return

end
