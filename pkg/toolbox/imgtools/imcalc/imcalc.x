define	HARMLESS	0.1d0

# IMCALC -- Evalauate an expression using images as variables
#
# This procedure evaluates an expression passed to the task as a string
# parameter. The syntax of the expression is similar to a fortran 
# statement. The variables in the expression are of two forms: im1, im2, 
# etc., which stand for the corresponding image in the input list or 
# x, y, z which stand for first, second or third dimension of the image.
#
# B.Simon	18-May-90	Original
# B.Simon	02-May-91	Revised to handle multiple data types
# B.Simon	09-Aug-91	Fix for unusual data types

procedure t_imcalc ()

#--
pointer	input		# List of input images
pointer	output		# Name of output image
pointer	equals		# Expression to evaluate
pointer	pixtype		# Pixel type of output image
double	nullval		# Value to substitute for undefined expression
bool	verbose		# Print diagnostic messages?

int	itype, type, etype, pixval[6]
int	npixel, nline, iline, inull, oldfrac, frac
real	rnull
pointer	sp, code, buffer

string	pixstr	"|old|short|ushort|int|real|double|"
data	pixval	/ 0, TY_SHORT, TY_USHORT, TY_INT, TY_REAL, TY_DOUBLE /

bool	clgetb()
double	clgetd()
int	strdic(), init_imvar(), next_imvar()
pointer	vex_compile()

extern	get_imvar

begin
	# Allocate dynamic memory

	call smark (sp)
	call salloc (input, SZ_FNAME, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (equals, SZ_FNAME, TY_CHAR)
	call salloc (pixtype, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgstr ("input", Memc[input], SZ_FNAME)
	call clgstr ("output", Memc[output], SZ_FNAME)
	call clgstr ("equals", Memc[equals], SZ_FNAME)
	call clgstr ("pixtype", Memc[pixtype], SZ_FNAME)

	nullval = clgetd ("nullval")
	verbose = clgetb ("verbose")

	inull = nullval
	rnull = nullval

	# Initialize image variable arrays

	itype = strdic (Memc[pixtype], Memc[pixtype], SZ_FNAME, pixstr)
	type = pixval[itype]

	nline = init_imvar (Memc[input], Memc[output], Memc[equals], type)

	# Compile the expression, producing reverse polish code

	code = vex_compile (Memc[equals])

	# Print percent done message

	if (verbose) {
	    oldfrac = 0
	    call printf ("Percent done:")
	    call flush (STDOUT)
	}

	# Evaluate the expression a line at a time

	repeat {
	    iline = next_imvar (buffer, npixel)
	    if (iline == EOF)
		break

	    call vex_eval (code, get_imvar, HARMLESS, etype)

	    switch (type) {
	    case TY_SHORT,TY_INT,TY_LONG,TY_USHORT:
		call vex_copyi (code, inull, Memi[buffer], npixel)
	    case TY_REAL:
		call vex_copyr (code, rnull, Memr[buffer], npixel)
	    case TY_DOUBLE:
		call vex_copyd (code, nullval, Memd[buffer], npixel)
	    default:
		call vex_copyr (code, rnull, Memr[buffer], npixel)
	    }

	    # Calculate percent done and print

	    if (verbose) {
		frac = (100 * iline) / nline
		if ((frac - oldfrac) >= 10) {
		    oldfrac = frac
		    call printf (" %d")
		    call pargi (frac)
		    call flush (STDOUT)
		}
	    }
	}

	if (verbose) {
	    call printf ("\n")
	    call flush (STDOUT)
	}

	# Free memory

	call free_imvar
	call vex_free(code)
	call sfree (sp)
end
