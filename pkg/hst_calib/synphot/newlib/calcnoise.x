include	"calcnoise.h"

#* HISTORY *
#* B.Simon	30-Mar-95	original

# CALCNOISE -- Calculate the mean noise level from a noise expression

procedure calcnoise (pcode, nflux, flux, mean)

int	pcode[ARB]	# i: pseudocode array
int	nflux		# i: length of flux array
real	flux[ARB]	# i: flux
real	mean[ARB]	# o: mean noise
#--
int	hi, top, icode, narg, bot
pointer temp, stack[SZ_STACK]
real	const

string	name  "calcnoise"
string	badexpr  "Expression too complex"

real    badvalmsg() # Print a message if an illegal value is encountered
extern	badvalmsg
errchk	asqrr, allnr, advzr

begin
	# Initialize stack and code indexes

	hi = 0
	top = 0
	icode = 1

	# Loop over code instructions until stop instruction

	while (pcode[icode] != N_DONE) {
	    # Allocate memory for intermediate result

	    top = top + 1
	    if (top > SZ_STACK)
		call synphoterr (badexpr, name)

	    if (top > hi) {
		hi = hi + 1
		call salloc (stack[top], nflux, TY_REAL)
	    }

	    # Evaluate a single instruction

	    switch (pcode[icode]) {
	    case N_SWAP:
		# Swap the two top arguments on a stack
		narg = 0
		bot = top - 2
		top = top - 1
		temp = stack[bot]
		stack[bot] = stack[top]
		stack[top] = temp
		
	    case N_NUM:
		# Push a constant embedded in the code on the stack
		narg = 0
		call get_const (pcode, icode, const)
		call amovkr (const, Memr[stack[top]], nflux)

	    case N_FLUX:
		# Push the flux on the stack
		narg = 0
		call amovr (flux, Memr[stack[top]], nflux)

	    case N_SQRT:
		# Square root
		narg = 1
		call asqrr (Memr[stack[top-1]], Memr[stack[top]], 
			    nflux, badvalmsg)
	    case N_LOG:
		# Natural logarithm
		narg = 1
		call allnr (Memr[stack[top-1]], Memr[stack[top]], 
			    nflux, badvalmsg)

	    case N_MAG:
		# Magnitude (sqrt (a**2 + b**2))
		narg = 2
		call amagr (Memr[stack[top-2]], Memr[stack[top-1]],
			    Memr[stack[top]], nflux)
	    case N_ADD:
		# Addition
		narg = 2
		call aaddr (Memr[stack[top-2]], Memr[stack[top-1]],
			    Memr[stack[top]], nflux)

	    case N_SUB:
		# Subtraction
		narg = 2
		call asubr (Memr[stack[top-2]], Memr[stack[top-1]],
			    Memr[stack[top]], nflux)
	    case N_MUL:
		# Multiplication
		narg = 2
		call amulr (Memr[stack[top-2]], Memr[stack[top-1]],
			    Memr[stack[top]], nflux)

	    case N_DIV:
		# Division
		narg = 2
		call advzr (Memr[stack[top-2]], Memr[stack[top-1]],
			    Memr[stack[top]], nflux, badvalmsg)
	    case N_EXP:
		# Exponentiation
		narg = 2
		call aexpr (Memr[stack[top-2]], Memr[stack[top-1]],
			    Memr[stack[top]], nflux)
	    case N_NEG:
		# Negation
		narg = 1
		call anegr (Memr[stack[top-1]], Memr[stack[top]], nflux)
	    }

	    # Remove evaluated arguments from stack

	    if (narg > 0) {
		bot = top - narg
		temp = stack[bot]
		stack[bot] = stack[top]
		stack[top] = temp

		top = bot
	    }

	    icode = icode + 1
	}

	if (top <= 0) {
	    call amovkr (0.0, mean, nflux)
	} else {
	    call amovr (Memr[stack[1]], mean, nflux)
	}

end

# GET_CONST -- Extract a string from the code and convert it to a real

procedure get_const (pcode, icode, const)

int	pcode[ARB]	# i: pseudocode array
int	icode		# u: pointer to current code
real	const		# o: constant read feom code
#--
int	ic, nc
pointer	sp, token

int	ctor()

begin
	# Allocate array to hold string version of constant

	call smark (sp)
	call salloc (token, SZ_TOKEN, TY_CHAR)

	# Copy string constant from code array

	ic = 0
	icode = icode + 1

	while (pcode[icode] != EOS) {
	    Memc[token+ic] = pcode[icode]
	    icode = icode + 1
	    ic = ic + 1
	}

	Memc[token+ic] = EOS

	# Convert string constant to real number

	ic = 1
	nc = ctor (Memc[token], ic, const)

	call sfree (sp)
end

# BADVALMSG -- Print a message if an illegal value is encountered

real procedure badvalmsg (value)

real	value		# i: Value which is outside of domain
#--
pointer	sp, errval

string	badvalue  "Illegal value when evaluating noise expression"

errchk	synphoterr

begin
	call smark (sp)
	call salloc (errval, SZ_TOKEN, TY_CHAR)

	call sprintf (Memc[errval], SZ_TOKEN, "%g")
	call pargr (value)

	call synphoterr (badvalue, Memc[errval])
	call sfree (sp)

	return (0.0)
end
