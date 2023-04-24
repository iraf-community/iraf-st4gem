include "addmasks.h"

define	MAX_STACK	32

# ADDEVAL -- Evaluate an expression for combining mask files
#
# B.Simon	16-Jul-93	Original

procedure addeval (pcode, iline, flag, order, nflag, 
		   inbuf, outbuf, lenbuf, nfile)

pointer	pcode		# i: Pseudocode to be evaluated
int	iline		# i: Mask line, for error reporting
int	flag[ARB]	# i: Array of flag values used in mask
int	order[ARB]	# i: Precedence of flag values
int	nflag		# i: Number of flags
pointer	inbuf[ARB]	# i: Array of input buffers
pointer	outbuf		# o: Output buffer
int	lenbuf		# i: Length of line
int	nfile		# i: Number of input files
#--
int	hi, top, nop, index
pointer	sp, op, stack, temp

string	badindex  "Variable in expression exceeds number of input images"
string	badexpr   "Expression too complex"
string	badcode   "Illegal code in expression evaluator"

begin
	# Allocate stack for intermediate results

	call smark (sp)
	call salloc (stack, MAX_STACK, TY_INT)
	
	hi = -1
	top = -1
	op = pcode

	while (Memi[op] != Y_DONE) {

	    # Allocate array on stack for result of calculation

	    top = top + 1
	    if (top == MAX_STACK)
		call error (1, badexpr)

	    if (top > hi) {
		hi = hi + 1
		call salloc (Memi[stack+hi], lenbuf, TY_INT)
	    }

	    # Perform one step in calculation

	    switch (Memi[op]) {
	    case Y_THEN:
		nop = 2
		call addthen (Memi[stack+top-2], Memi[stack+top-1], 
			      Memi[stack+top], lenbuf)
	    case Y_ELSE:
		nop = 3
		call addelse (Memi[stack+top-3], Memi[stack+top-2], 
			      Memi[stack+top-1], Memi[stack+top], lenbuf)
	    case Y_PUSH:
		op = op + 1
		index = Memi[op]
		if (index < 1 || index > nfile) {
		    call error (1, badindex)
		}

		nop = 0
		call addpush (flag, nflag, index, iline, inbuf[index], 
			     Memi[stack+top], lenbuf)
	    case Y_CONST:
		nop = 0
		op = op + 1
		call addconst (flag, nflag, op, Memi[stack+top], lenbuf)
	    case Y_OR:
		nop = 2
		call addor (Memi[stack+top-2], Memi[stack+top-1], 
			    Memi[stack+top], lenbuf)
	    case Y_AND:
		nop = 2
		call addand (Memi[stack+top-2], Memi[stack+top-1], 
			     Memi[stack+top], lenbuf)
	    case Y_NOT:
		nop = 1
		call addnot (nflag, Memi[stack+top-1], Memi[stack+top], lenbuf)
	    case Y_EQ:
		nop = 2
		call addeq (nflag, Memi[stack+top-2], Memi[stack+top-1], 
			    Memi[stack+top], lenbuf)
	    case Y_NE:
		nop = 2
		call addne (nflag, Memi[stack+top-2], Memi[stack+top-1], 
			    Memi[stack+top], lenbuf)
	    case Y_LT:
		nop = 2
		call addlt (order, nflag, Memi[stack+top-2], 
			    Memi[stack+top-1], Memi[stack+top], lenbuf)
	    case Y_GT:
		nop = 2
		call addgt (order, nflag, Memi[stack+top-2], 
			    Memi[stack+top-1], Memi[stack+top], lenbuf)
	    case Y_LE:
		nop = 2
		call addle (order, nflag, Memi[stack+top-2], 
			    Memi[stack+top-1], Memi[stack+top], lenbuf)
	    case Y_GE:
		nop = 2
		call addge (order, nflag, Memi[stack+top-2], 
			    Memi[stack+top-1], Memi[stack+top], lenbuf)
	    default:
		call error (1, badcode)
	    }

	    # Swap output array with first input array
	    # so that input arrays are beyond new top of stack

	    if (nop > 0) {
		temp = Memi[stack+top-nop]
		Memi[stack+top-nop] = Memi[stack+top]
		Memi[stack+top] = temp
		top = top - nop
	    }

	    op = op + 1
	}

	# Copy result to output buffer

	call addpop (flag, nflag, Memi[stack], outbuf, lenbuf)

	call sfree (sp)
end

# ADDTHEN -- If in1 != 0 then out = in2 else out = 0

procedure addthen (in1, in2, out, lenbuf)

pointer	in1		# i: First input array
pointer	in2		# i: Second input array
pointer	out		# o: Output array
int	lenbuf		# i: Length of arrays
#--
int	ibuf

begin
	do ibuf = 0, lenbuf-1 {
	    if (Memi[in1+ibuf] != 0) {
		Memi[out+ibuf] = Memi[in2+ibuf]
	    } else {
		Memi[out+ibuf] = 0
	    }
	}
end

# ADDELSE -- If in1 != 0 then out = in2 else out = in3

procedure addelse (in1, in2, in3, out, lenbuf)

pointer	in1		# i: First input array
pointer	in2		# i: Second input array
pointer	in3		# i: Third input array
pointer	out		# o: Output array
int	lenbuf		# i: Length of arrays
#--
int	ibuf

begin
	do ibuf = 0, lenbuf-1 {
	    if (Memi[in1+ibuf] != 0) {
		Memi[out+ibuf] = Memi[in2+ibuf]
	    } else {
		Memi[out+ibuf] = Memi[in3+ibuf]
	    }
	}
end

# ADDPUSH -- Convert mask line to bit flags, copy to out

procedure addpush (flag, nflag, index, iline, in, out, lenbuf)

int	flag[ARB]	# i: Array of flag values
int	nflag		# i: Number of flags
int	index		# i: Input image number, for reporting
int	iline		# i: Current line in mask, for reporting
pointer	in		# i: Buffer containing mask line
pointer	out		# o: Output buffer
int	lenbuf		# i: Length of buffer
#--
int	ibuf
int	mapflag()

begin
	do ibuf = 0, lenbuf-1 {
	    if (Memi[in+ibuf] == 0) {
		Memi[out+ibuf] = 0

	    } else {
		if (mapflag (flag, nflag, 
			     Memi[in+ibuf], Memi[out+ibuf]) == ERR) {

		    Memi[out+ibuf] = 0
		    call sendwarn (index, iline, ibuf+1, Memi[in+ibuf])
		}
	    }
	}
end

# ADDCONST -- Read constant from code, convert to bit flags, copy to out

procedure addconst (flag, nflag, op, out, lenbuf)

int	flag[ARB]	# i: Array of flag values
int	nflag		# i: Number of flags
pointer	op		# u: Pointer to start of constant
pointer	out		# o: Output buffer
int	lenbuf		# i: Length of buffer
#--
int	ic, junk, val, const
pointer	sp, str

string	badconst "Constant has illegal value (%d). Treated as if 0.\n"

int	ctoi(), mapflag()

begin
	call smark (sp)
	call salloc (str, SZ_LINE, TY_CHAR)

	# Copy constant from code into string
	# Code pointer is advanced to end of string

	for (ic = 0; ic < SZ_LINE && Memi[op] != EOS; ic = ic + 1) {
	    Memc[str+ic] = Memi[op]
	    op = op + 1
	}
	Memc[str+ic] = EOS

	# Convert string to array of bit flags

	ic = 1
	junk = ctoi (Memc[str], ic, val)
	if (mapflag (flag, nflag, val, const) == ERR) {
	    const = 0
	    call eprintf (badconst)
	    call pargi (val)
	}

	# Copy bit flags to output array

	call amovki (const, Memi[out], lenbuf)
	call sfree (sp)

end

# ADDOR -- Bitwise or: out = or(in1,in2)

procedure addor (in1, in2, out, lenbuf)

pointer	in1		# i: First input array
pointer	in2		# i: Second input array
pointer	out		# o: Output array
int	lenbuf		# i: Length of arrays
#--

begin
	call abori (Memi[in1], Memi[in2], Memi[out], lenbuf)
end

# ADDAND -- Bitwise and: out = and(in1,in2)

procedure addand (in1, in2, out, lenbuf)

pointer	in1		# i: First input array
pointer	in2		# i: Second input array
pointer	out		# o: Output array
int	lenbuf		# i: Length of arrays
#--

begin
	call aandi (Memi[in1], Memi[in2], Memi[out], lenbuf)
end

# ADDNOT -- Bitwise not: out = not(in1)

procedure addnot (nflag, in1, out, lenbuf)

int	nflag		# i: Number of flags 
pointer	in1		# i: First input array
pointer	out		# o: Output array
int	lenbuf		# i: Length of arrays
#--
int	tval
int	maptrue()

begin
	tval = maptrue (nflag)

	call anoti (Memi[in1], Memi[out], lenbuf)
	call aandki (Memi[out], tval, Memi[out], lenbuf)
end

# ADDEQ -- if in1 == in2 then out = true else out = false

procedure addeq (nflag, in1, in2, out, lenbuf)

int	nflag		# i: Number of flags 
pointer	in1		# i: First input array
pointer	in2		# i: Second input array
pointer	out		# o: Output array
int	lenbuf		# i: Length of arrays
#--
int	ibuf, tval
int	maptrue()

begin
	tval = maptrue (nflag)

	do ibuf = 0, lenbuf-1 {
	    if (Memi[in1+ibuf] == Memi[in2+ibuf]) {
		Memi[out+ibuf] = tval
	    } else {
		Memi[out+ibuf] = 0
	    }
	}
end

# ADDNE -- if in1 != in2 then out = true else out = false

procedure addne (nflag, in1, in2, out, lenbuf)

int	nflag		# i: Number of flags 
pointer	in1		# i: First input array
pointer	in2		# i: Second input array
pointer	out		# o: Output array
int	lenbuf		# i: Length of arrays
#--
int	ibuf, tval
int	maptrue()

begin
	tval = maptrue (nflag)

	do ibuf = 0, lenbuf-1 {
	    if (Memi[in1+ibuf] != Memi[in2+ibuf]) {
		Memi[out+ibuf] = tval
	    } else {
		Memi[out+ibuf] = 0
	    }
	}
end

# ADDLT -- if in1 < in2 then out = true else out = false

procedure addlt (order, nflag, in1, in2, out, lenbuf)

int	order[ARB]	# i: Precedence order of flags
int	nflag		# i: Number of flags 
pointer	in1		# i: First input array
pointer	in2		# i: Second input array
pointer	out		# o: Output array
int	lenbuf		# i: Length of arrays
#--
int	ibuf, tval
int	maptrue(), pixorder()

begin
	tval = maptrue (nflag)

	do ibuf = 0, lenbuf-1 {
	    if (pixorder (order, nflag, Memi[in1+ibuf], Memi[in2+ibuf]) < 0) {
		Memi[out+ibuf] = tval
	    } else {
		Memi[out+ibuf] = 0
	    }
	}
end

# ADDGT -- if in1 > in2 then out = true else out = false

procedure addgt (order, nflag, in1, in2, out, lenbuf)

int	order[ARB]	# i: Precedence order of flags
int	nflag		# i: Number of flags 
pointer	in1		# i: First input array
pointer	in2		# i: Second input array
pointer	out		# o: Output array
int	lenbuf		# i: Length of arrays
#--
int	ibuf, tval
int	maptrue(), pixorder()

begin
	tval = maptrue (nflag)

	do ibuf = 0, lenbuf-1 {
	    if (pixorder (order, nflag, Memi[in1+ibuf], Memi[in2+ibuf]) > 0) {
		Memi[out+ibuf] = tval
	    } else {
		Memi[out+ibuf] = 0
	    }
	}
end

# ADDLE -- if in1 <= in2 then out = true else out = false

procedure addle (order, nflag, in1, in2, out, lenbuf)

int	order[ARB]	# i: Precedence order of flags
int	nflag		# i: Number of flags 
pointer	in1		# i: First input array
pointer	in2		# i: Second input array
pointer	out		# o: Output array
int	lenbuf		# i: Length of arrays
#--
int	ibuf, tval
int	maptrue(), pixorder()

begin
	tval = maptrue (nflag)

	do ibuf = 0, lenbuf-1 {
	    if (pixorder (order, nflag, Memi[in1+ibuf], Memi[in2+ibuf]) <= 0) {
		Memi[out+ibuf] = tval
	    } else {
		Memi[out+ibuf] = 0
	    }
	}
end

# ADDGE -- if in1 >= in2 then out = true else out = false

procedure addge (order, nflag, in1, in2, out, lenbuf)

int	order[ARB]	# i: Precedence order of flags
int	nflag		# i: Number of flags 
pointer	in1		# i: First input array
pointer	in2		# i: Second input array
pointer	out		# o: Output array
int	lenbuf		# i: Length of arrays
#--
int	ibuf, tval
int	maptrue(), pixorder()

begin
	tval = maptrue (nflag)

	do ibuf = 0, lenbuf-1 {
	    if (pixorder (order, nflag, Memi[in1+ibuf], Memi[in2+ibuf]) >= 0) {
		Memi[out+ibuf] = tval
	    } else {
		Memi[out+ibuf] = 0
	    }
	}
end

# ADDPOP -- Convert bit flags to mask line, copy to out

procedure addpop (flag, nflag, in, out, lenbuf)

int	flag[ARB]	# i: Array of flag values
int	nflag		# i: Number of flags
pointer	in		# i: Buffer containing mask line
pointer	out		# o: Output buffer
int	lenbuf		# i: Length of buffer
#--
int	ibuf

begin
	do ibuf = 0, lenbuf-1 {
	    if (Memi[in+ibuf] == 0) {
		Memi[out+ibuf] = 0
	    } else {
		call unmapflag (flag, nflag, Memi[in+ibuf], Memi[out+ibuf])
	    }
	}
end
