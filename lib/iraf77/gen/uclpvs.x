include	<iraf77.h>

# UCLPV? -- Put an array of parameters into the CL.
# The way the cl currently stores arrays the parameters have to be retrieved
# individually by names of the form param[001], param[002], ....
# Note that due to cl restrictions nelems + felem cannot be greater than
# the actual size of the array.
#
# Phil Hodge,  8-Mar-1989  modified:  don't use ' ' in stridx

procedure uclpvs (f77par, buf, felem, nelems, istat)

% 	character*(*) f77par
short	buf[ARB]
int	felem			# first element in array
int	nelems			# number of elements
int	istat

char 	parnam[SZ_FNAME]
char 	outstr[SZ_FNAME]
char	blank
int	opbrkt, parindex, bufindex, stridx(), strlen(), nvals
# PIXEL	clput$t()
short	value

begin
	if (felem <= 0)  {
	   istat = ER_CLBADFELEM
	   return
	}
	if (nelems <= 0)  {
	   istat = ER_CLBADNELEMS
	   return
	}
	istat = ER_OK
	parindex = felem
	bufindex = 1
	nvals = 0

	# Convert character string to SPP string
	call f77upk (f77par, parnam, SZ_FNAME)

	# Build parameter's root (parnam[)
	blank = ' '
	opbrkt = stridx (blank, parnam)
	if (opbrkt == 0)
	   opbrkt = strlen (parnam) + 1
	parnam[opbrkt] = '['
	opbrkt = opbrkt + 1
	parnam[opbrkt] = EOS

	repeat {
	   call sprintf (outstr, SZ_FNAME, "%03d")
	      call pargi (parindex)
	   call strcat (outstr, parnam[opbrkt], SZ_FNAME)
	   call strcat ("]", parnam[opbrkt], SZ_FNAME)

	   value = buf[bufindex]
	   if (IS_INDEFS (value))  {
	      istat = ER_CLUNDF
	      return
	   }
	   call clputs (parnam, value)
	   parindex = parindex + 1
	   bufindex = bufindex + 1
	   nvals = nvals + 1
	   parnam[opbrkt] = EOS
	} until (nvals == nelems)

	return
end
