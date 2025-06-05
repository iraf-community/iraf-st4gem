include	<iraf77.h>

# UCLGV? -- Get an array of parameters from the CL.
# The way the cl currently stores arrays the parameters have to be retrieved
# individually by names of the form param[001], param[002], ....
# Note that due to cl restrictions nelems + felem cannot be greater than
# the actual size of the array.
#
# Phil Hodge,  3-Feb-1989  modified to delete the call to clgstr
# Phil Hodge,  8-Mar-1989  modified:  don't use ' ' in stridx

procedure uclgvs (f77par, felem, nelems, nvals, buf, istat)

% 	character*(*) f77par
short	buf[ARB]
int	felem			# first element in array
int	nelems			# number of elements
int	nvals			# mumber of elements actually retrieved
int	istat

char 	parnam[SZ_FNAME]
char 	outstr[SZ_FNAME]
char	blank
int	opbrkt, index, stridx(), strlen()
short	clgets()
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
	nvals = 0
	index = felem

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
	      call pargi (index)
	   call strcat (outstr, parnam[opbrkt], SZ_FNAME)
	   call strcat ("]", parnam[opbrkt], SZ_FNAME)

	   # get parameter value
	   value = clgets (parnam)
	   if (IS_INDEFS (value))  {
	      istat = ER_CLUNDF
	      return
	   }
	   buf[nvals + 1] = value
	   index = index + 1
	   nvals = nvals + 1
	   parnam[opbrkt] = EOS
	} until (nvals == nelems)

	return
end
