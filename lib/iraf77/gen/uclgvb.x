include	<iraf77.h>

# UCLGV? -- Get an array of parameters from the CL.
# The way the cl currently stores arrays the parameters have to be retrieved
# individually by names of the form param[001], param[002], ....
# Note that due to cl restrictions nelems + felem cannot be greater than
# the actual size of the array.

procedure uclgvb (f77par, felem, nelems, nvals, buf, istat)

% 	character*(*) f77par
bool	buf[ARB]
int	felem			# first element in array
int	nelems			# number of elements
int	nvals			# mumber of elements actually retrieved
int	istat

char 	parnam[SZ_FNAME]
char 	parstr[SZ_FNAME]
char 	outstr[SZ_FNAME]
int	opbrkt, index, junk, stridx(), strlen()
bool	clgetb()
bool	value

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
	opbrkt = stridx (' ', parnam)
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

	   # get parameter as a string - If not found we have run out of
	   # elements in the list unless it is the first attempt in which
	   # case the parameter does not exist.
	   iferr  {
	       call clgstr (parnam, parstr, SZ_FNAME) 
	   }  then  {
                 if (nvals == 0)  
                     istat = ER_CLNOTFND
	             return
	   }
	   value = clgetb (parnam)
	   if (IS_INDEFR (value))  {
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
