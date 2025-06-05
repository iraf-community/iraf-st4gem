include	<syserr.h>
include	<imhdr.h>
include <imio.h>
include <mach.h>
include	<iraf77.h>

# IMADCF -- Add a user field to the image header.  It is an error if the named
# field already exists. 
# **If a group parameter block for an STF format is
# requested then put it on the gpb section
# of the datafile too.**NOT IMPLEMENTED**

procedure imadcf (im, key, lenv, datatype, comment, htype)

pointer	im			# image descriptor
char	key[ARB]		# name of the new parameter
int	lenv			# keyword value length
int	datatype		# datatype of parameter
char	comment[ARB]		# comment describing new parameter
int     htype                   # type of header parameter

int	max_lenuserarea, diff
pointer	sp, keyname, rp, ua, op
int	idb_kwlookup(), idb_findrecord(), strlen()
errchk	syserrs

begin
	call smark (sp)
	call salloc (keyname, SZ_FNAME, TY_CHAR)
	call strcpy (key, Memc[keyname], SZ_FNAME)

	# Check for a redefinition.
	if ((idb_kwlookup (key) > 0) || (idb_findrecord (im, key, rp) > 0))
	    call syserrs (SYS_IDBREDEF, key)
	
	# Open the user area string for appending.  If the user area is not
	# empty the last character must be the newline record delimiter,
	# else the new record we add will be invalid.

	max_lenuserarea = (LEN_IMDES + IM_LENHDRMEM(im) - IMU + 1) * SZ_STRUCT
	ua = IM_USERAREA(im)

	for (rp=ua;  Memc[rp] != EOS;  rp=rp+1)
	    ;
	if (rp - ua + SZ_VALSTR + 1 >= max_lenuserarea)
	    call syserrs (SYS_IDBOVFL, key)

	if (rp > ua && Memc[rp-1] != '\n') {
	    Memc[rp] = '\n'
	    rp = rp + 1
	}

	# Append the new record with an uninitialized value field.  Keyword
	# value pairs are encoded in FITS format.

	do op = rp, rp + SZ_VALSTR		# blank fill card
	    Memc[op] = ' '

	# Add the "= 'value' / comment".
	call amovc (Memc[keyname], Memc[rp], strlen(Memc[keyname]))
	Memc[rp+9-1] = '='
	if (datatype == TY_CHAR) {
	    Memc[rp+11-1] = '\''
	    Memc[rp+20-1] = '\''
	}

	# Add the comment field.
	diff = lenv - 20
	if (diff > 0) {
	   Memc[rp+34-1+diff] = '/'
	   call amovc (comment, Memc[rp+36-1+diff],
	       min (SZ_VALSTR-36+1-diff, strlen(comment)))
	} else {
	   Memc[rp+32-1] = '/'
	   call amovc (comment, Memc[rp+34-1],
	       min (SZ_VALSTR-34+1, strlen(comment)))
	}
	# Terminate the card.
	Memc[rp+SZ_VALSTR] = '\n'
	Memc[rp+SZ_VALSTR+1] = EOS

	# see if we need to add the parameter to the gpb section of the 
	# data file. **NOT IMPLEMENTED**

	call sfree (sp)
end
