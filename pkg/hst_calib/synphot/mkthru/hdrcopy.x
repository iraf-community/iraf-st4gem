include <time.h>
include <tbset.h>

# HDRCOPY -- Copy the table header to the primary hdu

procedure hdrcopy (ifile, ofile)

char	ifile[ARB]	# i: input file name
char	ofile[ARB]	# i: output file name
#--
int	ic
long	aeon
pointer	sp, header, keyword, value, date, itp, otp

data	aeon  / 0 /

string	keylist  "instrume,compname,useafter,pedigree,descrip,comment"

int	word_fetch()
long	clktime()
pointer	tbtopn()

errchk	tbtopn, tbhadt

begin
	# Allocate memory for strings

	call smark (sp)
	call salloc (header, SZ_FNAME, TY_CHAR)
	call salloc (keyword, SZ_KEYWORD, TY_CHAR)
	call salloc (value, SZ_PARREC, TY_CHAR)
	call salloc (date, SZ_DATE, TY_CHAR)

	# Open input table and primary header of output file

	call strcpy (ofile, Memc[header], SZ_FNAME)
	call strcat ("[0]", Memc[header], SZ_FNAME)

	itp = tbtopn (ifile, READ_ONLY, NULL)
	otp = tbtopn (Memc[header], READ_WRITE, NULL)

	# Set dbtable keyword

	call tbhadt (otp, "DBTABLE", "CRTHROUGHPUT")

	# Set other keywords from input table or get from user if missing

	ic = 1
	while (word_fetch (keylist, ic, Memc[keyword], SZ_KEYWORD) > 0) {
	    iferr {
		call tbhgtt (itp, Memc[keyword], Memc[value], SZ_PARREC)
	    } then {
		call clgstr (Memc[keyword], Memc[value], SZ_PARREC)
	    }

	    call strupr (Memc[keyword])
	    call tbhadt (otp, Memc[keyword], Memc[value])
	}

	# Write history record with date

	call cnvdate (clktime (aeon), Memc[date], SZ_DATE)

	call sprintf (Memc[value], SZ_FNAME, "Converted to fits from %s on %s")
	call pargstr (ifile)
	call pargstr (Memc[date])

	call tbhadt (otp, "HISTORY", Memc[value])

	# Close tables

	call tbtclo (itp)
	call tbtclo (otp)

	call sfree (sp)
end

