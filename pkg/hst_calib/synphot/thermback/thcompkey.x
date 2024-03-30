include <tbset.h>

define	MAXLEN   8  # max length of thtype keyword
define  SZ_COMPID       18              # Length of component name
#include "libsynphot.h": convention precludes including this

# GET_TCOMP_KEYWORDS -- Get some keywords from the header of a thermal component file

procedure thcompkey(tfile, thlmode, Bff,temp,  thtype)     
           # Open thermal component table (filename tflist[icomp])
           # Obtain component name in the COMPNAME keyword
           # Find match to the mode string (in the THTYPE keyword)
           # Look up beam filling fraction (in corresponding BEAMFILL keyword)


char	tfile[ARB]  # i: the file
char	thlmode[SZ_COMPID]   # o: the component "thermal mode". This is the
			     # keyword that can be specified inthe obsmode for this component.
real	Bff		  # o: the visual area of the component
real	temp		  # o: the default temperature of the component
char    thtype[MAXLEN]	  # o: thermal type (legal values CLEAR, OPAQUE, THRU)


pointer	opnsyntab()
errchk	opnsyntab

#pointer tbtopen()

pointer	tp
real    tbhgtr()

begin

	# Extract table and column name for file name and open table

#	call breakcomp (tfile, Memc[tabname], Memc[thruname], SZ_FNAME)  -- DO LATER

	tp = opnsyntab(tfile)
#	tp = tbtopen(tfile, READ_ONLY, NULL)

	Bff = tbhgtr(tp, "BEAMFILL")
	temp = tbhgtr(tp, "DEFT")

	call tbhgtt(tp, "THTYPE", thtype, MAXLEN) 
	call strupr(thtype)

	call tbhgtt(tp, "THMODE", thlmode, SZ_COMPID)

# do we need to close the table? We'd rather not, actually.
# Will the caching things be clever enough to figure out that it's already open when we come back?

end

