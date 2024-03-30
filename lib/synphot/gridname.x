include <tbset.h>

# GRIDNAME -- Write names of files used in grid interpolation to log file

procedure gridname (list, place, log)

char	list[ARB]	# i: list of filenames
real	place		# i: place to interpolate at
int	log		# i: log file names are written on
#--
int	nfiles
pointer	sp, filelist

errchk	gridfiles

begin
	call smark (sp)
	call salloc (filelist, 2*(SZ_FNAME+1), TY_CHAR)

	call gridfiles (list, place, SZ_FNAME, 2, nfiles, Memc[filelist])

	# Print filenames

	call fprintf (log, "%s\n%s\n")
	call pargstr (Memc[filelist])
	call pargstr (Memc[filelist+SZ_FNAME+1])

	call sfree (sp)
end

# GRIDFILES -- Return names of files used in grid interpolation

procedure gridfiles (list, place, maxname, maxfile, nfiles, filelist)

char	list[ARB]		# i: list of filenames
real	place			# i: place to interpolate at
int	maxname			# i; maximum length of filename
int	maxfile			# i: maximum number of files
int	nfiles			# o: number of files
char	filelist[maxname,ARB]	# o: list of filenames
#--
int	ilo, ihi, nrow
pointer	 tp, cp

string	noextrap  "Cannot extrapolate from list of spectra"
string	nofiles   "Insufficient space to write file names"

int	tbpsta(), opnsyntab()
errchk	opnsyntab, syncolptr, tbeggt, synphoterr

begin
	# Check arguments

	tp = opnsyntab (list)
	nrow = tbpsta (tp, TBL_NROWS)

	ilo = int (place)
	ihi = min (ilo + 1, nrow)

	if (ilo < 1 || ilo > nrow)
	    call synphoterr (noextrap, "grid")

	if (maxfile < 2)
	    call synphoterr (nofiles, "grid")

	# Read spectra which bound interpolating point

	nfiles = 2

	call syncolptr (tp, "FILENAME", 1, cp)
	call tbegtt (tp, cp, ilo, filelist[1,1], maxname)
	call tbegtt (tp, cp, ihi, filelist[1,2], maxname)

end

