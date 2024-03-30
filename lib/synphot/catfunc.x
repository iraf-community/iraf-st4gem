include	<mach.h>
include <tbset.h>
define	MAXARGS		10
define	MAXSPEC		1024
define	INDEXSTR	"INDEX"
define	FILESTR		"FILENAME"

# CATFUNC -- Read a spectrum from a catalog which most closely matches

procedure catfunc (index, catfile, nwave, wave, spec)

char	index[ARB]	# i: index string used to search catalog
char	catfile[ARB]	# i: catalog filename template
int	nwave		# i: length of wavelength and spectrum arrays
real	wave[ARB]	# i: wavelengths at which spectrum is computed
real	spec[ARB]	# o: spectrum from catalog
#--
int	nmatch, match
pointer	tp

pointer	catopen()

begin
	# Open catalog file

	tp = catopen (index, catfile)

	# Find best match in catalog to argument 

	call catmatch (tp, index, 1, nmatch, match)

	# Read spectrum which matches index

	call catinterp (tp, index, nmatch, match, nwave, wave, spec)
			
end

# ICATFUNC -- Interpolate between spectra in a catalog

procedure icatfunc (index, catfile, nwave, wave, spec)

char	index[ARB]	# i: index string used to search catalog
char	catfile[ARB]	# i: catalog filename template
int	nwave		# i: length of wavelength and spectrum arrays
real	wave[ARB]	# i: wavelengths at which spectrum is computed
real	spec[ARB]	# o: spectrum from catalog
#--
int	nmatch
pointer	tp, sp, match

pointer	catopen()

begin
	# Allocate memory for temporary array

	call smark (sp)
	call salloc (match, MAXSPEC, TY_INT)

	# Open catalog file

	tp = catopen (index, catfile)

	# Find best match in catalog to argument 

	call catmatch (tp, index, MAXSPEC, nmatch, Memi[match])

	# Read spectrum which matches index

	call catinterp (tp, index, nmatch, Memi[match], nwave, wave, spec)

	call sfree (sp)
end

# CATNAME -- Write  names of spectra used in interpolation to a file

procedure catname (index, catfile, interp, log)

char	index[ARB]	# i: index string used to search catalog
char	catfile[ARB]	# i: catalog filename template
int	interp		# i: return entire bracket (YES) or closet file (NO)
int	log		# i: log file names are written on
#--
int	ifile, nfiles
pointer	sp, filelist, fname

begin
	# Allocate memory for temporary array

	call smark (sp)
	call salloc (filelist, MAXSPEC*(SZ_FNAME+1), TY_CHAR)

	call catfiles (index, catfile, interp, SZ_FNAME, 
		       MAXSPEC, nfiles, Memc[filelist])

	fname = filelist
	do ifile = 1, nfiles {
	    call fprintf (log, "%s\n")
	    call pargstr (Memc[fname])
	    fname = fname + SZ_FNAME + 1
	}

	call sfree (sp)
end

# CATFILES -- Return a list of files used in interpolation

procedure catfiles (index, catfile, interp, maxname, maxfile, nfiles, filelist)

char	index[ARB]		# i: index string used to search catalog
char	catfile[ARB]		# i: catalog filename template
int	interp			# i: return bracket (YES) or closest file (NO)
int	maxname			# i; maximum length of filename
int	maxfile			# i: maximum number of files
int	nfiles			# o: number of files
char	filelist[maxname,ARB]	# o: list of filenames
#--
int	nc, nmatch, ifiles
pointer	tp, cp, sp, catalog, path, match

string	filecol   FILESTR

int	fnldir(), gstrcpy()
pointer	catopen()

begin
	# Allocate memory for temporary array

	call smark (sp)
	call salloc (catalog, SZ_PATHNAME, TY_CHAR)
	call salloc (path, SZ_PATHNAME, TY_CHAR)
	call salloc (match, maxfile, TY_INT)

	# Open catalog file

	tp = catopen (index, catfile)

	# Find best match in catalog to argument 

	if (interp == NO) {
	    nmatch = 1
	} else {
	    nmatch = maxfile
	}

	call catmatch (tp, index, nmatch, nfiles, Memi[match])

	# Eliminate duplicates from list of filenames

	call catuniq (nfiles, Memi[match])

	# Extract directory name from catalog table

	call syntabname (tp, Memc[catalog], SZ_PATHNAME)
	nc = fnldir (Memc[catalog], Memc[path], SZ_PATHNAME)
	
	# Read catalog table and write matched filenames to output array

	call syncolptr (tp, filecol, 2, cp)

	do ifiles = 1, nfiles {
	    nc = gstrcpy (Memc[path], filelist[1,ifiles], maxname)
	    call tbegtt (tp, cp, Memi[match+ifiles-1], 
			 filelist[nc+1,ifiles], maxname)
	}

	call sfree (sp)
end

#--------------------------------------------------------------#
# The remaining routines are private to catfunc and icatfunc   #
#--------------------------------------------------------------#

# CATARGS -- Extract arguments from catalog argument string

procedure catargs (index, maxargs, maxch, argbuf, argname, argval, marg, narg)

char	index[ARB]	# i: index string used to search catalog
int	maxargs		# i: maximum number of arguments
int	maxch		# i: size of argument buffer
pointer	argbuf		# i: ptr to name buffer
pointer	argname[ARB]	# o: ptrs to argument names
real	argval[ARB]	# o: values of numeric arguments
int	marg		# o: number of numeric arguments
int	narg		# o: total number of arguents
#--
int	ic, jc, kc, mc, nc, iarg

string	toofew   "Arguments to cat function missing"
string	toomany  "Too many arguments to cat function"

int	word_fetch(), word_count(), ctor()

begin
	# Skip first argument, which is catalog name

	ic = 1
	nc = word_fetch (index, ic, Memc[argbuf], maxch)

	# Extract arguments from index string. Convert numeric arguments

	narg = word_count (index[ic])
	if (narg <= 0)
	    call synphoterr (toofew, index)

	if (narg > maxargs)
	    call synphoterr (toomany, index)

	jc = 0
	iarg = 1
	marg = 0
	repeat {
	    nc = word_fetch (index, ic, Memc[argbuf+jc], maxch-jc)
	    call strlwr (Memc[argbuf+jc])
	    if (nc == 0)
		break

	    kc = 1
	    argname[iarg] = argbuf +jc
	    mc = ctor (Memc[argbuf+jc], kc, argval[iarg])
	    if (mc != nc) {
		argval[iarg] = INDEFR
	    } else {
		marg = marg + 1
	    }

	    jc = jc + nc + 1
	    iarg = iarg + 1
	}

end

# CATDIST -- Compute the distance between the index and the argument list

real procedure catdist (index, narg, argname, argval, argsize)

char	index[ARB]	# i: list of catalog parameters
int	narg		# i: number of arguments
pointer	argname[narg]	# i: list of argument strings
real	argval[narg]	# i: list of argument values
real	argsize[narg]	# i: argument scaling factors
#--
int	ic, jc, nc, iarg
pointer	sp, name
real	value, dist, diff

string	errstr     "cat"
string	badargnum  "Incorrect number of arguments for function"
string	badargtyp  "Incorrect argument type for function"
string	badargsize "Bad argument scaling factor"

bool	strne()
int	word_fetch(), ctor()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (name, SZ_FNAME, TY_CHAR)

	# Strings must match exactly, otherwise relative squared distance
	# between catalog entry and argument list is evaluated

	ic = 1
	iarg = 1
	dist = 0.0
	while (word_fetch (index, ic, Memc[name], SZ_FNAME) > 0) {
	    if (iarg > narg)
		call synphoterr (badargnum, errstr)

	    # An indefinite value in the argument array signals
	    # that this is a string argument, which must match exactly

	    if (IS_INDEFR (argval[iarg])) {
		call strlwr (Memc[name])

		if (strne (Memc[name], Memc[argname[iarg]])) {
		    dist = MAX_REAL
		    break
		}

	    } else {
		jc = 1
		nc = ctor (Memc[name], jc, value)
		if (Memc[name+jc-1] != EOS)
		    call synphoterr (badargtyp, errstr)

		if (argsize[iarg] <= 0.0)
		    call synphoterr (badargsize, errstr)

		diff = (argval[iarg] - value) / argsize[iarg]
		dist = dist + diff * diff
	    }

	    iarg = iarg + 1
	}

	call sfree (sp)
	return (dist)
end

# CATINTERP -- Interpolate in a catalog of spectra

procedure catinterp (tp, index, nmatch, match, nwave, wave, spec)

pointer	tp		# i: catalog table descriptor
char	index[ARB]	# i: index string used to search catalog
int	nmatch		# i: number of rows matched
int	match[ARB]	# i: row indices
int	nwave		# i: length of wavelength and spectrum arrays
real	wave[ARB]	# i: wavelengths at which spectrum is computed
real	spec[ARB]	# o: spectrum from catalog
#--
int	ispec, jspec, nspec, ic, jc, nc, iarg, jarg, marg, narg
pointer	sp, name, rowidx, path, catalog, argbuf, argname, argval
pointer	fname, specptr, argptr, list, rowargs
pointer	fcol, icol, minbrak, maxbrak, newbrak, minspec, maxspec, newspec

string	filecol   FILESTR
string	indexcol  INDEXSTR

int	fnldir(), ctor(), word_fetch()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (name, SZ_FNAME, TY_CHAR)
	call salloc (rowidx, SZ_FNAME, TY_CHAR)
	call salloc (path, SZ_PATHNAME, TY_CHAR)
	call salloc (catalog, SZ_PATHNAME, TY_CHAR)
	call salloc (argbuf, SZ_LINE, TY_CHAR)
	call salloc (argname, MAXARGS, TY_POINTER)
	call salloc (argval, MAXARGS, TY_REAL)

	# Extract arguments from index string. Convert numeric arguments

	call catargs (index, MAXARGS, SZ_LINE, argbuf, Memi[argname], 
		      Memr[argval], marg, narg)

	# Copy catalog directory into filename path

	call syntabname (tp, Memc[catalog], SZ_PATHNAME)
	fname = path + fnldir (Memc[catalog], Memc[path], SZ_PATHNAME)
	
	# Get column pointers

	call syncolptr (tp, indexcol, 1, icol)
	call syncolptr (tp, filecol, 2, fcol)

	# Check for case of a single spectrum

	if (nmatch == 1) {
	    call tbegtt (tp, fcol, match[1], Memc[fname], SZ_FNAME)
	    call rdspec (Memc[path], nwave, wave, spec)

	    call sfree (sp)
	    return
	}

	# Allocate interpolation arrays

	call salloc (specptr, nmatch, TY_POINTER)
	call salloc (argptr, nmatch, TY_POINTER)
	call salloc (list, nmatch, TY_INT)

	# Read matched rows from catalog table

	do ispec = 0, nmatch-1 {
	    # Read spectra to interpolate between

	    call salloc (Memi[specptr+ispec], nwave, TY_REAL)

	    call tbegtt (tp, fcol, match[ispec+1], Memc[fname], SZ_FNAME)
			     
	    call rdspec (Memc[path], nwave, wave, Memr[Memi[specptr+ispec]])
			     

	    # Read numeric arguments associated with spectra

	    call salloc (rowargs, marg, TY_REAL)
	    Memi[argptr+ispec] = rowargs

	    ic = 1
	    iarg = 0
	    call tbegtt (tp, icol, match[ispec+1], Memc[rowidx], SZ_FNAME)
			     
	    while (word_fetch (Memc[rowidx], ic, Memc[name], SZ_FNAME) > 0) {
		if ( ! IS_INDEFR (Memr[argval+iarg])) {
		    jc = 1
		    nc = ctor (Memc[name], jc, Memr[rowargs])
		    rowargs = rowargs + 1
		}

		iarg = iarg + 1
	    }
	}


	# Perform linear interpolation on successive pairs of spectra

	jarg = 0
	nspec = nmatch

	do iarg = 0, narg-1 {
	    if (IS_INDEFR(Memr[argval+iarg]))
		next

	    ispec = 0
	    jspec = 0
	    while (ispec < nspec - 1) {
		minbrak = Memi[argptr+ispec]
		maxbrak = Memi[argptr+ispec+1]
		newbrak = Memi[argptr+jspec]

		minspec = Memi[specptr+ispec]
		maxspec = Memi[specptr+ispec+1]
		newspec = Memi[specptr+jspec]

		call funinterp (Memr[minbrak+jarg], Memr[maxbrak+jarg],
				Memr[argval+iarg], nwave, Memr[minspec],
				Memr[maxspec], Memr[newspec])

		call funinterp (Memr[minbrak+jarg], Memr[maxbrak+jarg],
				Memr[argval+iarg], nmatch, Memr[minbrak],
				Memr[maxbrak], Memr[newbrak])

		ispec = ispec + 2
		jspec = jspec + 1
	    }

	    jarg = jarg + 1
	    nspec = jspec
	}

	call amovr (Memr[newspec], spec, nwave)
	call sfree (sp)
end

# CATLIST -- Compute a list of categories the index falls in

procedure catlist (index, narg, argname, argval, maxlist, nlist, list)

char	index[ARB]	# i: argument string associated with catalog spectrum
int	narg		# i: number of arguments
pointer	argname[narg]	# i: list of argument strings
real	argval[narg]	# i: list of argument values
int	maxlist		# i: maximum number of catgories
int	nlist		# o: number of catgories matched
int	list[maxlist]	# o: list of categories matched
#--
int	ilist, ic, jc, nc, inc, iarg
pointer	sp, name
real	value

string	errstr    "cat"
string	badargnum "Incorrect number of arguments for function"
string	badargtyp "Incorrect argument type for function"

bool	strne()
int	ctor(), word_fetch()

begin
	# Allocate temporary arrays

	call smark (sp)
	call salloc (name, SZ_FNAME, TY_CHAR)

	# There are two categories for each numeric argument. The first
	# is for all values less than argval[iarg] and the second is for 
	# all values greater than argval[iarg]. Values equal to iarg get 
	# placed in both categories. The categories are numbered so that 
	# the high catgeory is 2 ** (jarg-1) above the low category, where
	# jarg denotes the j-th numeric argument. Categories are numbered
	# starting at zero, which is the category for values less than all
	# of argval[].

	ic = 1
	inc = 1
	iarg = 1
	nlist = 1
	call aclrr (list, maxlist)

	while (word_fetch (index, ic, Memc[name], SZ_FNAME) > 0) {
	    if (iarg > narg)
		call synphoterr (badargnum, errstr)

	    # If two string arguments do not match, no category
	    # matches the catalog spectrum

	    if (IS_INDEFR(argval[iarg])) {
		call strlwr (Memc[name])
		if (strne (Memc[name], Memc[argname[iarg]])) {
		    nlist = 0
		    break
		}

	    } else {
		# If the catalog spectrum is less than the argument list
		# zero (implicitly) gets added to all values in the list.
		# If it is greater inc gets added. If it is equal the list
		# is doubled in length and zero gets added to half the values
		# and inc to the other half.

		jc = 1
		nc = ctor (Memc[name], jc, value)
		if (Memc[name+jc-1] != EOS)
		    call synphoterr (badargtyp, errstr)
		
		if (value > argval[iarg]) {
		    do ilist = 1, nlist 
			list[ilist] = list[ilist] + inc

		} else if (value == argval[iarg]) {
		    do ilist = nlist, 1, -1 {
			list[2*ilist] = list[ilist] + inc
			list[2*ilist-1] = list[ilist]
		    }

		    nlist = 2 * nlist
		}

		inc = 2 * inc
	    }

	    iarg = iarg + 1
	}

	call sfree (sp)

end

# CATMATCH -- Find rows in catalog which match or bracket index string

procedure catmatch (tp, index, maxmatch, nmatch, match)

pointer	tp		# i: catalog table descriptor
char	index[ARB]	# i: index into catalog
int	maxmatch	# i: maximum number of rows to match
int	nmatch		# o: number of rows matched
int	match[ARB]	# o: row indices
#--
int	marg, narg, irow, jrow, nrow, ilist, jlist, nlist, imatch
pointer	sp, argbuf, argname, argval, argsize, rowidx
pointer	minrow, mindist, list, cp
real	dist

string	errstr    "cat"	
string	indexcol  INDEXSTR
string  nomatch   "No matching spectrum found in catalog"
string	extrap    "Spectrum not bounded. Using nearest neighbor"
string	toomany	  "Too many arguments to function"

int	tbpsta()
real	catdist()

begin
	# Allocate memory for tempoarary strings

	call smark (sp)
	call salloc (argbuf, SZ_LINE, TY_CHAR)
	call salloc (argname, MAXARGS, TY_POINTER)
	call salloc (argval, MAXARGS, TY_REAL)
	call salloc (argsize, MAXARGS, TY_REAL)
	call salloc (rowidx, SZ_FNAME, TY_CHAR)

	# Extract arguments from index string. Convert numeric arguments

	call catargs (index, MAXARGS, SZ_LINE, argbuf, Memi[argname], 
		      Memr[argval], marg, narg)

	# Allocate interpolation arrays

	nmatch = 2 ** marg

	call salloc (minrow, nmatch, TY_INT)
	call salloc (mindist, nmatch, TY_REAL)
	call salloc (list, nmatch, TY_INT)

	call amovki (0, Memi[minrow], nmatch)
	call amovkr (MAX_REAL, Memr[mindist], nmatch)

	# Compute distance scale factors

	call catscale (tp, narg, Memr[argsize])

	# Search for spectra to use in interpolation

	nrow = tbpsta (tp, TBL_NROWS)
	call syncolptr (tp, indexcol, 1, cp)

	do irow = 1, nrow {
	    call tbegtt (tp, cp, irow, Memc[rowidx], SZ_FNAME)

	    dist = catdist (Memc[rowidx], narg, Memi[argname], 
			    Memr[argval], Memr[argsize])

	    if (dist == MAX_REAL)
		next

	    call catlist (Memc[rowidx], narg, Memi[argname], Memr[argval], 
			  nmatch, nlist, Memi[list])

	    do ilist = 0, nlist-1 {
		jlist = Memi[list+ilist]

		if (dist < Memr[mindist+jlist]) {
		    Memr[mindist+jlist] = dist
		    Memi[minrow+jlist] = irow
		}
	    }
	}

	# Check to see if spectrum is completely bounded
	# If it is not, print a warning message

	for (imatch = 0; imatch < nmatch; imatch = imatch + 1)
	    if (Memi[minrow+imatch] == 0)
		break

	if (imatch < nmatch)
	    call synphotwarn (extrap, errstr)

	# Collapse list of rows down to the best match 
	# if not bounded or only one match was sought

	if (imatch < nmatch || maxmatch == 1) {
	    jrow = 0
	    dist = MAX_REAL
	    do imatch = 0, nmatch-1 {
		if (Memr[mindist+imatch] < dist) {
		    dist = Memr[mindist+imatch]
		    jrow = Memi[minrow+imatch]
		}
	    }

	    if (jrow == 0)
		call synphoterr (nomatch, errstr)

	    nmatch = 1
	    match[1] = jrow

	} else {
	    if (nmatch > maxmatch)
		call synphoterr (toomany, errstr)

	    call amovi (Memi[minrow], match, nmatch)
	}

	call sfree (sp)
end

# CATOPEN -- Open the catalog table

pointer procedure catopen (index, catfile)

char	index[ARB]	# i: index string used to search catalog
char	catfile[ARB]	# i: catalog filename template
#--
char	star
int	ic, jc, nc
pointer	sp, dir, catalog, ext, tp

data	star	/ '*' /

string	extensions   "fits,tab"
string	nowild       "No wildcard character in catalog filename template"
string	nodirectory  "Catalog directory was not found"
string	nocatalog    "Catalog file was not found"

int	word_fetch(), stridx(), strlen(), access()
pointer	opnsyntab()

begin
	# Allocate memory for tempoarary strings

	call smark (sp)
	call salloc (dir, SZ_FNAME, TY_CHAR)
	call salloc (ext, SZ_FNAME, TY_CHAR)
	call salloc (catalog, SZ_PATHNAME, TY_CHAR)

	# Extract directory name from argument list

	ic = 1
	nc = word_fetch (index, ic, Memc[dir], SZ_FNAME)

	# Replace star in template with directory name

	ic = stridx (star, catfile)
	if (ic == 0)
	    call synphoterr (nowild, catfile)

	call strcpy (catfile, Memc[catalog], ic-1)
	call strcat (Memc[dir], Memc[catalog], SZ_PATHNAME)

	if (access (Memc[catalog], 0, 0) == NO)
	    call synphoterr (nodirectory, Memc[dir])

	call strcat (catfile[ic+1], Memc[catalog], SZ_PATHNAME)

	jc = strlen (Memc[catalog])
	Memc[catalog+jc] = '.'
	jc = jc + 1

	# Try different extensions for the file

	ic = 1
	tp = NULL
	while (word_fetch (extensions, ic, Memc[ext], SZ_FNAME) > 0) {
	    call strcpy (Memc[ext], Memc[catalog+jc], SZ_FNAME)

	    if (access (Memc[catalog], READ_ONLY, 0) == YES) {
		# Open catalog file
		tp = opnsyntab (Memc[catalog])
	    }
	}

	if (tp == NULL) 
	    call synphoterr (nocatalog, Memc[catalog])

	call sfree (sp)
	return (tp)
end


# CATSCALE -- Compute argument scaling factors

procedure catscale (tp, narg, argsize)

pointer	tp		# i: table descriptor
int	narg		# i: number of arguments
real	argsize[ARB]	# i: argument scaling factors
#--
int	ic, jc, nc, iarg, pass, irow, nrow
pointer	sp, name, index, avg, cp
real	value

string	errstr     "cat"
string	scalekey   "SCALE"
string	indexcol   INDEXSTR
string	badargnum  "Incorrect number of arguments for function"

int	word_fetch(), tbpsta(), ctor()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (name, SZ_FNAME, TY_CHAR)
	call salloc (index, SZ_FNAME, TY_CHAR)
	call salloc (avg, narg, TY_REAL)
	call amovkr (0.0, Memr[avg], narg)

	# Look in table header for the scale factors first

	ifnoerr {
	    call rdtabhdt (tp, scalekey, Memc[index], SZ_FNAME)

	} then {
	    ic = 1
	    iarg = 1
	    while (word_fetch(Memc[index], ic, Memc[name], SZ_FNAME) > 0) {
		if (iarg > narg)
		    call synphoterr (badargnum, errstr)

		jc = 1
		nc = ctor (Memc[name], jc, argsize[iarg])
		iarg = iarg + 1
	    }

	    call sfree (sp)
	    return
	}

	# Otherwise set the scale factors to the average deviation
	# from the mean

	nrow = tbpsta (tp, TBL_NROWS)
	call syncolptr (tp, indexcol, 1, cp)

	# The first pass calculates the mean
	# The second pass calculates the deviation from the mean

	do pass = 1, 2 {
	    call amovkr (0.0, argsize, narg)

	    do irow = 1, nrow {
		ic = 1
		iarg = 1
		call tbegtt (tp, cp, irow, Memc[index], SZ_FNAME)

		while (word_fetch(Memc[index], ic, Memc[name], SZ_FNAME) > 0) {
		    if (iarg > narg)
			call synphoterr (badargnum, errstr)

		    jc = 1
		    nc = ctor (Memc[name], jc, value)
		    if (Memc[name+jc-1] == EOS) {
			if (pass == 1) {
			    argsize[iarg] = argsize[iarg] + value
			} else {
			    argsize[iarg] = argsize[iarg] + 
					    abs (value - Memr[avg+iarg-1])
			}
		    }

		    iarg = iarg + 1
		}
	    }

	    call adivkr (argsize, real(nrow), argsize, narg)
	    call amovr (argsize, Memr[avg], narg)
	}

	call sfree (sp)
end

# CATUNIQ -- Eliminate duplicates in list of row numbers

procedure catuniq (nmatch, match)

int	nmatch		# u: number of rownumbers
int	match[ARB]	# u: row numbers
#--
int	imatch, jmatch

begin
	call asrti (match, match, nmatch)

	jmatch = 1
	do imatch = 2, nmatch {
	    if (match[imatch] != match[imatch-1])
		jmatch = jmatch + 1

	    match[jmatch] = match[imatch]
	}

	nmatch = jmatch
end
