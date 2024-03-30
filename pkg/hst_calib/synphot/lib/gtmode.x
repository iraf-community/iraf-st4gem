define	SZ_EXTN		9

# GTMODE -- Get the next instrument mode from a text file
#
# The mode string passed as input may either contain a file name, indicated
# by a leading '@', or an instrument configuration mode. If the string contains
# a file name, each line the file contains a different instrument configuration
# mode. The grand throughput file name string is interpreted as a file name
# root and a suffix is added to to make the name unique. If the mode string
# is not a file, it and the grand throughput filename string are copied
# to the output variables as is. The variable `index' should be set to zero
# the first time this procedure is called. Each time the procedure is called
# index is incremented, unless no more lines remain in the file, in which
# case it is reset to zero. If the mode string is not a file, the second
# time it is called it is reset to zero, as there are no more instrument
# modes to process.
#
# Written by B.Simon on 6 June 88

procedure gtmode (modstr, gndstr, mode, gndtbl, index)

#	Input:
				# String containing mode file name
%	character*(*)	modstr
				# String containing grand thruput name root
%	character*(*)	gndstr

#	Output:
				# Instrument mode extracted from file
%	character*(*)	mode
				# Grand thruput table name
%	character*(*)	gndtbl

#	Input/Output:
				# Count on number of modes processed so far
int	index
#--
char	nl
int	fd, ic, junk
pointer	sp, modfile, sppmode, ldir, root, extn, gndfile

data	nl	/ '\n' /

int	open(), getline(), fnldir(), fnroot(), fnextn(), stridx()

begin
	# Allocate dynamic memory

	call smark (sp)
	call salloc (modfile, SZ_FNAME, TY_CHAR)
	call salloc (sppmode, SZ_LINE, TY_CHAR)
	call salloc (ldir, SZ_FNAME, TY_CHAR)
	call salloc (root, SZ_FNAME, TY_CHAR)
	call salloc (extn, SZ_EXTN, TY_CHAR)
	call salloc (gndfile, SZ_FNAME, TY_CHAR)

#	If the first character in the mode string is an '@'
#	the mode string is a file and each line of the file
#	is an instrument mode
	
	call f77upk (modstr, Memc[modfile], SZ_FNAME)
	if (Memc[modfile] == '@') {

	    # The first time this procedure is called,
	    # index equals 0. The file must be opened.
	    # On subsequent calls, one line is read from
	    # the mode file until end of file, when the
	    # file is closed.

	    if (index == 0)
		fd = open (Memc[modfile+1], READ_ONLY, TEXT_FILE)
	    if (getline (fd, Memc[sppmode]) != EOF) {
		index = index + 1

		ic = stridx (nl, Memc[sppmode])
		if (ic > 0)
		    Memc[sppmode+ic-1] = EOS
		call f77pak (Memc[sppmode], mode, SZ_LINE)

		call f77upk (gndstr, Memc[gndfile], SZ_FNAME)
		junk = fnldir (Memc[gndfile], Memc[ldir], SZ_FNAME)
		junk = fnroot (Memc[gndfile], Memc[root], SZ_FNAME)
		if (fnextn (Memc[gndfile], Memc[extn], SZ_EXTN) == 0)
		    call strcpy ("tab", Memc[extn], SZ_EXTN)

		call sprintf (Memc[gndfile], SZ_FNAME, "%s%s%03d.%s")
		    call pargstr (Memc[ldir])
		    call pargstr (Memc[root])
		    call pargi (index)
		    call pargstr (Memc[extn])
		call f77pak (Memc[gndfile], gndtbl, SZ_FNAME)

	    } else {

		index = 0
%		mode = ' '
%		gndtbl = ' '
		call close (fd)

	    }
	} else {

	    # If the mode string is not a file name, copy
	    # to output only on the first time through

	    if (index == 0) {
		index = 1
%		mode = modstr
%		gndtbl = gndstr
	    } else {
		index = 0
%		mode = ' '
%		gndtbl = ' '
	    }
	}

	call sfree (sp)
	return

end
