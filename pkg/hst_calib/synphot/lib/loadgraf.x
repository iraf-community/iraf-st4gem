include	<tbset.h>
include	"../graferr.h"

# LOADGRAF -- Load contents of graph table(s) into arrays
#
# This procedure dynamically allocates arrays to hold the contents of a
# graph table. These arrays must be freed by the calling procedure. If
# a file name template is passed to this procedure, the graph tables which
# match the template are combined into a single graph.
#
# B.Simon	18-Jul-88	First Code

procedure loadgraf (graph, nrow, nkeyword, ncomp, innode, outnode,
		    compname, keyword)

char	graph[ARB]	# i: Name of graph file or template
int	nrow		# o: Number of rows in graph
int	nkeyword	# o: Length of keyword string
int	ncomp		# o: Length of component name string
pointer	innode		# o: Array of innode values
pointer	outnode		# o: Array of outnode values
pointer	compname	# o: Array of component names
pointer	keyword		# o: Array of mnemonic keywords
#--
int	irow, mrow
pointer	sp, tp, file, errmsg, cp, flag, record

string	norows  "Graph table is empty (%s)"

int	tbpsta(), tbcigi()
pointer	tbtopn()

errchk	tbtopn, tbcgti, tbcgtt

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (file, SZ_FNAME, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	# Open the table and get the number of rows

	call lastfile (graph, Memc[file], SZ_FNAME)
	tp = tbtopn (Memc[file], READ_ONLY, NULL)
	nrow = tbpsta (tp, TBL_NROWS)

	# Print error message if no rows in table

	if (nrow == 0) {
	    call sprintf (Memc[errmsg], SZ_LINE, norows)
	    call pargstr (graph)
	    call error (SYNTAX, Memc[errmsg])
	}

	# Get the length of columns containing strings

	call tbcfnd (tp, "COMPNAME", cp, 1)
	ncomp = - tbcigi (cp, TBL_COL_DATATYPE)
	call tbcfnd (tp, "KEYWORD", cp, 1)
	nkeyword = - tbcigi (cp, TBL_COL_DATATYPE)

	# Allocate arrays to hold data

	mrow = nrow + 1
	call malloc (innode, mrow, TY_INT)
	call malloc (outnode, mrow, TY_INT)
	call malloc (compname, mrow*(ncomp+1), TY_CHAR)
	call malloc (keyword, mrow*(nkeyword+1), TY_CHAR)
	call malloc (flag, nrow, TY_BOOL)

	# Read the table columns into the arrays

	call tbcfnd (tp, "INNODE", cp, 1)
	call tbcgti (tp, cp, Memi[innode], Memb[flag], 1, nrow)

	call tbcfnd (tp, "OUTNODE", cp, 1)
	call tbcgti (tp, cp, Memi[outnode], Memb[flag], 1, nrow)

	call tbcfnd (tp, "COMPNAME", cp, 1)
	call tbcgtt (tp, cp, Memc[compname], Memb[flag], ncomp, 1, nrow)

	call tbcfnd (tp, "KEYWORD", cp, 1)
	call tbcgtt (tp, cp, Memc[keyword], Memb[flag], nkeyword, 1, nrow)

	# Justify all strings by removing whitespace and 
	# converting to lower case

	record = compname
	do irow = 1, nrow {
	    call strjust (Memc[record], Memc[record], ncomp)
	    record = record + ncomp + 1
	}

	record = keyword
	do irow = 1, nrow {
	    call strjust (Memc[record], Memc[record], nkeyword)
	    record = record + nkeyword + 1
	}

	# Close table and release memory

	call tbtclo (tp)
	call mfree (flag, TY_INT)
	call sfree (sp)
end
