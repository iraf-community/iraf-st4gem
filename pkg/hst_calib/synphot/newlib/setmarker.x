include	<gset.h>

define	MAXSYM   16
define	MARKSTR  "clear,solid,dashed,dotted,dotdash,point,fill,box,plus,\
cross,diamond,hline,vline,hebar,vebar,circle"

#* HISTORY *
#* B.Simon	10-Jun-94	Original
#* B.Simon	22-Jul-94	Added getmarker and linemarker

int procedure getmarker (mark)

char	mark[ARB]	# i: Marker name
#--
int	mtype, symbol[MAXSYM]
pointer	sp, umark

string	marklist MARKSTR
string	badmark  "Unrecognized marker type"

data	symbol   / GL_CLEAR, GL_SOLID, GL_DASHED, GL_DOTTED, GL_DOTDASH,
		   GM_POINT, GM_FILL, GM_BOX, GM_PLUS, GM_CROSS, 
	           GM_DIAMOND, GM_HLINE, GM_VLINE, GM_HEBAR, GM_VEBAR,
	           GM_CIRCLE /

int	word_match()

begin
	# Convert marker name to standard form

	call smark (sp)
	call salloc (umark, SZ_LINE, TY_CHAR)

	call strcpy (mark, Memc[umark], SZ_LINE)
	call strfix (Memc[umark])

	# Get corresponding symbol and set the mark type

	mtype = word_match (Memc[umark], marklist)
	if (mtype == 0 || mtype > MAXSYM)
	    call printerr_str (badmark, mark)

	call sfree (sp)
	return (symbol[mtype])

end

# LINEMARKER -- Return YES if marker is a line type

int procedure linemarker (mark)

char	mark[ARB]	# i: Marker name
#--
int	mtype, ltype, isline[MAXSYM]
pointer	sp, umark

string	marklist MARKSTR

data	isline   / 5 * YES, 11 * NO /

int	word_match()

begin
	# Convert marker name to standard form

	call smark (sp)
	call salloc (umark, SZ_LINE, TY_CHAR)

	call strcpy (mark, Memc[umark], SZ_LINE)
	call strfix (Memc[umark])

	# Get corresponding symbol and set the mark type

	mtype = word_match (Memc[umark], marklist)
	if (mtype == 0 || mtype > MAXSYM) {
	    ltype = NO
	} else {
	    ltype = isline[mtype]
	}

	call sfree (sp)
	return (ltype)

end

# SETMARKER -- Set the line or marker type for a plot

procedure setmarker (gp, mark)

pointer	gp		# i: Graphics descriptor
char	mark[ARB]	# i: Marker name
#--
int	symbol
int	getmarker()

begin
	symbol = getmarker (mark)
	call gseti (gp, G_PLTYPE, symbol)
end
