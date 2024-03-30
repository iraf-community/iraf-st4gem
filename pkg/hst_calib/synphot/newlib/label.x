include	<gset.h>

#* HISTORY *
#* B.Simon	09-Jun-94	original

# LABELI -- Print label containing integer

procedure labeli (gp, label, data)

pointer	gp		# i: graphics descriptor
char	label[ARB]	# i: text label
int	data		# i: data following label
#--
pointer	sp, label2

begin
	if (IS_INDEFI(data))
	    return

	call smark (sp)
	call salloc (label2, SZ_LINE, TY_CHAR)

	call sprintf (Memc[label2], SZ_LINE, " %d")
	call pargi (data)

	call labelz (gp, label, Memc[label2])
	call sfree (sp)
end

# LABELR -- Print label containing real number

procedure labelr (gp, label, data)

pointer	gp		# i: graphics descriptor
char	label[ARB]	# i: text label
real	data		# i: data following label
#--
pointer	sp, label2

begin
	if (IS_INDEFR(data))
	    return

	call smark (sp)
	call salloc (label2, SZ_LINE, TY_CHAR)

	call sprintf (Memc[label2], SZ_LINE, " %g")
	call pargr (data)

	call labelz (gp, label, Memc[label2])
	call sfree (sp)
end

# LABELT -- Print label containing text string

procedure labelt (gp, label, data)

pointer	gp		# i: graphics descriptor
char	label[ARB]	# i: text label
char	data[ARB]	# i: data following label
#--
pointer	sp, label2
bool	strne()

begin
	call smark (sp)
	call salloc (label2, SZ_LINE, TY_CHAR)

	call strcpy (data, Memc[label2], SZ_LINE)
	call strfix (Memc[label2])

	if (Memc[label2] != EOS && strne (Memc[label2], "none")) {
	    call sprintf (Memc[label2], SZ_LINE, " %s")
	    call pargstr (data)

	    call labelz (gp, label, Memc[label2])
	}

	call sfree (sp)
end

# LABELZ -- Low level procedure to print text label and move cursor

procedure labelz (gp, label1, label2)

pointer	gp		# i: graphics descriptor
char	label1[ARB]	# i: first text label
char 	label2[ARB]	# i: second text label
#--
real	xcur	# x cursor position
real	ycur	# y cursor position
real	yoff	# y cursor offset

data	xcur  / INDEFR /
data	ycur  / INDEFR /

pointer	sp, text
real	x1, x2, y1, y2
real	ggetr()

string	blank   ""

begin
	# Set cursor position if unset

	if (IS_INDEFR (xcur) || IS_INDEFR(ycur)) {
	    call ggview (gp, x1, x2, y1, y2)
	    call gseti (gp, G_WCS, 0)
	    xcur = x1 + 0.08 * (x2 - x1)
	    ycur = y2 - 0.09 * (y2 - y1)
	    yoff = 1.15 * ggetr (gp, "ch")
	}

	# Concatenate text strings

	call smark (sp)
	call salloc (text, SZ_LINE, TY_CHAR)

	call strcpy (label1, Memc[text], SZ_LINE)
	call strcat (label2, Memc[text], SZ_LINE)

	# Print text label and update cursor position

	ycur = ycur - yoff

	call gseti (gp, G_WCS, 0)
	call gseti (gp, G_TXQUALITY, GT_HIGH)
	call gtext (gp, xcur, ycur, Memc[text], blank)

	call sfree (sp)
	return

	# Unset cursor position

	entry labelu ()

	xcur = INDEFR
	ycur = INDEFR
	return

end
