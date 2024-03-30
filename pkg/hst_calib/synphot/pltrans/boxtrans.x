include	<gset.h>
include "../limit.h"

#* HISTORY *
#* B.Simon	22-Jul-94	original

procedure boxtrans (gp, inlimit, outlimit, xmode, ymode, xform, yform, title)

pointer	gp		# i: graphics descriptor
real	inlimit[4]	# i: input plot limits
real	outlimit[4]	# i: calculated plot limits
char	xmode[ARB]	# i: x axis mode
char	ymode[ARB]	# i: y axis mode
char	xform[ARB]	# i: x axis form
char	yform[ARB]	# i: y axis form
char	title[ARB]	# i: plot title
#--
int	i, junk, xmag, ymag, nexp
pointer sp, dir, mode1, mode2, xlabel, ylabel
real	dif, temp, limit[4]

string	dirlist  "left,right,bottom,top"
string  notauto  "Cannot autoscale plot, calculated limit is INDEF"
string	labelfmt "%s (%s)"

int	word_find(), is_magunit()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (dir, SZ_FNAME, TY_CHAR)
	call salloc (mode1, SZ_LINE, TY_CHAR)
	call salloc (mode2, SZ_LINE, TY_CHAR)
	call salloc (xlabel, SZ_LINE, TY_CHAR)
	call salloc (ylabel, SZ_LINE, TY_CHAR)

	# Adjust plot limits

	call amovr (outlimit, limit, 4)

	do i = 1, 4 {
	    if (IS_INDEFR(limit[i])) {
		junk = word_find (i, dirlist, Memc[dir], SZ_FNAME)
		call printerr_str (notauto, limit)
	    }
	}

	dif = limit[TOP] - limit[BOTTOM]
	if (dif == 0.0)
	    dif = limit[BOTTOM]

	if (IS_INDEFR(inlimit[BOTTOM]))
	    limit[BOTTOM] = limit[BOTTOM] - 0.05 * dif

	if (IS_INDEFR(inlimit[TOP]))
	    limit[TOP] = limit[TOP] + 0.05 * dif

	dif = limit[RIGHT] - limit[LEFT]
	if (dif == 0.0)
	    dif = limit[LEFT]

	if (IS_INDEFR(inlimit[LEFT]))
	    limit[LEFT] = limit[LEFT] - 0.05 * dif

	if (IS_INDEFR(inlimit[RIGHT]))
	    limit[RIGHT] = limit[RIGHT] + 0.05 * dif

	# Swap limits for magnitude units

	xmag = is_magunit (xform)
	if (xmag == YES) {
	    call splitexp (xmode, nexp, Memc[mode1], Memc[mode2], SZ_LINE)
	    if (nexp > 1)
		xmag = NO
	}

	if (xmag == YES) {
	    temp = limit[RIGHT]
	    limit[RIGHT] = limit[LEFT]
	    limit[LEFT] = temp
	}

	ymag = is_magunit (yform)

	if (ymag == YES) {
	    temp = limit[TOP]
	    limit[TOP] = limit[BOTTOM]
	    limit[BOTTOM] = temp
	}

	# Get axis labels

	call sprintf (Memc[xlabel], SZ_LINE, labelfmt)
	call pargstr (xmode)
	call pargstr (xform)

	call sprintf (Memc[ylabel], SZ_LINE, labelfmt)
	call pargstr (ymode)
	call pargstr (yform)

	# Set limits and plot axes

	call gswind (gp, limit[LEFT], limit[RIGHT], limit[BOTTOM], limit[TOP])

	call gseti (gp, G_TXQUALITY, GT_HIGH)
	call glabax (gp, title, Memc[xlabel], Memc[ylabel])

end
