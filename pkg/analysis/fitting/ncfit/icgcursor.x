include	"icfit.h"

# ICG_CURSOR -- Routine to take care of cursor reading. It's necessary
# because of the y magnitude scale option, which is not supported by
# gtools.

int procedure icg_cursor (ic, gt, cursor, wx, wy, wcs, key, cmd, siz) 

pointer	gt, ic
char	cursor[ARB], cmd[SZ_LINE]
real	wx, wy
int	wcs, key, siz

int	i

int	gt_gcur1()

begin
	i = gt_gcur1 (gt, cursor, wx, wy, wcs, key, cmd, siz)
	if (IC_YAXIS(ic) == IC_MAG) {
	    wy = 10. ** ((IC_MAG0(ic) - wy) / 2.5)
	}
	return (i)
end
