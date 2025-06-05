include <iraf77.h>
define  MAX_OPNF	25		# Maximun open templates.

# UIMCTP-- Close an image template.

procedure uimctp (imt, istat)

pointer	imt			# image template descriptor
int	i, istat
int	imtbuf[MAX_OPNF], imtgrn[MAX_OPNF], gcount[MAX_OPNF]
bool	flags[MAX_OPNF]
int	count

common	/tpcom/imtbuf,flags,count,imtgrn,gcount

begin
	istat = ER_OK
	#regroup pointer information in buffer
	do i = 1, count {
	   if (imtbuf[i] == imt) {
	      imtbuf[i] = imtbuf[count]
	      imtgrn[i] = imtgrn[count]
	      gcount[i] = gcount[count]
	      flags[i]  = flags[count]
	      count = count - 1
	   }
	}
	iferr(call fntclsb (imt))
	   istat = ER_IMTCLST
end
