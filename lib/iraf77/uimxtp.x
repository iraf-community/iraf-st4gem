include <imio.h>
include <iraf77.h>

define	MAX_OPNF	25
# UIMXTP -- Is the adaptation of IMTGETIM to the IRAF77 interface. 
#  Get the next image name from the image template.  
# Our job is to escape any [ in the filename part of the image
# name to avoid interpretation of these chars as image section characters by
# IMIO.  An extra function to IMTGETIM is the interpretation of the
# group template i.e 'file[*]' which this routine expand to all the
# groups in a multigroup STF data file.
#
# Phil Hodge,  8-Mar-1989  modified:  don't use "[" in stridx

procedure uimxtp (imt, f77out, istat)

pointer	imt			# image template descriptor
%	character*(*)  f77out   
int     istat			# error code return   

char    strbuf[SZ_FNAME,MAX_OPNF]
char	outstr[SZ_FNAME], extn[3]
char	image[SZ_FNAME], fulln[SZ_FNAME]
char	lbr_char
int	nchars
pointer	sp, buf, imid, immap()
int	fntgfnb(), imt_mapname(), strncmp()
errchk	fntgfnb
bool	first_time
int	i, j, ind, stridx()
data	first_time /true/
int	imtgrn[MAX_OPNF]	# holds the current group number for a file
int	gcount[MAX_OPNF]	# holds the number of groups for a file
int	imtbuf[MAX_OPNF]	# buffer to hold list descriptors
bool	flags[MAX_OPNF]		# false if more than 1 group per file
int	count			# actual number of lists open
common	/tpcom/imtbuf,flags,count,imtgrn,gcount
define  next_ 99
define  again_ 98

begin
	istat = ER_OK
	call smark (sp)
	call salloc (buf, SZ_PATHNAME, TY_CHAR)

	do i = 1, count {
	   if (imtbuf[i] == imt) {
	      first_time = flags[i]
	      j = i
	      break
	   }
	}
	  
	if (first_time) {
again_
	   # get the next name from the list
	   if (fntgfnb (imt, Memc[buf], SZ_PATHNAME) == EOF) {
	       outstr[1] = EOS
	       call sfree (sp)
	       istat = ER_EOF
	       return
	   }

	   nchars = imt_mapname (Memc[buf], strbuf[1,j] , SZ_FNAME)
	   call strcpy (strbuf[1,j], outstr, SZ_FNAME)
	   lbr_char = '['
	   ind = stridx (lbr_char, outstr)   # a possible '[*]'
	   call strcpy (outstr, fulln, SZ_FNAME)
	   fulln[ind] = EOS		# accept only "image.extn"
	   call iki_parse (fulln, image, extn)
	   if (extn[3] != 'h' || extn[1] == EOS)  # accept only '??h' extensions.
		goto again_

   	   # If the expanded filename has a wild card to indicate the
	   # expansion of the groups then open the image and get the
	   # number of groups.

	   if (strncmp ("[*]", outstr[ind], 3) == 0) {
	      iferr(imid = immap (fulln, READ_ONLY, 0)) {
	          istat = ER_IMOPOLD
		  return
	      }
	      gcount[j] = IM_CLSIZE(imid)
	      call imunmap(imid)
	      do i = 1, count {
	         if (imtbuf[i] == imt) {
	         flags[i] = false
	         break
	         }
	      }
	      imtgrn[j] = 1
	      goto next_
	   }
 	   call f77pak (outstr, f77out, SZ_FNAME)

	} else {
	   call strcpy (strbuf[1,j], outstr, SZ_FNAME)
	   lbr_char = '['
	   ind = stridx (lbr_char, outstr)   # a possible '[*]'
next_
	   call sprintf (outstr[ind+1], SZ_FNAME, "%d]" )
	       call pargi(imtgrn[j])
	   call sfree (sp)
   	   call f77pak (outstr, f77out, SZ_FNAME)
	   if (imtgrn[j] >= gcount[j]) {
	      do i = 1, count {
	         if (imtbuf[i] == imt) {
	         flags[i] = true    # Ready to read next item in list
	         break
	         }
	      }
	   } else 
	     imtgrn[j] = imtgrn[j] + 1
	}
	return
end
