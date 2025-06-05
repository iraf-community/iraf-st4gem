# MAKE_IMNAME	build image name from root + index
# NOTE: Throw away code

procedure make_imname (rootname, index, count, imname)

char	rootname[SZ_PATHNAME]	# root file name
int	index			# the image number
int	count			# the group number
char	imname[SZ_PATHNAME]	# receives image name

int 	junk, itoc()
char	charindex[SZ_FNAME]	
char	charcount[SZ_FNAME]	

begin
	call strcpy (rootname, imname, SZ_PATHNAME)
	#
	junk = itoc (index, charindex, SZ_FNAME)	
	#
	call strcat ("[", imname, SZ_PATHNAME)
	call strcat (charindex, imname, SZ_PATHNAME)
	# see if name is of the form name[i] or name [i/j]
	if (count >= 0 )  {
	   # build group count number into file name
	   call strcat ("/", imname, SZ_PATHNAME)
	   junk = itoc (count, charcount, SZ_FNAME)	
	   call strcat (charcount, imname, SZ_PATHNAME)
	}
	call strcat ("]", imname, SZ_PATHNAME)
	return

end
