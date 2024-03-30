# PARSEDAT -- Parse the data string for fitspec and fitband

procedure parsedat(data, sphotlist, nsphot, mxspec, photlist, nphot, mxphot)

char	data[ARB]		# i: Data string = "s= sfile p= pfile"
char	sphotlist[SZ_LINE,ARB]	# o: specphot data files
int	nsphot			# o: number of specphot files
int	mxspec			# i: maximum number of specphot files
char	photlist[SZ_LINE,ARB]	# o: photometry data files
int	nphot			# o: number of phot files
int	mxphot			# i: maximum number of phot files
#--

int	ip, word
int	strsearch(), ctowrd()
pointer	sp, pfile, spfile
char	junk[SZ_LINE]

string	nodata		"No data to fit in BANDCHI2"

begin

	# Allocate memory
	call smark ( sp )
	call salloc( pfile, SZ_FNAME, TY_CHAR )
	call salloc( spfile, SZ_FNAME, TY_CHAR )

	call strlwr( data )

	# Parse the data string
	ip = 1
	for(word=1; ctowrd( data, ip, junk, SZ_FNAME ) > 0; word=word+1 ) {
	   if( word == 2 )
	      call strcpy( junk, Memc[spfile], SZ_FNAME)
	   else if (word == 4 )
	      call strcpy( junk, Memc[pfile], SZ_FNAME)
	   else if (word > 4 )
	      call printf( "Too many words in DATA string, continuing\n" )
	}

	if ( strsearch( Memc[spfile], "none" ) <= 0 )
	   call loadlist( Memc[spfile], nsphot, sphotlist, mxspec )

	if (strsearch(Memc[pfile], "none") <= 0)
	   call loadlist( Memc[pfile], nphot, photlist, mxphot )

	if (nsphot <= 0 && nphot <= 0 )
	   call error( 1, nodata )

end

