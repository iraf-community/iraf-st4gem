# GTPHOT -- Fortran interface to getphotx

procedure gtphot (modstr, grftbl, cmptbl, path, phot)

%       character*(*)   modstr
%       character*(*)   grftbl
%       character*(*)   cmptbl
%       character*(*)   path
real	phot[4]			# photometric parameters (output)
#--
int	mxpath, mxmode
pointer	sp, mode, pathx, graphtab, comptab

begin
	# Get length of fortran strings

	mxpath = len(path)
	mxmode = len(modstr)

	# Allocate memory for temporary spp strings

	call smark (sp)
	call salloc (mode, mxmode, TY_CHAR)
	call salloc (pathx, mxpath, TY_CHAR)
	call salloc (graphtab, SZ_FNAME, TY_CHAR)
	call salloc (comptab, SZ_FNAME, TY_CHAR)

	# Convert strings from fortran to spp

	call f77upk (modstr, Memc[mode], mxmode)
	call f77upk (grftbl, Memc[graphtab], SZ_FNAME)
	call f77upk (cmptbl, Memc[comptab], SZ_FNAME)

	# Call spp interface procedure

	call getphotx (Memc[mode], Memc[graphtab], Memc[comptab], 
		       Memc[pathx], mxpath, phot)

	# Convert output string from spp to fortran

	call f77pak (Memc[pathx], path, mxpath)
	call sfree (sp)
end
