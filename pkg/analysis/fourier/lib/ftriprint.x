include <imhdr.h>		# for IM_HDRFILE
include "../fourier.h"

# ft_ri_print -- print names of real & imaginary parts
# This routine prints two lines.  Each line begins with the prefix string,
# then has " r:  " or " i:  " followed by the name of the real or imaginary
# part.  If the real or imaginary part does not exist, the string "[none]"
# will be printed instead of the name.
#
# ft_strmrg merges two versions of the image name, the user-supplied name
# FT_NAME_R (and FT_NAME_I) and the IMIO name IM_HDRFILE(im).  The former
# has already been modified to include the "r" (and "i") but may lack an
# extension, while the latter will not include any group number or image
# section, if those were specified.  FT_NAME_R will be the complete name
# if the user specified an extension; if not, ft_strmrg merges in the
# extension from IM_HDRFILE(im).

# Phil Hodge, ?  Subroutine created.
# Phil Hodge, 15-Jan-1995  Use FT_NAME_R or I instead of IM_HDRFILE.

procedure ft_ri_print (ft, prefix)

pointer ft		# i: pointer to ft struct
char	prefix[ARB]	# i: prefix string
#--
pointer sp
pointer imgname		# scratch for an image name

begin
	call smark (sp)
	call salloc (imgname, SZ_LINE, TY_CHAR)

	call printf ("%s r:  ")
	    call pargstr (prefix)
	if (FT_REAL(ft) == YES) {
	    call ft_strmrg (FT_NAME_R(ft), IM_HDRFILE(FT_REPT(ft)),
			Memc[imgname], SZ_LINE)
	    call printf ("%s\n")
		call pargstr (Memc[imgname])
	} else {
	    call printf ("[none]\n")
	}

	call printf ("%s i:  ")
	    call pargstr (prefix)
	if (FT_IMAG(ft) == YES) {
	    call ft_strmrg (FT_NAME_I(ft), IM_HDRFILE(FT_IMPT(ft)),
			Memc[imgname], SZ_LINE)
	    call printf ("%s\n")
		call pargstr (Memc[imgname])
	} else {
	    call printf ("[none]\n")
	}

	call sfree (sp)
end

# ft_strmrg -- merge image names
# If s1 lacks an extension, the extension is taken from s2.  s2 will
# include the extension but may lack a group number or image section.

procedure ft_strmrg (s1, s2, s3, maxch)

char	s1[ARB]		# i: image name, possibly without extension
char	s2[ARB]		# i: image name from IM_HDRFILE
char	s3[maxch]	# o: merged image name
int	maxch		# i: allocated size of s3
#--
int	n, len1, len2
int	strlen(), ft_strcmp()
bool	streq()

begin
	if (streq (s1, s2)) {
	    call strcpy (s1, s3, maxch)
	    return
	}

	len1 = strlen (s1)
	len2 = strlen (s2)

	# Find n such that s1[1:n] is equal to s2[1:n].  In general,
	# n could be zero, but it won't be in this situation.
	n = ft_strcmp (s1, s2)

	if (n == len1) {

	    # s1 is a substring of s2.
	    call strcpy (s2, s3, maxch)

	} else if (n == len2) {

	    # s2 is a substring of s1.
	    call strcpy (s1, s3, maxch)

	} else {

	    # Copy the name that includes the extension, then append the
	    # group number and/or image section.
	    call strcpy (s2, s3, maxch)
	    call strcat (s1[n+1], s3, maxch)
	}
end

# ft_strcmp -- compare strings
# This routine returns N such that s1 and s2 agree up to the Nth character.
# N will be zero if the strings do not have the same first character.

int procedure ft_strcmp (s1, s2)

char	s1[ARB], s2[ARB]	# i: strings to be compared
#--
int	i

begin
	i = 1

	while (s1[i] == s2[i] && s1[i] != EOS)		# implies s2[i] != EOS
	    i = i + 1

	return (i - 1)
end
