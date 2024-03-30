include	<gset.h>

# MARK2INT -- Changes a gio character plot mark to corresponding integer 
# value.  If an unknown mark is found, the previous mark value is returned.

int procedure mark2int( mark )

char	mark[ARB]		# i: GIO mark marcro

int	markval
int	strsearch()

save

begin
	
	call strlwr( mark )


	# Find the matching mark

	if ( strsearch( "gl_clear", mark ) > 0 )
	      markval = GL_CLEAR
	else if ( strsearch( "gl_solid", mark ) > 0 )
	      markval = GL_SOLID
	else if ( strsearch( "gl_dashed", mark ) > 0 )
	      markval = GL_DASHED
	else if ( strsearch( "gl_dotted", mark ) > 0 )
	      markval = GL_DOTTED
	else if ( strsearch( "gl_dotdash", mark ) > 0 )
	      markval = GL_DOTDASH
	else if ( strsearch( "gm_point", mark ) > 0 )
	      markval = GM_POINT 
	else if ( strsearch( "gm_fill", mark ) > 0 )
	      markval = GM_FILL
	else if ( strsearch( "gm_box", mark ) > 0 )
	      markval = GM_BOX
	else if ( strsearch( "gm_plus", mark ) > 0 )
	      markval = GM_PLUS
	else if ( strsearch( "gm_cross", mark ) > 0 )
	      markval = GM_CROSS
	else if ( strsearch( "gm_diamond", mark ) > 0 )
	      markval = GM_DIAMOND
	else if ( strsearch( "gm_hline", mark ) > 0 )
	      markval = GM_HLINE
	else if ( strsearch( "gm_vline", mark ) > 0 )
	      markval = GM_VLINE
	else if ( strsearch( "gm_hebar", mark ) > 0 )
	      markval = GM_HEBAR
	else if ( strsearch( "gm_vebar", mark ) > 0 )
	      markval = GM_VEBAR
	else if ( strsearch( "gm_circle", mark ) > 0 )
	      markval = GM_CIRCLE

	return ( markval )
end
