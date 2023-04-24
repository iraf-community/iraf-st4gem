#---------------------------------------------------------------------------
.help psiparams.com 4Nov94 source
.ih
NAME
psiparams -- Store input parameters from the task level.
.endhelp
#---------------------------------------------------------------------------
int     area_color              # Default color for filled areas.
char    bold_font[SZ_LINE]      # PostScript font for GIO Bold 
char    greek_font[SZ_LINE]     # PostScript font for GIO Greek
char    italic_font[SZ_LINE]    # PostScript font for GIO Italic
int     line_color              # Default color for lines.
char	lut_gr[SZ_LINE]		# Graphics LUT.
char	lut_im[SZ_LINE]		# Image LUT.
int     marker_color            # Default color for markers.
char    out_file[SZ_LINE]       # Output file name.
int     proportional            # Use variable-spaced text?
char    roman_font[SZ_LINE]     # PostScript font for GIO ROMAN
int     text_color              # Default color for text.

common /psiparams/ 	area_color,
			bold_font,
			greek_font,
			italic_font,
			line_color,
			lut_gr,
			lut_im,
			marker_color,
			out_file,
			proportional,
			roman_font,
			text_color	
#---------------------------------------------------------------------------
# End of psiparams.com
#---------------------------------------------------------------------------
