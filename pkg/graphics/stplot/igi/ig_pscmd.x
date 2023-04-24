include <gset.h>
include "igi.h"

# 21 May 1997 - WJH 
#
# Create the Postscript string which will define the 
#   PDFMARK command desired.  
#
# The input string for this command 'psmark' will be completely
#   defined by the user and passed along verbatim to the Postscript.
#
#  This function will only output to those devices which use the
#   taskname 'tn=psikern' in their graphcap entry.
#  
#  *********************** 
#  It will be the user's responsibility to pass error-free PDFMARK
#   commands through PSCMD.
#  ***********************
#
# This function is based on MGOSTR, but hard-wires
# the necessary settings for this particular command
# Since this is creating a POSTSCRIPT command, and not a 
# text string to be printed on the page, the following settings
# are fixed for this command:
#      XPOS,YPOS = 0.0,0.0
#      ANGLE = 0.0
#      SIZE = 1.0
#      JUSTIFICATION = 6 
# This command can be 
# given anywhere in the IGI script to place a PDFMARK on the
# plotting page.
  
procedure ig_pscmd (igs)

pointer	igs		# igi parameters structure

pointer	sp, line, psmark

begin
	call lcmdcat (igs, YES)

	call smark (sp)
	call salloc (line, SZ_LINE, TY_CHAR)
	call salloc (psmark, SZ_LINE, TY_CHAR)

	call igstarg (igs, Memc[line], SZ_LINE)

	if (DEBUG_OUTPUT(igs) == YES) {
	    call eprintf ("PostScript command:  %s \n")
		call pargstr (Memc[line])
	}
	
	call ii_pscmd (igs, Memc[line], Memc[psmark])

	call cmdcat (igs, YES)
	call sfree  (sp)
end


procedure ii_pscmd (igs, label, mark)

pointer	igs		# igi parameters structure
char	label[ARB]
char	mark[ARB]

char    gp_tn[SZ_LINE]
int     gp_out

pointer		gp
real        wl,wr,wb,wt
real        vl,vr,vb,vt

int ggets()
int strsearch()

begin
	call gseti (GIO_GP(igs), G_CLIP, NO)

	call setltype (igs, SOLID_LINE)

	# The following code comes from 'MGOSTR' with the 
	# settings assumed to be those shown in the call.		
	#call mgostr (igs, 0.0, 0.0, label, 1.0, 0.0, 6, 4)

	gp = GIO_GP(igs)

    call ggview (gp, vl,vr,vb,vt)
    call ggwind (gp, wl, wr, wb, wt)
            
	# Don't clip at viewport boundary;  allow labels outside axes
	call gseti (gp, G_CLIP, NO)
    
    # Only send this string along if the output device is Postscript
    # Search the graphcap entry for the output device for the 
    #  value of the taskname 'tn'.  This will be 'psikern' for Postscript
    #  devices which support this pass-through feature.
    gp_out = ggets (gp, "tn", gp_tn, SZ_LINE)

    if (strsearch (gp_tn, "psik") > 0) { 
	    # Send the string '<mark>' directly to Postscript
	    #  	The escape sequence will be picked up in 'psitx.x'
	    # 	and sent to 'ps_out' without being modified.
	    call sprintf(mark, SZ_LINE, "%s\%s%sfU%s")
		    call pargstr(TEXT_ESCAPE)
		    call pargstr(TEXT_ESCAPE)
		    call pargstr(TEXT_ESCAPE)
		    call pargstr(label)

	    # Use GIO (possibly hardware) characters
	    call giostr (gp, mark, 0.0, 0.0, 6, 1.0, 0.0)
    }	
	call gamove (GIO_GP(igs), 0.0, 0.0)
	call gflush (GIO_GP(igs))
end
