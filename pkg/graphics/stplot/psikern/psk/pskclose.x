include "psk.h"

#---------------------------------------------------------------------------
.help psk_close 1May92 plot
.hi
NAME
psk_close -- Close the metacode spool file and dispose of it.
.ih
USAGE
call psk_close (do_dispose,font_list)
.ih
ARGUMENTS
.ls do_disponse (boolean)
TRUE to actually dispose of the current file.  False for do not dispose.
.le
.ls font_list (pointer)
The font list to write out.
.le
.ih
DESCRIPTION
Close the metacode spool file and dispose of it to a host system
metacode translation task.  Delete the spool file when the OS command
completes, unless it has already been deleted by the task run.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psk_close (do_dispose,font_list)

bool	do_dispose	# I:  Actually send file to printer?
pointer font_list	# IO: The font list.

# Declarations.
int i          # Generic.

pointer fname  # File name.
pointer sp     # Stack pointer.

include "psk.com"

# Function prototypes.
int     oscmd ()

errchk  psk_flush, close, oscmd

begin
        
        # Allocate memory.
        call smark (sp)
        call salloc (fname, SZ_FNAME, TY_CHAR)
        
        if (ps_debug)
            call eprintf ("psk_close: close device\n")
        
        if (ps_fd != NULL) {
            
            if (ps_debug)
                call eprintf ("psk_close: Sending frame and flush commands.\n")
            
            # Finish the frame(page) but include an extra restore since no
            # more output will appear in this file.
            call psk_frame
            call psk_flush
            call psk_out ("%%Trailer")
            call psk_flush
            call psk_out ("R")
            call psk_out ("R")
            
            # Write out the fonts.  Don't use the buffer commands because all
            # the fonts must appear on one line- regardless.
            call psk_flush
            call fprintf (ps_fd, "%%DocumentFonts:")
            while (font_list != NULL) {
                call psk_popfont (font_list, ps_output, MAX_CHAR)
                call fprintf (ps_fd, " %s")
                call pargstr (ps_output)
            }
            call fprintf (ps_fd, "\n")
            
            # Write out the number of pages and close the file.
            call sprintf (ps_output, MAX_CHAR, "%%%%Pages: %d")
            call pargi (ps_frame)
            call psk_out (ps_output)
            call psk_flush
            call psk_out ("%%EOF")
            call psk_flush
            call close (ps_fd)
            ps_fd = NULL
        }
        
        # Send the dispose command to the host system.
	if (do_dispose) {
	    if (ps_debug) {
		call eprintf ("psk_close: dispose: %s\n")
		call pargstr (ps_dispose)
	    }
	    
	    if (ps_dispose[1] != EOS)
		if (oscmd (ps_dispose, "", "", "") != OK)
		    call eprintf ("psk_close: Warning: PSK graphics output dispose error\n")
	}
	
        # Delete the metacode or raster file if so indicated in the graphcap
        # entry for the device.
        if (ps_delete) {
            if (ps_debug) {
                call eprintf ("psk_close: delete metafile %s\n")
                call pargstr (ps_fname)
            }
            if (ps_oneperfile) {
                do i = 1, ps_frame {
                    call psk_mkfname (ps_fname, i, Memc[fname], SZ_FNAME)
                    iferr (call delete (Memc[fname]))
                        ;
                }
            } else iferr (call delete (ps_fname))
                ;
        }
        
        # Release memory.
        call sfree (sp)
        
end
#---------------------------------------------------------------------------
# End of psk_close
#---------------------------------------------------------------------------
