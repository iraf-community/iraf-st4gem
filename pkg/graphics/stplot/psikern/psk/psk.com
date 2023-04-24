.help psk.com 1May92 plot
.ih
NAME
psk.com -- The common for the PSK kernel.  
.ih
DESCRIPTION
The may keep states needed to decide what type of PostScript and at what
point in a PostScript file the output is in.  At this point it is not clear
who or what belongs here, but it will be kept around for future expansion.
.ih
HISTORY
.nf
  25Jul91 - Copied from sys$gio/sgikern/sgk.com and modified for psikern.
            Jonathan D. Eisenhamer, STScI.
.fe
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

int     ps_buffer_eos           # EOS position in the ps_buffer.
int     ps_fd                   # file descriptor of output file
int     ps_frame                # Number of frames plotted.
 
bool    ps_buffer_output        # TRUE to buffer the output.
bool    ps_debug                # print kernel debugging messages
bool    ps_delete               # delete metacode file after dispose
bool    ps_oneperfile           # store each frame in a new file

char    ps_buffer[MAX_CHAR]     # Output buffer.
char    ps_dispose[SZ_OSCMD]    # host dispose command
char    ps_fname[SZ_PATHNAME]   # metafile filename
char    ps_output[MAX_CHAR]     # Temporary output string.

common  /pskcom/ ps_buffer_eos, 
                 ps_fd, ps_frame, ps_buffer_output,
                 ps_debug, ps_delete, ps_oneperfile,
                 ps_buffer, ps_dispose, ps_fname, ps_output
#---------------------------------------------------------------------------
# End of psk.com
#---------------------------------------------------------------------------
