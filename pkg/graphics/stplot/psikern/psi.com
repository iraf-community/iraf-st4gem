.help psicom 1May92 st4gem.stplot
.ih
NAME
psi common.  
.ih
DESCRIPTION
A common is necessary since there is no graphics descriptor in the
argument list of the kernel procedures.  The stdgraph data structures
are designed along the lines of FIO: a small common is used to hold the time
critical data elements, and an auxiliary dynamically allocated descriptor is
used for everything else.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

pointer g_kt                            # kernel transform graphics descriptor
pointer g_tty                           # graphcap descriptor
int     g_ndraw                         # Number of drawing operations performed
int     g_nframes                       # Which frame we're on now.
int     g_maxframes                     # Maximum frames per metafile.
int     g_maxgrcolor                    # Current maximum graphics color.
int	g_defmaxgrcolor			# Default maximum graphics color.
int     g_nseg                          # Number of line segments drawn.
bool    g_debug                         # TRUE for debugging output.
bool	g_bpage				# TRUE if beginning a new page.
char    g_device[SZ_GDEVICE]            # force output to named device
char    g_output[SZ_LINE]               # A genericly used string for
                                        # many of the routines.

common	/psicom/	g_kt, g_tty, g_ndraw, g_nframes, g_maxframes, 
	                g_maxgrcolor, g_nseg, g_debug, g_bpage,
			g_device, g_output, g_defmaxgrcolor
#---------------------------------------------------------------------------
# End of psi.com
#---------------------------------------------------------------------------

