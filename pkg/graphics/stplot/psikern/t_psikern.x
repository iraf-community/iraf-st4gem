include <error.h>
include <gki.h>
include <fset.h>
include <fio.h>

# Default color for all graphics modes.
define  DEF_COLOR       1

#---------------------------------------------------------------------------
.help t_psikern 9Nov94 source
.ih
NAME
psikern -- The PostScript GIO Kernel
.endhelp
#---------------------------------------------------------------------------

procedure t_psikern()

# Declarations.
int deb[LEN_GKIDD]  # ?
int debug           # YES if the gki code should be printed in text form.
int dev[LEN_GKIDD]  # ?
int fd              # File descriptor.
int gkiunits        # YES if debug should use GKI units instead of NDC.
int list            # Input list of metacode files.
int verbose         # YES if verbose mode should be used.

pointer devname     # Device name.
pointer fname       # File name
pointer gki         # ?
pointer sp          # Stack pointer

#int statbuf, statmax # status values
#int	extstat, nruns
#pointer	extbuf

include "psiparams.com"

# Function prototypes.
bool clgetb()
int clpopni(), clgeti(), clgfil(), open(), btoi(), gki_fetch_next_instruction()
#int fstati()
#long fstatl()

begin
        
        # Allocate memory.
        call smark (sp)
        call salloc (fname, SZ_FNAME, TY_CHAR)
        call salloc (devname, SZ_FNAME, TY_CHAR)
	
#	call eprintf("Allocated memory in PSIKERN...\n")

        # Open list of metafiles to be decoded.
        list = clpopni ("input")
        
        # Get parameters.
        call clgstr ("device", Memc[devname], SZ_FNAME)
        if (clgetb ("generic")) {
            roman_font[1] = EOS
            greek_font[1] = EOS
            bold_font[1] = EOS
            italic_font[1] = EOS
            out_file[1] = EOS
	    lut_gr[1] = EOS
	    lut_im[1] = EOS
            proportional = INDEFI
            area_color = INDEFI
            line_color = INDEFI
            marker_color = INDEFI
            text_color = INDEFI
            debug = NO
            verbose = NO
            gkiunits = NO
        } else {
            call clgstr ("output", out_file, SZ_LINE)
            call clgstr ("roman_font", roman_font, SZ_LINE)
            call clgstr ("greek_font", greek_font, SZ_LINE)
            call clgstr ("bold_font", bold_font, SZ_LINE)
            call clgstr ("italic_font", italic_font, SZ_LINE)
            proportional = btoi (clgetb ("proportional"))
	    call clgstr ("graphics_lut", lut_gr, SZ_LINE)
	    call clgstr ("image_lut", lut_im, SZ_LINE)
            line_color = clgeti ("linecolor")
            marker_color = clgeti ("markercolor")
            text_color = clgeti ("textcolor")
            area_color = clgeti ("areacolor")
            debug   = btoi (clgetb ("debug"))
            verbose = btoi (clgetb ("verbose"))
            gkiunits = btoi (clgetb ("gkiunits"))
        }
        
        # Open the graphics kernel.
        call psi_open (Memc[devname], dev)
        call gkp_install (deb, STDERR, verbose, gkiunits)
        
        # Process a list of metacode files, writing the decoded metacode
        # instructions on the standard output.
        while (clgfil (list, Memc[fname], SZ_FNAME) != EOF) {
            
            # Open input file.
            iferr (fd = open (Memc[fname], READ_ONLY, BINARY_FILE)) {
                call erract (EA_WARN)
                next
            }
	    # Check on the status of the file buffer
		#statbuf = fstati(fd, F_FILESIZE)
           	#statmax = fstati(fd, F_BLKSIZE)
 
		#iferr(call malloc(extbuf, statbuf, TY_SHORT) ){ 
		#	call eprintf("Not enough memory for extended buffer...\n")
		#	extstat = NO
		#}
		
		#if(extstat == YES) {
		#	call fseti(fd, F_BUFPTR, extbuf)
		#	call fseti(fd, F_BUFSIZE, statbuf)
		#	call fseti(fd, F_BLKSIZE, statbuf)
		#} else	
		#	call fseti(fd, F_ADVICE, SEQUENTIAL)

		
 		#statbuf = fstati(fd, F_BLKSIZE)
           	#statmax = fstati(fd, F_BUFSIZE)

		#debug = YES
           # Process the metacode instruction stream.
            while (gki_fetch_next_instruction (fd, gki) != EOF) {
	     if (debug == YES)
                    call gki_execute (Mems[gki], deb)
                call gki_execute (Mems[gki], dev)
            }
 
            call close (fd)
        }

        # Close up and unsuck memory.
        call gkp_close ()
        call psi_close ()
        call clpcls (list)
        call sfree (sp)
        
end
#---------------------------------------------------------------------------
# End of t_psikern
#---------------------------------------------------------------------------
