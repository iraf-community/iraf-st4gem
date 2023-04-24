include <chars.h>
include <gki.h>
include "psk.h"

# Define the default line width and change in line width.
define DEFAULT_LWSLOPE 10.

#---------------------------------------------------------------------------
.help psk_open 1May92 plot
.ih
NAME
psk_open -- Open the metacode file.  
.ih
USAGE
call psk_open (device, tty)
.ih
ARGUMENTS
.ls device (char[ARB])
Device name to open (NOT USED)
.le
.ls tty (pointer)
Pointer to graphcap descriptor
.le
.ih
DESCRIPTION
Parse the DD string from the graphcap entry for the device to get the 
file template and OS dispose command.  Generate a unique file name and
open the metacode file as a NEW_FILE. Save the dispose command for later.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psk_open (device, tty, out_file)

char    device[ARB]             # device name (NOT USED)
pointer tty                     # pointer to graphcap descriptor
char    out_file[ARB]           # Alternate output file.

# Declarations.
int i               # Generic.
int len_nodeprefix  # Length of the node prefix
int off             # Offset into the nodename/dispose commands
int op              # Pointer into the strings

char cap[2]         # Capability string from the graphcap.

pointer ddstr       # Parsed dispose string.
pointer devname     # Device name.
pointer fname       # Specific form of the spool file.
pointer ip          # Generic pointer.
pointer sp          # Stack pointer.
pointer spool       # Spool file name.
pointer raw_ddstr   # The graphcap version of the dispose command.
pointer tempfn      # Temporary file name.
pointer val         # Graphcap parameter value.

include "psk.com"

# Function prototypes.
int     open(), ttygets(), gstrcpy(), strlen()
bool    streq(), strne(), ttygetb ()

errchk  open, ttygets, ttygeti, ttygetb

begin
        
        # Allocate memory.
        call smark (sp)
        call salloc (raw_ddstr, SZ_DDSTR, TY_CHAR)
        call salloc (ddstr, SZ_DDSTR, TY_CHAR)
        call salloc (devname, SZ_FNAME, TY_CHAR)
        call salloc (spool, SZ_FNAME, TY_CHAR)
        call salloc (fname, SZ_PATHNAME, TY_CHAR)
        call salloc (tempfn, SZ_PATHNAME, TY_CHAR)
        call salloc (val, SZ_FNAME, TY_CHAR)
        
        # The DB flag may be set in the graphcap entry for an SGI device to
        # print debug messages during execution.
        ps_debug = ttygetb (tty, "DB")
        
        # The DD string is used to pass device dependent information to the
        # graphics device driver.
        if (ttygets (tty, "DD", Memc[raw_ddstr], SZ_DDSTR) <= 0)
            call error (1, "psikern: missing DD parameter in graphcap")
        
        # Expand any $(XX) graphcap parameter refferences in the DD string.
        op = ddstr
        for (ip=raw_ddstr;  Memc[ip] != EOS;  ip=ip+1)
            if (Memc[ip] == '$' && Memc[ip+1] == '(' && Memc[ip-1] != '\\') {
                
                # Graphcap parameter substitution.
                call strcpy (Memc[ip+2], cap, 2)
                if (ttygets (tty, cap, Memc[val], SZ_FNAME) <= 0) {
                    call eprintf ("Warning: graphcap field `%s' not found\n")
                    call pargstr (cap)
                } else {
                    for (off=val;  Memc[off] == '#';  off=off+1)
                        ;
                    for (;  Memc[off] != EOS;  off=off+1) {
                        Memc[op] = Memc[off]
                        op = op + 1
                    }
                }
                ip = ip + 4
                
            } else {
                
                # Ordinary character.
                Memc[op] = Memc[ip]
                op = op + 1
            }
        Memc[op] = EOS
        
        # Parse the DD string into the node/device name, temp file name,
        # and host dispose command.
        
        # Get node and device name (e.g., "node!device,...").
        len_nodeprefix = 0
        ip = ddstr
        for (op=devname;  Memc[ip] != EOS;  ip=ip+1)
            if (Memc[ip] == ',') {
                if (Memc[ip-1] == '\\') {
                    Memc[op-1] = ','
                    ip = ip - 1
                } else {
                    ip = ip + 1
                    break
                }
            } else {
                if (Memc[ip] == FNNODE_CHAR)
                    len_nodeprefix = op - devname + 1
                Memc[op] = Memc[ip]
                op = op + 1
            }
        Memc[op] = EOS
        
        # Get spoolfile root name.
        op = spool + gstrcpy (Memc[devname], Memc[spool], len_nodeprefix)
        for (;  Memc[ip] != EOS;  ip=ip+1)
            if (Memc[ip] == ',') {
                if (Memc[ip-1] == '\\') {
                    Memc[op-1] = ','
                    ip = ip - 1
                } else {
                    ip = ip + 1
                    break
                }
            } else {
                Memc[op] = Memc[ip]
                op = op + 1
            }
        Memc[op] = EOS
        
        # Get OS pathname of spoofile.
        if (streq (Memc[spool], "STDOUT")) {
            call strcpy (Memc[spool], ps_fname, SZ_PATHNAME)
        } else {
            call mktemp (Memc[spool], Memc[tempfn], SZ_PATHNAME)
            call fmapfn (Memc[tempfn], ps_fname, SZ_PATHNAME)
            call strupk (ps_fname, ps_fname, SZ_PATHNAME)
            
            # Get pathname of spoolfile on the remote node.  The call to
            # ki_fmapfn () is currently necessary to translate the filename for
            # the remote node, but may be replaced by the usual fmapfn () in a
            # future version of the kernel interface.
            call ki_fmapfn (Memc[tempfn], ps_fname, SZ_PATHNAME)
            call strupk (ps_fname, ps_fname, SZ_PATHNAME)
        }

        # Now that we have a filename, check to see if one was
        # specified in the argument list.  If so, then just use that
        # file.
        if (strlen (out_file) > 0)
            call strcpy (out_file, ps_fname, SZ_FNAME)
        
        if (ps_debug) {
            call eprintf ("psk_open: open device %s, outfile = %s\n")
            call pargstr (Memc[devname])
            call pargstr (ps_fname)
        }
        
        # Copy OS command for disposing of metacode file into common, replacing
        # all $F sequences in the command by the OS pathname of the spool file.
        op = gstrcpy (Memc[devname], ps_dispose, len_nodeprefix) + 1
        for (;  Memc[ip] != EOS;  ip=ip+1)
            if (Memc[ip] == '$' && Memc[ip-1] == '\\') {
                
                # Escape a $.
                ps_dispose[op-1] = '$'
                
            } else if (Memc[ip] == '$' && Memc[ip+1] == 'F') {
                
                # Filename substitution.
                for (i=1;  ps_fname[i] != EOS;  i=i+1) {
                    ps_dispose[op] = ps_fname[i]
                    op = op + 1
                }
                ip = ip + 1
                
            } else {
                
                # Ordinary character.
                ps_dispose[op] = Memc[ip]
                op = op + 1
            }
        
        ps_dispose[op] = EOS
        
        # Remove (delete) metacode file after issuing OS dispose command?
        ps_delete = ttygetb (tty, "RM")
        
        # Store each frame in a new file?
        ps_oneperfile = ttygetb (tty, "NF")
        
        # Should the output be buffered.  Generally yes if the output is a
        # file, else it should be no, such as when sending to a pipe.
        ps_buffer_output = ttygetb (tty, "BO")
        
        # Open a new metacode file.
        if (ps_oneperfile && strne(ps_fname, "STDOUT"))
            call psk_mkfname (ps_fname, ps_frame, Memc[fname], SZ_FNAME)
        else
            call strcpy (ps_fname, Memc[fname], SZ_FNAME)
        
        if (ps_debug) {
            call eprintf ("psk_open: open frame %2d, outfile = %s\n")
            call pargi (ps_frame)
            call pargstr (Memc[fname])
        }
        ps_fd = open (Memc[fname], NEW_FILE, TEXT_FILE)
        ps_buffer_eos = 1
        
        # Initialize the page counter.
        ps_frame = 0
        
        # Release memory.
        call sfree (sp)
        
end
#---------------------------------------------------------------------------
# End of psk_open
#---------------------------------------------------------------------------
