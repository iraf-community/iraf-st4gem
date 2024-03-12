#---------------------------------------------------------------------------
.help t_fparse 31Mar92 develop
.ih
NAME
t_fparse -- Break the pathname into its components.
.ih
USAGE
task fparse = t_fparse
call t_fparse
.ih
DESCRIPTION
t_fparse seperates a filename into its component parts and returns them
in parameters.
.ih
BUGS
This routine uses the "illegal" call to imparse.

This routine ignores what is called a "ksection".  Not clear what this is.

Though the iraf system parsing routines don't handle wildcards, the capability
is here.  The assumption is that directory specifications can never have
wildcards and anything after the last '.' is an extension, wildcards included.
.endhelp
#---------------------------------------------------------------------------
procedure t_fparse()

# Declarations
int     cl_index                # Cluster index.
int     cl_size                 # Cluster size.

pointer cluster                 # Cluster.
pointer directory               # The directory part.
pointer extension               # Extension.
pointer input                   # The input.
pointer ksection                # Kernel section.
pointer root                    # Rootname.
pointer section                 # Section.
pointer sp                      # Stack pointer.

#string  wildcards       "*?"

# Function prototypes.
int     strlen()
bool    clgetb()

begin
        
        call smark(sp)
        call salloc (cluster, SZ_LINE, TY_CHAR)
        call salloc (directory, SZ_LINE, TY_CHAR)
        call salloc (extension, SZ_LINE, TY_CHAR)
        call salloc (input, SZ_LINE, TY_CHAR)
        call salloc (ksection, SZ_LINE, TY_CHAR)
        call salloc (root, SZ_LINE, TY_CHAR)
        call salloc (section, SZ_LINE, TY_CHAR)
        
        # Get the filename.
        call clgstr ("input", Memc[input], SZ_LINE)

        # Parse the name.
	call fparse (Memc[input], Memc[directory], SZ_LINE, Memc[root],
		     SZ_LINE, Memc[extension], SZ_LINE, cl_index, cl_size,
		     Memc[section], SZ_LINE, Memc[ksection], SZ_LINE)

        # Write the parameters back.
        call clpstr ("directory", Memc[directory])
        call clpstr ("root", Memc[root])
        call clpstr ("extension", Memc[extension])
        call clputi ("cl_index", cl_index)
        call clputi ("cl_size", cl_size)
        call clpstr ("section", Memc[section])
	call clpstr ("ksection", Memc[ksection])

        # If verbose, write to standard output.
        if (clgetb ("verbose")) {
            call printf ("%s %s %s %d %d %s %s\n")
            if (strlen (Memc[directory]) > 0)
                call pargstr (Memc[directory])
            else
                call pargstr ("\"\"")
            if (strlen (Memc[root]) > 0)
                call pargstr (Memc[root])
            else
                call pargstr ("\"\"")
            if (strlen (Memc[extension]) > 0)
                call pargstr (Memc[extension])
            else
                call pargstr ("\"\"")
            call pargi (cl_index)
            call pargi (cl_size)
            if (strlen (Memc[section]) > 0)
                call pargstr (Memc[section])
            else
                call pargstr ("\"\"")

            if (strlen (Memc[ksection]) > 0)
                call pargstr (Memc[ksection])
            else
                call pargstr ("\"\"")
        }

        # That's all folks.
        call sfree(sp)
end
#---------------------------------------------------------------------------
# End of t_fparse
#---------------------------------------------------------------------------
