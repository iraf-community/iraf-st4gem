# MKFITSNAME -- Create fits file name from input file name

procedure mkfitsname (ifile, ofile, maxch)

char	ifile[ARB]	# i: input file name
char	ofile[ARB]	# o: output file name
int	maxch		# i: max length of ouput file name
#--
char	dot
int	nc, try

data	dot  /'.'/

bool	strne()
int	strldx()

begin
	try = 0

	repeat {
	    try = try + 1

	    # Copy the file name root

	    nc = strldx (dot, ifile)
	    if (nc == 0) {
		call strcpy (ifile, ofile, maxch)
	    } else {
		call strcpy (ifile, ofile, nc-1)
	    }

	    # Add suffix to avoid collision with input name

	    if (try == 2)
		call strcat ("_new", ofile, maxch)

	    # Add the fits extension

	    call strcat (".fits", ofile, maxch)

	} until (try == 2 || strne (ifile, ofile))
end

