#* HISTORY *
#* B.Simon	30-Mar-95	original

# RDNOISE -- Read and compile a noise expression

procedure rdnoise (obsmode, grftable, cmptable, exptime, nread, pcode, maxcode)

char	obsmode[ARB]	# i: observation mode
char	grftable[ARB]	# i: graph table name
char	cmptable[ARB]	# i: component lookup table name
real	exptime		# i: exposure time
int	nread		# i: number of readouts
int	pcode[ARB]	# o: pseudocode array
int	maxcode		# i: max code length
#--
pointer	sp, noise

string	sepstr   " & "
string	keyword  "NOISE"

begin
	# Allocate memory for temporary string

	call smark (sp)
	call salloc (noise, maxcode, TY_CHAR)

	# Read the noise parameter from the throughput  table

	call rdthruhdr (obsmode, grftable, cmptable, keyword, 
			sepstr, Memc[noise], maxcode)

	# Compile the noise expression into pseudocode

	call compnoise (Memc[noise], exptime, nread, pcode, maxcode)

	# Free memory

	call sfree (sp)
end
