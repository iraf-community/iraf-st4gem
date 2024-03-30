# SPLITMODE -- splits a modestring with a dash (-) into two strings

# jun 1989 keith horne @ stsci
# jul 1989 kdh @ stsci - require ' - ' in color indices.
# Sep 1989 Dave Bazell - SPP version

procedure splitmode( mode, nmode, mode1, mode2 )

char	mode[ARB]	# i: Input mode string
int	nmode		# o: Number of modes output (1 or 2)
char	mode1[ARB]	# o: Mode string before the ' - ', if any
char	mode2[ARB]	# o: Mode string after the  ' - '

char	modsep
data	modsep /'-'/

int 	ic, jc, dash, nl
int	strlen(), stridx()

begin

	if( strlen(mode) == 0 ) return

	# split color indices into mode1 and mode2

	dash = stridx(modsep,mode)
	nl = stridx('\n', mode)
	if (nl > 0)
	   mode[nl] = EOS

	nmode = 1
	call strcpy(mode, mode1, SZ_LINE)
	if ( dash <= 0 )	# No dash, only 1 mode
	   return

	# Terminate mode1 with EOS if dash surrounded by spaces
	else if (mode1[dash-1] == ' ' && mode1[dash+1] == ' ' )
	   mode1[dash-1] = EOS

	# Copy second mode into mode2
	nmode = 2
	jc = 0
	for (ic = dash + 2;mode[ic] != EOS && mode[ic] != '\n';ic = ic + 1) {
	   jc = jc + 1
	   mode2[jc] = mode[ic]
	}
	mode2[jc+1] = EOS
end
