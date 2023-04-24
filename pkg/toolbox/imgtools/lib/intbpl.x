include	<plio.h>
include	<plset.h>

# INTBPL -- interpret a string command and return the corresponding rop
#
# D. Giaretta, 28-Jun-1988	STScI	Original code

int procedure intbpl( comm, default, defop, rop)

char	comm[ARB]	# i: input command
int	default		# i: default rop
char	defop		# i: default op to go with default
int	rop		# o: rop created
#--


int	ncmd
int	strdic(), ctoi()
int	ip, nowhite(), strsearch(), stridxs(), strlen(), val, pval, valop
char	sval[SZ_LINE], cmdstr[SZ_LINE], cmdstr2[SZ_LINE]
string  cmds	"oPIX_CLRoPIX_SEToPIX_SRCoPIX_DSToPIX_NOT(PIX_SRC)o\
PIX_NOT(PIX_DST)oPIX_SRC&PIX_DSToPIX_SRC|PIX_DSToPIX_SRC^PIX_DSTo\
PIX_SRC&PIX_NOT(PIX_DST)oPIX_SRC|PIX_NOT(PIX_DST)o\
PIX_NOT(PIX_SRC)&PIX_DSToPIX_NOT(PIX_SRC)|PIX_DSTo\
PIX_NOT(PIX_SRC&PIX_DST)oPIX_NOT(PIX_SRC|PIX_DST)oPIX_NOT(PIX_SRC^PIX_DST)"
string	pixval	"PIX_VALUE"

begin

	rop = 0

	if ( nowhite( comm, comm, SZ_LINE) == 0)
	    return ( NULL)

	ip = 1
	if ( ctoi( comm, ip, val) != 0 ) {
	    rop = val
	    return ( 1)
	}

	call strupr( comm)

	pval = strsearch( comm, pixval) - strlen(pixval)
	if ( pval > 0 ) { 
	    if ( pval > 1 ) {
	        valop = comm[pval-1]
	    	call strcpy( comm, cmdstr, pval-2 )
	    } else {
		valop = 'B'
		call strcpy( EOS, cmdstr, SZ_LINE)
	    }
	    call strcpy( comm[pval], sval, SZ_LINE)
	} else {
	    valop = EOS
	    call strcpy( comm, cmdstr, strlen(comm) )
	    call strcpy( EOS, sval, SZ_LINE)
	}

	ncmd = strdic (cmdstr, cmdstr2, SZ_LINE, cmds) + 1
	switch (ncmd) {
	case 1:
	    if ( cmdstr[1] != EOS )
		return ( NULL )
	case 2:
	    rop = PIX_CLR
	case 3:
	    rop = PIX_SET
	case 4:
	    rop = PIX_SRC
	case 5:
	    rop = PIX_DST
	case 6:
	    rop = PIX_NOTSRC
	case 7:
	    rop = PIX_NOTDST
	case 8:
	    rop = PIX_SRC_AND_DST
	case 9:
	    rop = PIX_SRC_OR_DST
	case 10:
	    rop = PIX_SRC_XOR_DST
	case 11:
	    rop = PIX_SRC_AND_NOTDST
	case 12:
	    rop = PIX_SRC_OR_NOTDST
	case 13:
	    rop = PIX_NOTSRC_AND_DST
	case 14:
	    rop = PIX_NOTSRC_OR_DST
	case 15:
	    rop = PIX_NOT_SRC_AND_DST
	case 16:
	    rop = PIX_NOT_SRC_OR_DST
	case 17:
	    rop = PIX_NOT_SRC_XOR_DST
	default:
	    return ( NULL )
	}

	if ( sval[1] == EOS )
	    return ( 1)
	else {
	    ip = stridxs ("(", sval) + 1
	    if( ctoi( sval, ip, val) < 1)
	        return ( NULL )
	}

	switch (valop) {
	case 'B':
	    if ( defop == '+' )
	    	rop = default + PIX_VALUE(val)
	    else if ( defop == '|' )
	    	rop = or(default, PIX_VALUE(val) )
	    else if ( defop == '&' ) 
		# Note that here we must ensure the op part gets through
	    	rop = and( default, PIX_SET+PIX_VALUE(val) )
	    else if ( defop == '^' )
	    	rop = xor(default, PIX_VALUE(val) )
	    else
		rop = PIX_SET + PIX_VALUE(val)
	case '+':
	    rop = rop + PIX_VALUE(val)
	case '|':
	    rop = or ( rop, PIX_VALUE(val) )
	case '&':
		# Note that here we must ensure the op part gets through
	    rop = and( rop, PIX_SET+PIX_VALUE(val) )
	case '^':
	    rop = xor( rop, PIX_VALUE(val) )
	default:
	    return ( NULL )
	}

	return ( 1)

end
