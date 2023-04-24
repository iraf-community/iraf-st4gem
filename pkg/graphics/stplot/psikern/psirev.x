#---------------------------------------------------------------------------
.help psi_reverse 1May92 plot
.ih
NAME
psi_reverse -- Reverse a segmented text string.
.ih
USAGE
call psi_reverse (intext, outtext)
.ih
ARGUMENTS
.ls intext (char[ARB])
The input text string to reverse.  See DESCRIPTION for format.
.le
.ls outtext (char[ARB])
The reversed string.
.le
.ih
DESCRIPTION
The input text string is a "segmented" string consisting of the following
format:

    font char char ... EOS font char char ... EOS EOS

Where the font is the GKI font index (GT_ROMAN etc) for the following
null-terminated string.  The whole string ends with two EOS.  There can be
any number of font changes.  This routine reverses the order of the
characters in each segment, and also reverses the order of the segments.
.ih
SEE ALSO
t_psikern, stx_segment
.endhelp
#---------------------------------------------------------------------------

procedure psi_reverse (intext, outtext)

char    intext[ARB]                  # I:  The input segmented string.
char    outtext[ARB]                 # O:  The output reversed string.

# Declarations
int     fontp, inp, outp             # Pointers into the strings.

begin
        
        # Look for the double end-of-string marker and start from there.
        inp = 1
        while (! (intext[inp] == EOS && intext[inp+1] == EOS))
            inp = inp + 1
        inp = inp - 1
        
        # Now reverse the strings.
        outp = 2
        fontp = 1
        while (inp > 1) {
            if (intext[inp] == EOS) {
                outtext[fontp] = intext[inp+1]
                outtext[outp-1] = EOS
                fontp = outp
            } else
                outtext[outp] = intext[inp]
            outp = outp + 1
            inp = inp - 1
        }
        outtext[fontp] = intext[1]
        outtext[outp] = EOS
        outtext[outp+1] = EOS
        
end
#---------------------------------------------------------------------------
# End of psi_reverse
#---------------------------------------------------------------------------
