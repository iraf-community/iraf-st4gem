.help psk_mkfname 1May92 plot
.ih
NAME
psk_mkfname -- Make the name of file N of a multiframe set.
.ih
USAGE
call psk_mkfname (root, num, outstr, maxch)
.ih
ARGUMENTS
.ls root (char[ARB])
The root filename.
.le
.ls num (int)
The file number.
.le
.ls outstr (char[maxch])
The constructed file name.
.le
.ls maxch (int)
The maximum length of the constructed filename.
.le
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psk_mkfname (root, num, outstr, maxch)

char    root[ARB]               # I:  root filename
int     num                     # I:  file number
char    outstr[maxch]           # O:  receives new filename
int     maxch                   # I:  Maximum length of outstr.

begin
        call sprintf (outstr, maxch, "%s.%d")
        call pargstr (root)
        call pargi (num)
end
#---------------------------------------------------------------------------
# End of psk_mkfname
#---------------------------------------------------------------------------
