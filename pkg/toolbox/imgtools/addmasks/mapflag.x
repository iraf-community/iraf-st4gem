# MAPFLAG -- Map a bad pixel value onto a set of binary flags
#
# B.Simon	04-Apr-91	First Code
# B.Simon	16-Jul-93	Modified to only map one pixel

int procedure mapflag (flag, nflag, inpix, outpix)

int	flag[ARB]	# i: List of flag values
int	nflag		# i: Number of flag values
int	inpix		# i: Input pixel
int	outpix		# o: Output pixel
#--
int	iflag, temp, bitval

begin
	if (nflag == 0) {
	    outpix = inpix
	    temp = 0

	} else {
	    outpix = 0
	    temp = inpix
	    bitval = 2 ** (nflag - 1)
	    
	    # Break each pixel value into its component flags
	    # Store each flag as a single bit
	    
	    do iflag = nflag, 1, -1 {
		if (temp == 0) {
		    break
		    
		} else if (temp >= flag[iflag]) {
		    temp = temp - flag[iflag]
		    outpix = or (outpix, bitval)
		}
		
		bitval = bitval / 2
	    }
	}

	if (temp != 0) {
	    return (ERR)
	} else {
	    return (OK)
	}

end
