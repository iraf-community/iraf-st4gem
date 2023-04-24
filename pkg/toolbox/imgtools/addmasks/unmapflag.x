# UNMAPFLAG -- Map a set of bit flags back onto their original values
#
# B.Simon	04-Apr-91	First Code
# B.Simon	16-Jul-93	Modified to only map one pixel

procedure unmapflag (flag, nflag, inpix, outpix)

int	flag[ARB]	# i: List of flag values
int	nflag		# i: Number of flag values
int	inpix		# i: Input pixel
int	outpix		# o: Output pixel
#--
int	iflag, bitval

begin

	if (nflag == 0) {
	    outpix = inpix

	} else {
	    bitval = 1
	    outpix = 0
	    
	    # Add the flag value whenever the corresponding bit is set
	    
	    do iflag = 1, nflag {
		if (inpix < bitval) {
		    break
		    
		} else if (and (inpix, bitval) != 0) {
		    outpix = outpix + flag[iflag]
		}
		
		bitval = 2 * bitval
	    }
	}

end
