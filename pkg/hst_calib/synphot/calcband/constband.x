include	"../calcspec/mac.h"

# CONSTBAND -- Create a constant passband

int procedure constband( script, iw, nwv, sp)

char script[ARB]	# i: script containing commands and arguments 
int iw			# io: index of next word
int nwv			# i: number of wavelengths
real sp[ARB]		# o: constant passband

char cval[SIZE_STR]	# char value of constant

int ip, i, nchar
int ctor(), ctowrd()

real value

begin

	nchar = ctowrd(script, iw, cval, SIZE_STR)
	ip = 1
	nchar = ctor(cval, ip, value)

	do i=1,nwv
	   sp[i] = value 

end
