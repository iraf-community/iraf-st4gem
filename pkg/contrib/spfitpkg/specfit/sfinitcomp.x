include "specfit.h"

procedure sfinitcomp()


include "specfit.com"

begin

	call strcpy("linear", compkeys[1,1], SZ_CNAME)
	call strcpy("powerlaw", compkeys[1,2], SZ_CNAME)
	call strcpy("blackbody", compkeys[1,3], SZ_CNAME)
	call strcpy("gaussian", compkeys[1,4], SZ_CNAME)
	call strcpy("logarith", compkeys[1,5], SZ_CNAME)
	call strcpy("labsorp", compkeys[1,6], SZ_CNAME)
	call strcpy("eabsorp", compkeys[1,7], SZ_CNAME)
	call strcpy("recomb", compkeys[1,8], SZ_CNAME)
	call strcpy("extcor", compkeys[1,9], SZ_CNAME)
	call strcpy("usercont", compkeys[1,10], SZ_CNAME)
	call strcpy("userline", compkeys[1,11], SZ_CNAME)
	call strcpy("userabs", compkeys[1,12], SZ_CNAME)
	call strcpy("logabs", compkeys[1,13], SZ_CNAME)
	call strcpy("lorentz", compkeys[1,14], SZ_CNAME)
	call strcpy("dampabs", compkeys[1,15], SZ_CNAME)
	call strcpy("bpl", compkeys[1,16], SZ_CNAME)
	call strcpy("ffree", compkeys[1,17], SZ_CNAME)
	call strcpy("tauabs", compkeys[1,18], SZ_CNAME)
	call strcpy("extdrude",compkeys[1,19],SZ_CNAME)
	call strcpy("disk",compkeys[1,20],SZ_CNAME)
	call strcpy("ccmext",compkeys[1,21],SZ_CNAME)
	#Add new component name here
#	call strcpy("NewName",compkeys[1,21],SZ_CNAME)


	ncpar[1] = 2
	ncpar[2] = 2
	ncpar[3] = 2
	ncpar[4] = 4
	ncpar[5] = 4
	ncpar[6] = 3
	ncpar[7] = 2
	ncpar[8] = 4
	ncpar[9] = 1
	ncpar[10] = 4
	ncpar[11] = 4
	ncpar[12] = 4
	ncpar[13] = 3
	ncpar[14] = 4
	ncpar[15] = 3
	ncpar[16] = 4
	ncpar[17] = 2
	ncpar[18] = 3
	ncpar[19] = 7
	ncpar[20] = 3
	ncpar[21] = 2
	#Add new number of parameters here
#	ncpar[21] = ?


end
