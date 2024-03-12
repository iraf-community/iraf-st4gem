#include	<floatingpoint.h>
#include	<signal.h>

void igfpe_()
{
	sigfpe(FPE_FLTINEX_TRAP, SIGFPE_IGNORE);
	sigfpe(FPE_FLTDIV_TRAP, SIGFPE_IGNORE);
	sigfpe(FPE_FLTUND_TRAP, SIGFPE_IGNORE);
	sigfpe(FPE_FLTOVF_TRAP, SIGFPE_IGNORE);
	sigfpe(FPE_FLTOPERR_TRAP, SIGFPE_IGNORE);
	ieee_handler("set", "all", SIGFPE_IGNORE);
}
