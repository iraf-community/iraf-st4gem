include	<pkg/gtools.h>
include	"icfit.h"

# IC_SHOW -- Show the values of the parameters.

procedure ic_show (ic, file, gt, nl)

pointer	ic			# ICFIT pointer
char	file[ARB]		# Output file
pointer	gt			# GTOOLS pointer
pointer	nl			# NLFIT pointer
#--

pointer	sp, str, coeff, error, flags
int	fd, npar, i

int	open(), nl_stati()
real	nl_statr()
long	clktime()
errchk	open()

begin
	fd = open (file, APPEND, TEXT_FILE)
	call smark (sp)
	call malloc (str, SZ_LINE, TY_CHAR)
	npar = nl_stati (nl, "npar")
	call salloc (coeff, npar, TY_REAL)
	call salloc (error, npar, TY_REAL)
	call salloc (flags, npar, TY_BOOL)
	call nl_gcoeff (nl, Memr[coeff], npar)
	call nl_gerrors (nl, Memr[error], npar)
	call nl_gflags (nl, Memb[flags], npar)

	call cnvtime (clktime(0), Memc[str], SZ_LINE)
	call fprintf (fd, "\n%s\n")
	    call pargstr (Memc[str])
	call gt_gets (gt, GTTITLE, Memc[str], SZ_LINE)
	call fprintf (fd, "%s\n")
	    call pargstr (Memc[str])
	call gt_gets (gt, GTYUNITS, Memc[str], SZ_LINE)
	if (Memc[str] != EOS) {
	    call fprintf (fd, "fit units = %s\n")
	        call pargstr (Memc[str])
	}
	call ic_gstr (ic, "function", Memc[str], SZ_LINE)
	call fprintf (fd, "function = %s\n")
	    call pargstr (Memc[str])

	do i = 1, npar {
	    call fprintf (fd, "(%s)  coeff%d = %g \t(%g)\n")
	    if (Memb[flags+i-1])	
	        call pargstr ("var")
	    else
	        call pargstr ("fix")
	    call pargi (i)
            call pargr (Memr[coeff+i-1])
            call pargr (Memr[error+i-1])
	}

	call fprintf (fd, "chi-sq        = %g\n")
	    call pargr (nl_statr (nl, "chisq"))
	call fprintf (fd, "rms           = %g\n")
	    call pargr (nl_statr (nl, "rms"))
	call fprintf (fd, "grow          = %g\n")
	    call pargr (IC_GROW(ic))
	call fprintf (fd, "naverage      = %d\n")
	    call pargi (IC_NAVERAGE(ic))
	call fprintf (fd, "low_reject    = %g\n")
	    call pargr (IC_LOW(ic))
	call fprintf (fd, "high_reject   = %g\n")
	    call pargr (IC_HIGH(ic))
	call fprintf (fd, "niterate      = %d\n")
	    call pargi (IC_NITERATE(ic))
	call fprintf (fd, "sample        = %s\n")
	    call pargstr (Memc[IC_SAMPLE(ic)])

	call sfree (sp)
	call mfree (str, TY_CHAR)
	call close (fd)
end
