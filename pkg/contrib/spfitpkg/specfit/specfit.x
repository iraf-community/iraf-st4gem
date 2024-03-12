###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	specfit
#
#  Description:	SPECFIT is an IRAF task for fitting complex continua with
#		mulitcomponents including emission lines.
#
#  Arguments:	See parameter file specfit.par
#
#  Returns:	Output is to user specified log files
#
#  Notes:	Information shared in common blocks defined in "specfit.com".
#		Most of the routines are written in SPP.
#
#  History:	May 1989	Gerard Kriss
#				Wrote major portions in SPP, based on SPLOT
#		Oct 1991	gak	Added v0, v1 for errors on 1D images
#		Feb 1992	hcf	Added option to compute expected
#					variance from model and use it
#					in computing Chisq.
#               June 1994       grimes  Added the ability to save current
#                                       condition of the model when
#                                       interrupted (^C)
#		June 1995	grimes  Added Knox Long's gridfit
#
###########################################################################

include	<imhdr.h>
include	<fset.h>
include <xwhen.h>
include "specfit.h"


procedure t_specfit ()

bool	plotdata, gridfit
char	typefit[SZ_FNAME]
int	nc, nfree, tfit
int	i, j, fd
int	logfiles	# list of log files

real	chisq
real	fpar[MAXFREE]	# Array of freely varying parameters

pointer	logfile, specfile, database, plotfile
pointer	ansrname, fitname, results, str, ranges

bool	clgetb()
int	clpopnu(), clgeti(), clgstr(), open(), clgfil(), clgwrd(), strlen()
real	clgetr()

extern xinterrupt  #interrupt function

int error_junk


pointer locpr()

include	"specfit.com"

begin

# Initialize global common variables
	call sfinit()

# Get task parameters.

	call malloc(specfile, SZ_FNAME, TY_CHAR)
	nc = clgstr("spectra", Memc[specfile], SZ_FNAME)

	debug = clgetb("debug")
	interact = clgetb("interact")
	err_from_model = clgetb("errors_from_model")
	plotdata = clgetb("plotdata")
	gridfit = clgetb("gridfit")
	tfit = clgwrd("type_of_fit",typefit,SZ_FNAME,"|marquadt|simplex|alternate|numrecipe|")

	logfiles = clpopnu("logfiles")

	call malloc(database, SZ_FNAME, TY_CHAR)
	call malloc(fitname, SZ_FNAME, TY_CHAR)
	call malloc(results, SZ_FNAME, TY_CHAR)
	call malloc(ansrname, SZ_FNAME, TY_CHAR)
	call malloc(plotfile, SZ_FNAME, TY_CHAR)
	nc = clgstr("database", Memc[database], SZ_FNAME)
	nc = clgstr("initial_fit", Memc[fitname], SZ_FNAME)
	nc = clgstr("final_fit", Memc[results], SZ_FNAME)
	nc = clgstr("flux_intervals", Memc[ansrname], SZ_FNAME)
	nc = clgstr("plot_file", Memc[plotfile], SZ_FNAME)

	itr = clgeti("max_iterations")
	tolerance = clgetr("fit_tolerance")
	v0 = clgetr("v0")
	v1 = clgetr("v1")

#Open log files and write a header.
	call fseti(STDOUT, F_FLUSHNL, YES)

	call malloc(logfile, SZ_FNAME, TY_CHAR)
	call malloc(str, SZ_FNAME, TY_CHAR)


	i = 0
	while ( clgfil(logfiles, Memc[logfile], SZ_FNAME) != EOF ) {
		fd = open(Memc[logfile], APPEND, TEXT_FILE)
		call sysid(Memc[str], SZ_FNAME)
		call fprintf(fd, "%s\n")
			call pargstr(Memc[str])
		call fprintf(fd, "  Task SPECFIT\n\n")

		i = i + 1
		logfd[i] = fd
	}
	nlogfd = i

	call mfree(logfile, TY_CHAR)
	call clpcls(logfiles)

# Start recording info in the log files
	for ( j = 1; j <= nlogfd; j = j + 1)
	{
		call fprintf(logfd[j],
"***************************************************************************\n")
	}

# Read in the spectrum
	call getspec(Memc[specfile])
#call printf("npts=%d\n")
#call pargi(npts)

# Add Object name and file name to the log
	for ( j = 1; j <= nlogfd; j = j + 1)
	{
		call fprintf(logfd[j], "** OBJECT: %s  FILE: %s **\n\n")
			call pargstr(specname)
			call pargstr(Memc[specfile])
		if ( err_from_model ) {  
		   call fprintf(logfd[j],"Using expected variance from model %g\n")
		   call pargb(err_from_model)
		}
		else {
		   call fprintf(logfd[j],"Using errors from input file\n")
	        }
	}

# Define sample ranges
	call malloc(fitsp, npts, TY_INT)
	call malloc(infit, npts, TY_INT)
	call malloc(ranges, 12 * MAXSAMPLE, TY_CHAR)
	nc = clgstr("Sample_ranges", Memc[ranges], 12 * MAXSAMPLE - 1)
	call decode_ranges( Memc[ranges] )
	call mfree(ranges, TY_CHAR)
	
# Get initial values for the fit components
	call sf_dbread(Memc[database], Memc[fitname], nfree, fpar)

# Get flux intervals for computing the results
	if ( strlen(Memc[ansrname]) > 0 ) {
		call sf_dbread2(Memc[database], Memc[ansrname])
	} else {
		ncont = -1
	}

#Set the error handler for interrupts to our own handler
        call xwhen(X_INT,locpr(xinterrupt),old_interrupt)


# If mode is interactive, enter the interactive graphics I/O sequence
	if ( interact ) {
		call sfinteract(nfree, fpar, chisq)
	} else {
# Not interactive, so do the fit and output the results
		if ( tfit == 1 && !gridfit) {
			call sfdomarquadt(nfree, fpar, chisq)
		} else if ( tfit == 2 && !gridfit) {
			call sfdosimplex(nfree, fpar, chisq)
		} else if ( tfit == 3 && !gridfit) {
			call sfalternate(nfree, fpar, chisq)
		} else if ( tfit == 4 && !gridfit) {
			call sfdonumrec(nfree, fpar, chisq)
		} else if (gridfit) {
			call sfdogrid(nfree,fpar,chisq,tfit)
		}
	}
	call sf_dbwrite(Memc[database], Memc[results])

#set error handler back to the built-in interrupts handler
        call xwhen(X_INT,old_interrupt,error_junk)


# Write final summary to the log files
	call sflogit(nfree, fpar, chisq)

# If output for plots requested, generate numbers for the components and
# output an ASCII file for use with smongo.
	if ( plotdata ) {
		call sfplot(Memc[plotfile])
	}

# Close the log files
	for ( i = 1; i <= nlogfd; i = i + 1)
	{
		call sysid(Memc[str], SZ_FNAME)
		call fprintf(logfd[i], "%s   --  SPECFIT is DONE.\n")
			call pargstr(Memc[str])
		call close (logfd[i])
	}

# Free allocated storage
	call mfree(str, TY_CHAR)
	call mfree(specfile, TY_CHAR)
	call mfree(plotfile, TY_CHAR)
	call mfree(database, TY_CHAR)
	call mfree(fitname, TY_CHAR)
	call mfree(results, TY_CHAR)
	call mfree(ansrname, TY_CHAR)
	call mfree(fitsp, TY_INT)
	call mfree(infit, TY_INT)
	call mfree(spectrum, TY_REAL)
	call mfree(lambdas, TY_REAL)
	call mfree(errors, TY_REAL)
	call mfree(cw, TY_REAL)
	call mfree(cf, TY_REAL)
	call mfree(pw, TY_REAL)
	call mfree(pf, TY_REAL)
	call mfree(aw, TY_REAL)
	call mfree(af, TY_REAL)


end

