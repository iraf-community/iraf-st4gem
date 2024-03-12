######################################################################
#	Synopsis:	procedure interrupt(code,nxt_handler)
#
#	Description:	INTERRUPT This intercepts the interrupt error
#			caused when the user types ^C.  It will then 
#			output the current status of the model to a
#			file so that any calculations that have already
#			been made will not be lost.
#			
#
#	Arguments:	int code		-identifies the error
#			int nxt_handler		-what will be the next
#						 procedure to be called
#	Returns:	none
#	
#	Notes
#
#	History:	June 1994	J Grimes	created
#
######################################################################


include "specfit.h"
include <xwhen.h>
include <config.h>

procedure xinterrupt(code,nxt_handler)

int code
int nxt_handler
int  nfree,nc
int fpar[MAXFREE]
pointer results, database

int clgstr()

include "specfit.com"

begin
	
# Insures that when this function is exited the normal error procedure will be
# used and cause the program to be interrupted 
	nxt_handler=old_interrupt
	call setpar(nfree,fpar)
	call update(nfree,fpar)
	
#Get the file names then write the data to them
	call malloc(database, SZ_FNAME, TY_CHAR)
	call malloc(results, SZ_FNAME, TY_CHAR)
	nc=clgstr("final_fit",Memc[results],SZ_FNAME)
	nc=clgstr("database",Memc[database],SZ_FNAME)

	call eprintf("\nCurrent status of the model being appended \n\tto file sf%s before program termination\n")
	call pargstr(Memc[results])

	call sf_dbwrite(Memc[database],Memc[results])	
        call mfree(results,TY_CHAR)
        call mfree(database,TY_CHAR)

#flush the error from the stream (helps prevent recursion)
	call flush(STDERR)
end




