include	<fset.h>
include	<ctype.h>
include "function.h"

define	INC		10

# DONE_MESSAGE -- Print a percent done message

procedure done_message (text, iter, niter)

char	text		# i: text message printed as leader
int	iter		# i: current iteration
int	niter		# i: total number of iterations
#--
pointer	sp, msg
int	percent, oldper

begin
	# Allocate space for message

	call smark (sp)
	call salloc (msg, SZ_LINE, TY_CHAR)

	if (iter == 0) {
	    # Print a leading text string indicating what we are doing

	    call sprintf (Memc[msg], SZ_LINE, "%s. Percent done:")
	    call pargstr (text)

	    call print_message (Memc[msg])

	} else {
	    # Print the percentag value whenever a milestone is passed

	    percent = iter * 100 / niter
	    oldper = (iter - 1) * 100 / niter

	    if (percent / INC != oldper / INC) {
		call sprintf (Memc[msg], SZ_LINE, " %d")
		call pargi (percent)

		call print_message (Memc[msg])
	    }

	    if (iter == niter)
		call print_message ("\n")
	}

	call flush (STDERR)
	call sfree (sp)
end

# OBJ_MESSAGE -- Print a diagnostic message about the current object

procedure obj_message (ox, oy, flux, shp)

real	ox		# i: x position of object
real	oy		# i: y position of object
real	flux		# i: object flux
pointer	shp		# i: object shape descriptor
#--
pointer	sp, msg, name
string	format  "%s object at (%7.2f,%7.2f) with %0.3g counts\n"

begin
	# Allocate space for message

	call smark (sp)
	call salloc (msg, SZ_LINE, TY_CHAR)
	call salloc (name, SZ_FNAME, TY_CHAR)

	# Read shape name from shape descriptor

	if (shp == NULL) {
	    call strcpy ("point", Memc[name], SZ_FNAME)
	} else {
	    call strcpy (FUN_STR(shp,1), Memc[name], SZ_FNAME)
	}

	# Capitalize first letter of name

	if (IS_LOWER(Memc[name]))
	    Memc[name] = TO_UPPER(Memc[name])

	# Create message containing shape name, position, and flux

	call sprintf (Memc[msg], SZ_LINE, format)
	call pargstr (Memc[name])
	call pargr (ox)
	call pargr (oy)
	call pargr (flux)

	# Print the message

	call print_message (Memc[msg])

	call sfree (sp)
end

# ORDER_MESSAGE -- Print a diagnostic message for each spectral order

procedure order_message (order, nwave, wave)

int	order		# i: spectral order
int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: wavlength array
#--
pointer	sp, msg
int	iwave

string	format  "Spectral order %d. Bottom: %7.2f  Middle: %7.2f  Top: %7.2f\n"

begin
	# Allocate space for message

	call smark (sp)
	call salloc (msg, SZ_LINE, TY_CHAR)

	iwave = (nwave + 1) / 2

	call sprintf (Memc[msg], SZ_LINE, format)
	call pargi (order)
	call pargr (wave[1])
	call pargr (wave[iwave])
	call pargr (wave[nwave])

	call print_message (Memc[msg])
	call sfree (sp)
end

# PRINT_MESSAGE -- Print a diagnostic message to STDERR

procedure print_message (msg)

char	msg[ARB]	# i: message text
bool	verbose		# i: verbose message switch
#--
int	show		# print message if show is YES

data	show  / NO /

int	btoi()

begin
	if (show == NO)
	    return

	call putline (STDERR, msg)
	return

	# Set the message flag according to the value of verbose

	entry set_message (verbose)

	show = btoi (verbose)
	call fseti (STDERR, F_FLUSHNL, show)

	return
end
