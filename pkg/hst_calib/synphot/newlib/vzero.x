include	<ctype.h>

define	MAXLIST		50
define	DASH		'-'
define	STEP		'x'
define	COMMA		','

#* HISTORY *
#* B.Simon	28-Apr-94	original

# GETVZERO -- Get the current value of vzero from the list

procedure getvzero (value)

real	value		# o: current value of vzero
#--
include	"vzero.com"

string	undefval  "Value of variable is undefined"

begin
	value = vzcur
	if (IS_INDEFR(value))
	    call printerr_int (undefval, 0)

end

# NRVZERO -- Compute number of vzero values in a single range

int procedure nrvzero (irange)

int	irange		# i: index to range
#--
include	"vzero.com"

string	zerostep  "Step size cannot be zero"

int	count

begin
	if (vzero[3,irange] == 0.0)
	    call printerr_real (zerostep, vzero[3,irange])

	# The additional 1.0e-5 is to handle roundoff problems

	count = ((vzero[2,irange] - vzero[1,irange]) + 
		 vzero[3,irange]) / vzero[3,irange] + 1.0e-5

	return (count)
end

# NUMVZERO -- Get number of vzero values in list

int procedure numvzero ()

#--
include	"vzero.com"

int	irange, count
int 	nrvzero()

begin
	count = 0
	do irange = 1, vzlast
	    count = count + max (nrvzero (irange), 1)

	return (count)
end

# NXTVZERO -- Increment current value of vzero

int procedure nxtvzero (value)

real	value		# o: New value of vzero
#--
include	"vzero.com"

int	status
int	nrvzero()

begin
	# Reset to first value if out of range

	status = OK
	if (vzindex > vzlast) {
	    vzcount = 0
	    vzindex = 0
	}

	# If count is zero, move to next range
	# Otherwise, add step to current value

	if (vzcount <= 0) {
	    vzindex = vzindex + 1
	    if (vzindex > vzlast) {
		vzcur = INDEFR
		vzcount = 0
		status = EOF

	    } else {
		vzcur = vzero[1,vzindex]
		vzcount = nrvzero (vzindex)
	    }

	} else {
	    vzcur = vzcur + vzero[3,vzindex]
	}

	# Decrement count

	vzcount = vzcount - 1

	# Return current value and status

	value = vzcur
	return (status)

end

# RDVZERO - Read a list of ranges for vzero from a string

procedure rdvzero (list)

char	list[ARB]	# i: string containing list of vzero values
#--
include	"vzero.com"

bool	num
char	sep
int	ic, ifield, irange
real	value
 
string	badrange  "Illegal range syntax"
string	toolong   "Too many ranges in list"

int	ctor()

begin
	ic = 1
	num = true
	ifield = 1
	irange = 1

	call undefvzero

	# Get the first, last and step values for each range in
	# the list from the input string

	repeat {
	    if (num) {
		# Next expected token in list is a number
		# If not found, this is an error unless it is the
		# first field, since the list should not end on a
		# separator

		if (ctor (list, ic, value) == 0) {
		    if (ifield == 1 && irange == 1) {
			break
		    } else {
			call printerr_str (badrange, list)
		    }
		}

		vzero[ifield,irange] = value

	    } else {
		# Get separator and update range accordingly

		while (IS_WHITE(list[ic]))
		    ic = ic + 1

		sep = list[ic]
		ic = ic + 1

		switch (sep) {
		case DASH:
		    if (ifield == 1) {
			ifield = 2
		    } else {
			call printerr_str (badrange, list)
		    }

		case STEP:
		    if (ifield == 2) {
			ifield = 3
		    } else {
			call printerr_str (badrange, list)
		    }

		case COMMA, EOS:
		    # Set defaults if range is incomplete

		    if (ifield == 1) {
			vzero[2,irange] = vzero[1,irange]
			vzero[3,irange] = 1.0
		    } else if (ifield == 2) {
			vzero[3,irange] = 
			    sign (1.0, vzero[2,irange] - vzero[1,irange])
		    }

		    irange = irange + 1
		    ifield = 1

		    if (sep == EOS) {
			vzlast = irange - 1
			break

		    } else if (irange > MAXLIST) {
			call printerr_str (toolong, list)
		    }

		default:
		    call printerr_str (badrange, list)
		}
	    }

	    num = ! num
	}

end

# UNDEFVZERO -- Set the range list to undefined

procedure undefvzero ()

#--
include "vzero.com"

begin
	vzindex = 0
	vzcount = 0
	vzcur = INDEFR

	# Install a fake range to force nxtvzero() to loop once

	vzlast = 1
	vzero[1,1] = INDEFR
	vzero[2,1] = INDEFR
	vzero[3,1] = 1.0

end
