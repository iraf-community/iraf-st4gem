include	<error.h>

# SYNSORT -- Put synthetic photometry arrays in ascending order of wavelength

int procedure synsort1 (nwave, wave)

int	nwave		# i: number of wavelengths
real	wave[ARB]	# u: wavelength array
#--
int	dir, status
int	wavedir()

begin
	status = OK
	iferr {
	    dir = wavedir (nwave, wave)
	} then {
	    call synsortmsg ()
	    return (ERR)
	}

	# Return error if throughput is not monotonic, 
	# invert array if it is monotonic decreasing

	if (dir == 0) {
	    status = ERR

	} else if (dir == -1) {
	    call flip (nwave, wave)
	}

	return (status)
end

int procedure synsort2 (nwave, wave, thruput)

int	nwave		# i: number of wavelengths
real	wave[ARB]	# u: wavelength array
real	thruput[ARB]	# u: throughput array
#--
int	dir, status
int	wavedir()

begin
	status = OK
	iferr {
	    dir = wavedir (nwave, wave)
	} then {
	    call synsortmsg ()
	    return (ERR)
	}

	# Return error if throughput is not monotonic, 
	# invert array if it is monotonic decreasing

	if (dir == 0) {
	    status = ERR

	} else if (dir == -1) {
	    call flip (nwave, wave)
	    call flip (nwave, thruput)
	}

	return (status)
end

int procedure synsort3 (nwave, wave, thruput, thruerr)

int	nwave		# i: number of wavelengths
real	wave[ARB]	# u: wavelength array
real	thruput[ARB]	# u: throughput array
real	thruerr[ARB]	# u: throughput error array
#--
int	dir, status
int	wavedir()

begin
	status = OK
	iferr {
	    dir = wavedir (nwave, wave)
	} then {
	    call synsortmsg ()
	}

	# Return error if throughput is not monotonic, 
	# invert array if it is monotonic decreasing

	if (dir == 0) {
	    status = ERR

	} else if (dir == -1) {
	    call flip (nwave, wave)
	    call flip (nwave, thruput)
	    call flip (nwave, thruerr)
	}

	return (status)
end

int procedure synsort4 (nwave, wave, thruput, thruerr, waverr)

int	nwave		# i: number of wavelengths
real	wave[ARB]	# u: wavelength array
real	thruput[ARB]	# u: throughput array
real	thruerr[ARB]	# u: throughput error array
real	waverr[ARB]	# u: wavelength error array
#--
int	dir, status
int	wavedir()

begin
	status = OK
	iferr {
	    dir = wavedir (nwave, wave)
	} then {
	    call synsortmsg ()
	    return (ERR)
	}

	# Return error if throughput is not monotonic, 
	# invert array if it is monotonic decreasing

	if (dir == 0) {
	    status = ERR

	} else if (dir == -1) {
	    call flip (nwave, wave)
	    call flip (nwave, thruput)
	    call flip (nwave, thruerr)
	    call flip (nwave, waverr)
	}

	return (status)
end

procedure synsortmsg ()

#--
int	code
pointer	sp, errmsg

int	errget()

begin
	call smark (sp)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	code = errget (Memc[errmsg], SZ_LINE)
	call eprintf ("%s\n")
	call pargstr (Memc[errmsg])

	call sfree (sp)
end
