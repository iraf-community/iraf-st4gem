define	LEN_PARSTRUCT	7

define	PAR_NMODE	Memi[$1]
define	PAR_NPHOT	Memi[$1+1]
define	PAR_MODELEN	Memi[$1+2]
define	PAR_PHOTLEN	Memi[$1+3]
define	PAR_MODEPTR	Memi[$1+4]
define	PAR_PHOTPTR	Memi[$1+5]
define	PAR_DATAPTR	Memi[$1+6]

define	PAR_MODE	Memc[PAR_MODEPTR($1)+PAR_MODELEN($1)*($2-1)]
define	PAR_PHOT	Memc[PAR_PHOTPTR($1)+PAR_PHOTLEN($1)*($2-1)]
define	PAR_DATA	Memr[PAR_DATAPTR($1)+PAR_NPHOT($1)*($2-1)+($3-1)]

#* HISTORY *
#* B.Simon	10-Oct-94	original

# OPN_PARMATRIX -- Create a matrix of the parameter results

pointer procedure opn_parmatrix (list, obslen, namelist, photvar)
				 
pointer	list		# i: list of observation modes
int	obslen		# i: maximum length of obsmode
char	namelist[ARB]	# i: list of possible photometry parameters
int	photvar[ARB]	# i: flags of parameters to save
#--
int	ic, iphot, nphot, photlen, imode, nmode, namelen
pointer	sp, name,par

int	numlist(), nxtlist(), word_fetch()

begin
	# Allocate memory for temporary string

	call smark (sp)
	call salloc (name, SZ_FNAME, TY_CHAR)

	# Get sizes of arrays to be allocated

	ic = 1
	iphot = 0
	nphot = 0
	photlen = 0
	nmode = numlist (list)

	repeat {
	    namelen = word_fetch (namelist, ic, Memc[name], SZ_FNAME)
	    if (namelen == 0)
		break

	    iphot = iphot + 1
	    if (photvar[iphot] == YES) {
		photlen = max (namelen, photlen)
		nphot = nphot + 1
	    }
	} 

	# Allocate data structure

	call malloc (par, LEN_PARSTRUCT, TY_STRUCT)
	call malloc (PAR_MODEPTR(par), nmode*(obslen+1), TY_CHAR)
	call malloc (PAR_PHOTPTR(par), nphot*(photlen+1), TY_CHAR)
	call malloc (PAR_DATAPTR(par), nmode*nphot, TY_REAL)

	PAR_NMODE(par) = nmode
	PAR_NPHOT(par) = nphot
	PAR_MODELEN(par) = obslen + 1
	PAR_PHOTLEN(par) = photlen + 1

	# Initialize date array to zero

	call amovkr (0.0, PAR_DATA(par,1,1), nmode*nphot)

	# Fill mode array with obsmode strings

	imode = 1
	while (nxtlist (list, Memc[name], SZ_FNAME) != EOF) {
	    call strcpy (Memc[name], PAR_MODE(par,imode), obslen)
	    imode = imode + 1
	}

	# Fill photname array with flagged photmetry parameter names

	ic = 1
	iphot = 1
	nphot = 1
	while (word_fetch (namelist, ic, Memc[name], SZ_FNAME) > 0) {
	    if (photvar[iphot] == YES) {
		call strupr (Memc[name])
		call strcpy (Memc[name], PAR_PHOT(par,nphot), photlen)
		nphot = nphot + 1
	    }
	    iphot = iphot + 1
	}

	call sfree (sp)
	return (par)
end

# UPD_PARMATRIX -- Update a value in the matrix of photometry parameters

procedure upd_parmatrix (par, imode, iphot, value)

pointer	par		# i: parameter marix descriptor
int	imode		# i: mode index
int	iphot		# i: photometry parameter index
real	value		# i: value to store
#--

begin
	# iphot must be the index to the i'th flagged parameter

	PAR_DATA(par,imode,iphot) = value

end

# WRT_PARMATRIX -- Write photometry parameters to STDOUT and free descriptor

procedure wrt_parmatrix (par)

pointer	par		# i: parameter matrix descriptor
#--
int	modelen, iphot, imode

begin
	modelen = max (PAR_MODELEN(par), 10)

	# Print photometry parameters in groups of three

	for (iphot = 1; iphot <=  PAR_NPHOT(par) - 2; iphot = iphot + 3) {
	    call printf ("%*s%15s%15s%15s\n")
	    call pargi (modelen)
	    call pargstr ("# OBSMODE")
	    call pargstr (PAR_PHOT(par,iphot))
	    call pargstr (PAR_PHOT(par,iphot+1))
	    call pargstr (PAR_PHOT(par,iphot+2))

	    do imode = 1, PAR_NMODE(par) {
		call printf ("%*s%15.5g%15.5g%15.5g\n")
		call pargi (modelen)
		call pargstr (PAR_MODE(par,imode))
		call pargr (PAR_DATA(par,imode,iphot))
		call pargr (PAR_DATA(par,imode,iphot+1))
		call pargr (PAR_DATA(par,imode,iphot+2))
	    }

	    call printf ("\n")
	}

	# Print remaining photometry parameters

	if (iphot ==  (PAR_NPHOT(par)-1)) {
	    call printf ("%*s%15s%15s\n")
	    call pargi (modelen)
	    call pargstr ("# OBSMODE")
	    call pargstr (PAR_PHOT(par,iphot))
	    call pargstr (PAR_PHOT(par,iphot+1))

	    do imode = 1, PAR_NMODE(par) {
		call printf ("%*s%15.5g%15.5g\n")
		call pargi (modelen)
		call pargstr (PAR_MODE(par,imode))
		call pargr (PAR_DATA(par,imode,iphot))
		call pargr (PAR_DATA(par,imode,iphot+1))
	    }

	    call printf ("\n")

	} else if (iphot == PAR_NPHOT(par)) {
	    call printf ("%*s%15s\n")
	    call pargi (modelen)
	    call pargstr ("# OBSMODE")
	    call pargstr (PAR_PHOT(par,iphot))

	    do imode = 1, PAR_NMODE(par) {
		call printf ("%*s%15.5g\n")
		call pargi (modelen)
		call pargstr (PAR_MODE(par,imode))
		call pargr (PAR_DATA(par,imode,iphot))
	    }

	    call printf ("\n")
	}

	# Realease memory held by descriptor

	call mfree (PAR_MODEPTR(par), TY_CHAR)
	call mfree (PAR_PHOTPTR(par), TY_CHAR)
	call mfree (PAR_DATAPTR(par), TY_REAL)
	call mfree (par, TY_STRUCT)
end
