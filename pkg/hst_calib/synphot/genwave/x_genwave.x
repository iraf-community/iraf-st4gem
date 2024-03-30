include <tbset.h>

#* HISTORY *
#* E.Medeiros	01-Jan-89	original
#* B.Simon	02-Jun-94	rewritten for new genwave

# X_GENWAVE -- Generate a wavelength set and write to output table

procedure x_genwave()

#--
pointer	output		# Wavelength set table name
real 	minwave		# Minimum of wavelength range (Angstroms)
real	maxwave		# Maximum of wavelength range (Angstroms)
real	dwave		# Wavelength interval (Angstroms/pixel)
real	dvel		# Velocity interval (km/s/pixel)
pointer	wavecol		# Wavelength set table column name

int	nwave
pointer	sp, tp, cp, wave

string	blank    " "
string	flipped  "The maximum wavelength must be greater than the minimum"
string	novalues "You must set either dwave or dvelocity"
string	bothvals "You cannot set both dwave and dvelocity"

pointer	tbtopn()
real	clgetr()

begin
	# Allocate memory for strings

	call smark (sp)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (wavecol, SZ_FNAME, TY_CHAR)

	# Read and check task parameters

	call clgstr ("output", Memc[output], SZ_FNAME)
	minwave = clgetr ("minwave")
	maxwave = clgetr ("maxwave")

	dwave = clgetr ("dwave")
	dvel = clgetr ("dvelocity")
	call clgstr ("wavecol", Memc[wavecol], SZ_FNAME)

	if (minwave >= maxwave)
	    call error (1, flipped)

	if (IS_INDEFR(dwave) && IS_INDEFR(dvel)) 
	    call error (1, novalues)

	if (! IS_INDEFR(dwave) && ! IS_INDEFR(dvel))
	    call error (1, bothvals)

	# Calculate wavelength set

	call genwave (minwave, maxwave, dwave, dvel, nwave, wave)

	# Write wavelength set to output table

	tp = tbtopn (Memc[output], NEW_FILE, NULL) 
	call tbcdef (tp, cp, Memc[wavecol], "angstroms", 
		     blank, TY_REAL, 1, 1) 

	call tbtcre (tp)
                      
	call tbcptr (tp, cp, Memr[wave], 1, nwave)
	call tbtclo (tp)

	# release memory            

	call mfree (wave, TY_REAL)
	call sfree (sp)

end
