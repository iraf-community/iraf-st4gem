###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	call fspec(x, val)
#
#  Description:	A procedure to evaluate multi-component spectra for SPECFIT
#
#  Arguments:	real	x	- abscissa at which to evaluate function
#
#  Returns:	real	val	- value of the function at x
#
#  Notes:	Information shared in common blocks defined in "specfit.com".
#
#  History:	May 1989	Gerard Kriss
#		August 15, 1989	Gerard Kriss	added Fe2 lines to empirical
#		August 17	gak		split out routine cspec()
#		March 11 1990	gak		added recomb, multigauss
#		March 12 1990	gak		added extcor
#		Jan. 6   1991	gak		added line profile to empirical
#		Feb. 19 1991	gak		added h2abs
#		Mar 3 1991	gak		No h2abs, multigauss, empirical
#						Added usercont,userline,userabs,
#							blackbody
#		April 12 1991	gak		Added modified Lorentz emission
#		April 16 1991	gak		Added damped absorption
#		April 17 1991	gak		Added multiple userabs option
#		April 18 1991	gak		Fixed omission of abs(dlambda)
#						in lorentz()
#		May 12 1991	gak		Added multiple usercont
#		Sep 30 1991	gak		modified usercont to interpolate
#						  in a grid of models
#		Oct  1 1991	gak		same for userabs and userline
#		Oct  2 1991	gak		added broken power law
# **********	Oct 21 1991	gak		Temporarily, labsorp additive
#		Apr 14 1992	gak		Added free-free: ffree()
#		Feb 17 1994	gak		Added tauabs()
#		Jun 24 1994	grimes		Added to rmodel so skip '#' in	
#						reading in user models
#		Jun 29 1994	grimes		Added ability to only plot
#						specific components
#		August 3 1994	grimes		Added extdrude
#		January 1995	grimes		Added ccmext
###########################################################################

include "specfit.h"

define	CW	Memr[cw+($1-1)+cpts*($2-1)]    # These define entries
define	CF	Memr[cf+($1-1)+cpts*($2-1)]    # into allocated storage
define	PW	Memr[pw+($1-1)+ppts*($2-1)]    # for user-supplied models
define	PF	Memr[pf+($1-1)+ppts*($2-1)]
define	AW	Memr[aw+($1-1)+apts*($2-1)]
define	AF	Memr[af+($1-1)+apts*($2-1)]

procedure fspec(x, val)
real	x, val

int	i, naa, ma
real	cval
real	mabsorp[MAXPAR/4]
real	aabsorp


include	"specfit.com"

begin
	val = 0.
	aabsorp = 0.
	naa = 0
	ma = 0

# Loop over all components.  Save absorption values for later.
	for ( i = 1; i <= ncomp; i = i + 1) {
		if ( inplot[i] == YES) {
			call cspec( i, x, cval)

			# Continuous absorption is always multiplicative; e.g.,
			# 9 == extinction, 7 == eabsorp.
			if ( comtyp[i] == 7 || comtyp[i] == 9 || comtyp[i] == 19 || comtyp[i] == 21) {

				ma = ma + 1
				mabsorp[ma] = cval
				# To allow for overlying emission, e.g. airglow, 
				# multiply absorption components as they come 
				# along.  components listed later will NOT be 
				# absorbed.
				val = val * cval

			} else if ( comtyp[i] == 6 || comtyp[i] == 18 ||
		     	comtyp[i] == 12 || comtyp[i] == 13 || comtyp[i] == 15 ) {

				ma = ma + 1
				mabsorp[ma] = cval
				# To allow for overlying emission, e.g. airglow, 
				# multiply absorption components as they 
				# come along.  components listed later will NOT 
				# be absorbed.
				val = val * cval

			# Temporarily assume that other absorption components simply
			# ADD rather than multiply
			#	aabsorp = aabsorp + (1. - cval)

			} else {
				val = val + cval
			}
		}
	}

# Now multiply in the absorption components
#	for ( i = 1; i <= ma; i = i + 1) {
#		val = val * mabsorp[i]
#	}
#	val = val * (1. - aabsorp)

end

procedure cspec( cmpn, x, val)
int	cmpn
real	x, val

real	c0, c1, norm, alpha, mean, sigma, skew, z, temp, ckey, alpha1, alpha2
real	linear(), powerlaw(), bb(), gaussian(), logarith(), labsorp(),
	logabs(), eabsorp(), recomb(), extcor(), usercont(), userline(),
	userabs(), lorentz(), dampabs(), bpl(), ffree(), tauabs(), extdrude(), 
	disk(), ccmext()


include	"specfit.com"

begin
		switch ( comtyp[cmpn] ) {
		case 1:
			c0 = par0[ parptr[1,cmpn] ]
			c1 = par0[ parptr[2,cmpn] ]
			val = linear(x, c0, c1)
		case 2:
			norm = par0[ parptr[1,cmpn] ]
			alpha = par0[ parptr[2,cmpn] ]
			val = powerlaw(x, norm, alpha)
		case 3:
			norm = par0[ parptr[1,cmpn] ]
			temp = par0[ parptr[2,cmpn] ]
			val = bb(x, norm, temp)
		case 4:
			norm = par0[ parptr[1,cmpn] ]
			mean = par0[ parptr[2,cmpn] ]
			sigma = par0[ parptr[3,cmpn] ]
			skew = par0[ parptr[4,cmpn] ]
			val = gaussian(x, norm, mean, sigma, skew)
		case 5:
			norm = par0[ parptr[1,cmpn] ]
			mean = par0[ parptr[2,cmpn] ]
			sigma = par0[ parptr[3,cmpn] ]
			skew = par0[ parptr[4,cmpn] ]
			val = logarith(x, norm, mean, sigma, skew)
		case 6:
			norm = par0[ parptr[1,cmpn] ]
			mean = par0[ parptr[2,cmpn] ]
			sigma = par0[ parptr[3,cmpn] ]
			val = labsorp(x, norm, mean, sigma)
		case 7:
			norm = par0[ parptr[1,cmpn] ]
			mean = par0[ parptr[2,cmpn] ]
			val = eabsorp(x, norm, mean)
		case 8:
			norm = par0[ parptr[1,cmpn] ]
			mean = par0[ parptr[2,cmpn] ]
			sigma = par0[ parptr[3,cmpn] ]
			skew = par0[ parptr[4,cmpn] ]
			val = recomb(x, norm, mean, sigma, skew)
		case 9:
			norm = par0[ parptr[1,cmpn] ]
			val = extcor(x, norm)
		case 10:
			norm = par0[ parptr[1,cmpn] ]
			mean = par0[ parptr[2,cmpn] ]
			z = par0[ parptr[3,cmpn] ]
			ckey = par0[ parptr[4,cmpn] ]
			val = usercont(x, norm, mean, z, ckey)
		case 11:
			norm = par0[ parptr[1,cmpn] ]
			mean = par0[ parptr[2,cmpn] ]
			z = par0[ parptr[3,cmpn] ]
			ckey = par0[ parptr[4,cmpn] ]
			val = userline(x, norm, mean, z, ckey)
		case 12:
			norm = par0[ parptr[1,cmpn] ]
			mean = par0[ parptr[2,cmpn] ]
			z = par0[ parptr[3,cmpn] ]
			ckey = par0[ parptr[4,cmpn] ]
			val = userabs(x, norm, mean, z, ckey)
		case 13:
			norm = par0[ parptr[1,cmpn] ]
			mean = par0[ parptr[2,cmpn] ]
 			sigma = par0[ parptr[3,cmpn] ]
			val = logabs(x, norm, mean, sigma)
		case 14:
			norm = par0[ parptr[1,cmpn] ]
			mean = par0[ parptr[2,cmpn] ]
			sigma = par0[ parptr[3,cmpn] ]
			skew = par0[ parptr[4,cmpn] ]
			val = lorentz(x, norm, mean, sigma, skew)
		case 15:
			norm = par0[ parptr[1,cmpn] ]
			mean = par0[ parptr[2,cmpn] ]
			sigma = par0[ parptr[3,cmpn] ]
			val = dampabs(x, norm, mean, sigma)
		case 16:
			norm = par0[ parptr[1,cmpn] ]
			mean = par0[ parptr[2,cmpn] ]
			alpha1 = par0[ parptr[3,cmpn] ]
			alpha2 = par0[ parptr[4,cmpn] ]
			val = bpl(x, norm, mean, alpha1, alpha2)
		case 17:
			norm = par0[ parptr[1,cmpn] ]
			mean = par0[ parptr[2,cmpn] ]
			val = ffree(x, norm, mean)
		case 18:
			norm = par0[ parptr[1,cmpn] ]
			mean = par0[ parptr[2,cmpn] ]
			sigma = par0[ parptr[3,cmpn] ]
			val = tauabs(x, norm, mean, sigma)
		case 19: # General extinction
			z = par0[parptr[1,cmpn] ] #ebv
			norm = par0[ parptr[2,cmpn] ] #gamma
			mean = par0[ parptr[3,cmpn] ] #w0
			sigma = par0[ parptr[4,cmpn] ] #C1
			alpha1 = par0[parptr[5,cmpn] ] #c2
			alpha2 = par0[parptr[6,cmpn] ] #c3
			skew = par0[parptr[7,cmpn] ] 	 #c4			
			val = extdrude(x,mean,sigma,alpha1,alpha2, skew,norm,z)

		case 20: #Disk
			z = par0[parptr[1,cmpn]]
			norm = par0[parptr[2,cmpn]]
			mean = par0[parptr[3,cmpn]]
			val = disk(x,z,norm,mean)

		case 21: #ccmext
			z = par0[parptr[1,cmpn]]
			norm = par0[parptr[2,cmpn]]
			val = ccmext(x,z,norm)

		default:
			Call error(1,"CSPEC: Unknown component type.")
		}
end

real procedure linear(x, c0, c1)
real	x, c0, c1

begin
	return (c0 + c1 * (x - 1000.))
end

real procedure powerlaw(x, norm, alpha)
real	x, norm, alpha

real	val

begin
	val = norm * (x / 1000.)**(-alpha)
	return (val)
end

real procedure bpl(x, norm, mean, alpha1, alpha2)
real	x, norm, mean, alpha1, alpha2

real	val

begin
	if ( x <= mean ) {
		val = norm * (x / mean)**(-alpha1)
	} else {
		val = norm * (x / mean)**(-alpha2)
	}
	return (val)
end

real procedure bb(x, norm, temp)
real	x, norm, temp

real	f0, val

begin
	# Black body in F lambda, normalized at 5500 A
	f0 = exp(26159.6/temp) - 1.
	val = norm * f0 * (5500./x)**5 / ( exp(1.438779e8/x/temp) - 1.)
	return (val)
end

real procedure gaussian(x, norm, mean, fwhm, skew)
real	x, norm, mean, fwhm, skew

real	sigma, delta, val, term
real	exp(), abs()

begin
	val = 0.0
	sigma = mean * ( fwhm / C ) / 2.354820044	# Convert from FWHM
	delta = abs( (x - mean) / sigma )
	if ( delta > 10.0 ) {
		val = 0.0
	} else {
		term = - delta**2 / 2.

		if ( skew == 1.0 ) {
			val = norm * exp( term ) / sigma / ROOT2PI
		} else {
			if ( x <= mean ) {
				term = exp ( term )
			} else {
				term = exp( term / skew**2 )
			}
			val = 2. * norm * term / sigma / ROOT2PI / (1. + skew) 
		}
	}
	return (val)
end

real procedure logarith(x, norm, mean, sigma, skew)
real	x, norm, mean, sigma, skew

real	val, alpha, alphar, fmax
real	alog()
begin
#
# This produces line profiles of the form f = fmax * (x/x0)**(-alpha)
	alpha = LN2 / alog(1. + sigma / C / 2.)
	if ( alpha <= 1.) alpha = 1.0001
	fmax = (alpha - 1.) * norm / mean / 2.
	if ( skew == 1.0 ) {
		if ( x >= mean )
			val = fmax * ( x / mean ) ** (-alpha)
		if ( x < mean )
			val = fmax * ( x / mean ) ** (alpha)
	} else {
		alphar = LN2 / alog(1. + skew * sigma / C / 2.)
		fmax = (alpha - 1.) *
			norm / mean / (1. + (alpha - 1.) / (alphar - 1.) )
		if ( x <= mean ) {
			val = fmax * ( x / mean ) ** (alpha)
		} else {
			val = fmax * ( x / mean ) ** (-alphar)
		}
	}
	return (val)
end

real procedure labsorp(x, ew, mean, fwhm)
real	x, ew, mean, fwhm

real	delta, val, term, norm, sigma
real	exp(), abs()

begin

# The absorption line is expressed as a Gaussian with FWHM in km/s and
# area in equivalent width in Angstroms.

	sigma = mean * ( fwhm / C ) / 2.354820044	# Convert from FWHM
	norm = ew / sigma / ROOT2PI
	delta = abs( (x - mean) / sigma )
	if ( delta > 10.0 ) {
		val = 1.0
	} else {
		term = -delta**2 / 2.
		val = 1. - norm * exp(term)
	}
	return (val)
end

real procedure tauabs(x, tau0, mean, fwhm)
real	x, tau0, mean, fwhm

real	delta, val, term, sigma
real	exp(), abs()

begin

# The absorption line is expressed as a Gaussian in optical depth
# with FWHM in km/s.

	sigma = mean * ( fwhm / C ) / 2.354820044	# Convert from FWHM
	delta = abs( (x - mean) / sigma )
	if ( delta > 10.0 ) {
		val = 1.0
	} else {
		term = -delta**2 / 2.
		val = exp( -tau0 * exp(term) )
	}
	return (val)
end

real procedure logabs(x, norm, mean, sigma)
real	x, norm, mean, sigma

real	val, alpha
real	alog(), exp()
begin
#
# This produces line profiles of the form f = fmax * (x/x0)**(-alpha)
# Ignore skew term for now, and produce symmetric profiles
	alpha = LN2 / alog(1. + sigma / C / 2.)
	if ( alpha <= 1.) alpha = 1.0001
	if ( x >= mean )
		val = norm * ( x / mean ) ** (-alpha)
	if ( x < mean )
		val = norm * ( x / mean ) ** (alpha)
	return ( exp( -val ) )
end

real procedure eabsorp(x, norm, mean)
real	x, norm, mean

real	val, tau
real	exp()

begin
	if ( x <= mean ) {
		tau = norm * (x / mean)**3
		val = exp( - tau )
		return (val)
	} else {
		return (1.)
	}
end

real procedure recomb(x, norm, lambda, temp, fwhm)
real	x, norm, lambda, temp, fwhm

real	delta, val, sigma
real	exp()

begin
	val = 0.0
	sigma = lambda * fwhm / C / 2.354820044
	delta = 1.440e8 * (1. / x - 1. / lambda) / temp
	if ( delta < 0.0 ) {
		val = norm * exp( -(x - lambda)**2 / sigma**2 / 2. )
	} else {
		# norm is assumed to be in flambda
		val = norm * (lambda / x)**2 * exp( - delta )
	}
	return (val)
end

real procedure usercont(x, norm, off, z, ckey)
real	x, norm, off, z, ckey

int	i, inp, cmodin1, cmodin2, ctemp, j
real	f1, f2, xx
real	val

int	open(), fscan()
real	interp()

pointer cwt, cft

include	"specfit.com"

begin

# Read in the continuum model file names
	if ( cfirst == 1 ) {
		inp = open("continuum.ls", READ_ONLY, TEXT_FILE)
		cnmod = 1
		for ( j = 1; j <= nlogfd; j = j + 1) {
		    call fprintf(logfd[j],
			"\nFiles read from continuum.ls for usercont:\n")
		}
		while ( fscan(inp) != EOF ) {
			call gargwrd(cmodname[1,cnmod],SZ_FNAME)
			call gargr(cmodkey[cnmod])

# List file names and keys to the log files
			for ( j = 1; j <= nlogfd; j = j + 1) {
			    call fprintf(logfd[j],"    %s %14.6g\n")
			    call pargstr(cmodname[1,cnmod])
			    call pargr(cmodkey[cnmod])
			}

			cnmod = cnmod + 1
		}
		for ( j = 1; j <= nlogfd; j = j + 1) {
			call fprintf(logfd[j],"\n")
		}
		cnmod = cnmod - 1
		call close(inp)
		cfirst = 0

		# Allocate storage space
		call malloc(cwt, PTSMODMAX*cnmod, TY_REAL)
		call malloc(cft, PTSMODMAX*cnmod, TY_REAL)
	
		ctemp = 0
		for ( i = 1; i <=cnmod; i = i + 1) {
			cpts = PTSMODMAX
			call rmodel(cmodname[1,i], cpts, Memr[cwt+PTSMODMAX*(i-1)],Memr[cft+PTSMODMAX*(i-1)])
			if (cpts > ctemp)
				{ctemp = cpts}
			}
		cpts = ctemp


		call malloc(cw, cpts*cnmod, TY_REAL)
		call malloc(cf, cpts*cnmod, TY_REAL)

		# Read in all the models
		for ( i = 1; i <= cnmod; i = i + 1 ) {
			for ( j = 1; j <= cpts; j = j + 1 ) {
				CW(j,i) = Memr[cwt+(j-1) + PTSMODMAX*(i-1)]
				CF(j,i) = Memr[cft+(j-1) + PTSMODMAX*(i-1)]
			}
		}

	call mfree(cwt, TY_REAL)
	call mfree(cft, TY_REAL)

	}

# Find appropriate files for interpolation, but do not extrapolate
	if ( ckey <= cmodkey[1] ) {		# Force use of model 1
		cmodin1 = 1
		cmodin2 = 1
	} else if ( ckey >= cmodkey[cnmod] ) {	# Force use of model cnmod
		cmodin1 = cnmod
		cmodin2 = cnmod
	} else {
		for ( i = 1; i < cnmod && ckey > cmodkey[i]; i = i + 1 ) {}
		i = i - 1
		cmodin1 = i
		cmodin2 = i + 1
		# Enforce limits
		if ( cmodin1 < 1 ) cmodin1 = 1
		if ( cmodin1 > cnmod ) cmodin1 = cnmod
		if ( cmodin2 < 1 ) cmodin2 = 1
		if ( cmodin2 > cnmod ) cmodin2 = cnmod
	}

# Interpolate the continuum intensity.
	xx = (x + off) / (1. + z)
	if ( cmodin1 == cmodin2 ) {
		val = interp(xx, CW(1,cmodin1), CF(1,cmodin1), cpts)
	} else {
		f1 = interp(xx, CW(1,cmodin1), CF(1,cmodin1), cpts)
		f2 = interp(xx, CW(1,cmodin2), CF(1,cmodin2), cpts)
		val = f1 + (ckey - cmodkey[cmodin1]) *
			(f2 - f1) / (cmodkey[cmodin2] - cmodkey[cmodin1])
	}
	return ( norm * val )
end

real procedure userline(x, norm, off, z, pkey)
real	x, norm, off, z, pkey

int	i, inp, pmodin1, pmodin2, ptemp, j
real	f1, f2, xx
real	val

int	open(), fscan()
real	interp()

pointer pwt,pft

include	"specfit.com"

begin

# Read in the line profile file names
	if ( pfirst == 1 ) {
		inp = open("profile.ls", READ_ONLY, TEXT_FILE)
		pnmod = 1
		for ( j = 1; j <= nlogfd; j = j + 1) {
		    call fprintf(logfd[j],
			"\nFiles read from profile.ls for userline:\n")
		}
		while ( fscan(inp) != EOF ) {
			call gargwrd(pmodname[1,pnmod],SZ_FNAME)
			call gargr(pmodkey[pnmod])

# List file names and keys to the log files
			for ( j = 1; j <= nlogfd; j = j + 1) {
			    call fprintf(logfd[j],"    %s %14.6g\n")
			    call pargstr(pmodname[1,pnmod])
			    call pargr(pmodkey[pnmod])
			}

			pnmod = pnmod + 1
		}
		for ( j = 1; j <= nlogfd; j = j + 1) {
			call fprintf(logfd[j],"\n")
		}
		pnmod = pnmod - 1
		call close(inp)
		pfirst = 0


                # Find the largest size of data (currently very ineffecient)
                call malloc(pwt, PTSMODMAX*pnmod, TY_REAL)
                call malloc(pft, PTSMODMAX*pnmod, TY_REAL)
 
                ptemp = 0
                for ( i = 1; i <= pnmod; i = i + 1 ) {
                        ppts = PTSMODMAX
                     	call rmodel(pmodname[1,i], ppts, Memr[pwt+PTSMODMAX*(i-1)],Memr[pft+PTSMODMAX*(i-1)])
           		if (ppts > ptemp)
                                { ptemp = ppts }
                }
                ppts = ptemp

                #Allocate storage space
                call malloc(pw, ppts*pnmod, TY_REAL)
                call malloc(pf, ppts*pnmod, TY_REAL)

                for ( i = 1; i <=pnmod; i = i + 1) {
                        for ( j = 1; j <=ppts; j = j + 1) {
                                PW(j,i) = Memr[pwt + (j - 1) + PTSMODMAX * (i - 1)]
                                PF(j,i) = Memr[pft + (j - 1) + PTSMODMAX * (i - 1)]
                        }
                }
 
                call mfree(pwt, TY_REAL)
                call mfree(pft, TY_REAL)



	}

# Find appropriate files for interpolation, but do not extrapolate
	if ( pkey <= pmodkey[1] ) {		# Force use of model 1
		pmodin1 = 1
		pmodin2 = 1
	} else if ( pkey >= pmodkey[pnmod] ) {	# Force use of model pnmod
		pmodin1 = pnmod
		pmodin2 = pnmod
	} else {
		for ( i = 1; i < pnmod && pkey > pmodkey[i]; i = i + 1 ) {}
		i = i - 1
		pmodin1 = i
		pmodin2 = i + 1
		# Enforce limits
		if ( pmodin1 < 1 ) pmodin1 = 1
		if ( pmodin1 > pnmod ) pmodin1 = pnmod
		if ( pmodin2 < 1 ) pmodin2 = 1
		if ( pmodin2 > pnmod ) pmodin2 = pnmod
	}

# Interpolate the line profile.
	xx = (x + off) / (1. + z)
	if ( pmodin1 == pmodin2 ) {
		val = interp(xx, PW(1,pmodin1), PF(1,pmodin1), ppts)
	} else {
		f1 = interp(xx, PW(1,pmodin1), PF(1,pmodin1), ppts)
		f2 = interp(xx, PW(1,pmodin2), PF(1,pmodin2), ppts)
		val = f1 + (pkey - pmodkey[pmodin1]) *
			(f2 - f1) / (pmodkey[pmodin2] - pmodkey[pmodin1])
	}
	return ( norm * val )
end

real procedure userabs(x, norm, off, z, akey)
real	x, norm, off, z, akey

int	i, inp, amodin1, amodin2, atemp, j
real	f1, f2, xx
real	val


int	open(), fscan()
real	interp()

pointer awt,aft

include	"specfit.com"

begin

# Read in the absorption model file names
	if ( afirst == 1 ) {
		inp = open("absorption.ls", READ_ONLY, TEXT_FILE)
		anmod = 1
		for ( j = 1; j <= nlogfd; j = j + 1) {
		    call fprintf(logfd[j],
			"\nFiles read from absorption.ls for userabs:\n")
		}
		while ( fscan(inp) != EOF ) {
			call gargwrd(amodname[1,anmod],SZ_FNAME)
			call gargr(amodkey[anmod])

# List file names and keys to the log files
			for ( j = 1; j <= nlogfd; j = j + 1) {
			    call fprintf(logfd[j],"    %s %14.6g\n")
			    call pargstr(amodname[1,anmod])
			    call pargr(amodkey[anmod])
			}

			anmod = anmod + 1
		}
		for ( j = 1; j <= nlogfd; j = j + 1) {
			call fprintf(logfd[j],"\n")
		}
		anmod = anmod - 1
		call close(inp)
		afirst = 0

		# Find the largest size of data (currently very ineffecient)
		call malloc(awt, PTSMODMAX*anmod, TY_REAL)
		call malloc(aft, PTSMODMAX*anmod, TY_REAL)

		atemp = 0
		for ( i = 1; i <= anmod; i = i + 1 ) {
			apts = PTSMODMAX
			call rmodel(amodname[1,i], apts, Memr[awt+PTSMODMAX*(i-1)],Memr[aft+PTSMODMAX*(i-1)])
			if (apts > atemp)
				{ atemp = apts }
		}
		apts = atemp			
		
		#Allocate storage space
		call malloc(aw, apts*anmod, TY_REAL)
		call malloc(af, apts*anmod, TY_REAL)

		for ( i = 1; i <=anmod; i = i + 1) {
			for ( j = 1; j <=apts; j = j + 1) {
				AW(j,i) = Memr[awt + (j - 1) + PTSMODMAX * (i - 1)]
				AF(j,i) = Memr[aft + (j - 1) + PTSMODMAX * (i - 1)] 
			}
		}

                call mfree(awt, TY_REAL)
                call mfree(aft, TY_REAL)



	}

# Find appropriate files for interpolation, but do not extrapolate
	if ( akey <= amodkey[1] ) {		# Force use of model 1
		amodin1 = 1
		amodin2 = 1
	} else if ( akey >= amodkey[anmod] ) {	# Force use of model anmod
		amodin1 = anmod
		amodin2 = anmod
	} else {
		for ( i = 1; i < anmod && akey > amodkey[i]; i = i + 1 ) {}
		i = i - 1
		amodin1 = i
		amodin2 = i + 1
		# Enforce limits
		if ( amodin1 < 1 ) amodin1 = 1
		if ( amodin1 > anmod ) amodin1 = anmod
		if ( amodin2 < 1 ) amodin2 = 1
		if ( amodin2 > anmod ) amodin2 = anmod
	}

# Interpolate the optical depth in the absorption model
	xx = (x + off) / (1. + z)
	if ( amodin1 == amodin2 ) {
		val = interp(xx, AW(1,amodin1), AF(1,amodin1), apts)
	} else {
		f1 = interp(xx, AW(1,amodin1), AF(1,amodin1), apts)
		f2 = interp(xx, AW(1,amodin2), AF(1,amodin2), apts)
		val = f1 + (akey - amodkey[amodin1]) *
			(f2 - f1) / (amodkey[amodin2] - amodkey[amodin1])
	}
#
# This is the normal return statment.
#
	return ( exp( -norm * val ) )
#
# The following modification scales the user absorption so that it
# 'partially covers' the source.  Used specifically for NGC 4151.
# norm=1 ==> full covering
#
#	return ( exp(-val) * norm + 1. - norm )

end

real procedure lorentz(x, norm, mean, fwhm, alpha)
real	x, norm, mean, fwhm, alpha

real	fwhma, denom, val
real	abs()

begin
	# Coded so that for alpha=2, norm is the total flux in the line
	val = 0.0
	fwhma = mean * fwhm / C		# Convert from velocity to wavelength
	denom = (abs( x - mean ) ) ** alpha + (fwhma / 2.)**2
	if ( denom < 1.e-15  ) {
		denom = 1.e-15
	}

	val = norm * fwhma / denom / TWOPI

	return (val)
end

real procedure dampabs(x, nhf, mean, gamma)
real	x, nhf, mean, gamma

real	nu, nu0
real	denom, val, tau
real	exp()

begin
	val = 0.0
	nu = C18 / x
	nu0 = C18 / mean
	denom = (nu - nu0) ** 2 + (gamma / FOURPI)**2
	if (denom < 1.e-15) denom = 1.e-15
	tau = PIE2MC * nhf * gamma / denom / FOURPI2
	val = exp( -tau )
	return (val)
end

real procedure ffree(x, norm, temp)
real	x, norm, temp

real	 val
real	exp()

begin
	val = 0.0
	# norm is assumed to be in flambda at 5500 A
	val = norm * (5500./x)**2 * exp(-1.438779e8/x/temp)
	return (val)
end

# Generic routine to read user-input model files
procedure rmodel(modname, npts, w, f)
char	modname[ARB]
int	npts
real	w[ARB], f[ARB]

char 	line[SZ_LINE], in
int	inp, i, j, junk
int	open(), fscan(), getc(), getline()

begin
# call fprintf(STDOUT,"Opening file %s\n")
# call pargstr(modname)
	inp = open(modname, READ_ONLY, TEXT_FILE)
	i = 1

	#Gets input ignores lines leading with a `#`
	while ( getc(inp,in) != EOF && i <= npts ) {
		if (in!='#') {
			call ungetc(inp,in)
			junk = fscan(inp)
				call gargr(w[i])
				call gargr(f[i])
			i = i + 1
		} else {
		junk = getline(inp,line)
		}
	}

	# If i still < npts, fill with zeroes out to npts
	for ( j = i; j <= npts; j = j + 1 ) {
		w[j] = 0.0
		f[j] = 0.0
	}
	# Set npts equal to the number of points actually read
	npts = i - 1
# call printf("%s  %d\n")
# call pargstr(modname)
# call pargi(npts)
	call close(inp)
end



real procedure interp(xx, x, y, npts)
real	xx
real	x[ARB], y[ARB]
int	npts

int	i
real	f

begin
	i = (xx - x[1]) / (x[2] - x[1]) + 1.0
	if ( i < 1 ) {
		f = 0.0
	}
	else if ( i > npts - 1 ) {
		f = 0.0
	} else {
	    if ( x[i+1]-x[i] != 0 ) {	
	        f = y[i] + (xx - x[i]) * (y[i+1]-y[i]) / (x[i+1]-x[i])
	    } else {
		f = 0.0
	    }
	}

	# if ( xx > 5234. && xx < 5235. ) {
	# 	call printf("%d %d  %f %f %e %e\n")
	# 	call pargi(i)
	# 	call pargi(npts)
	# 	call pargr(xx)
	# 	call pargr(x[i])
	# 	call pargr(y[i])
	# 	call pargr(f)
	# }
	return f
end



real procedure extdrude(x, gam, c1, c2, c3, c4, w0, ebv)
real	x,gam,c1,c2,c3,c4, w0, ebv

real 	ext, dx, xx, dru, fuv, av


begin


	av = 3.14

	xx = 10000. / x
	dru = xx**2 / ( gam * gam * xx * xx + (( xx - w0 )**2) )
	dx = xx - 5.9
	fuv = 0.5392 * (dx)**2 + 0.0564 * ( dx**3 )
	ext = av + c1 + c2 * xx + c3 * dru
	if ( xx > 5.9 )
		ext = c4 * fuv + ext 	

	ext = 10. ** ( -.4 * ext * ebv )

	return (ext)

end


# real procedure disk(x, A, Beta, Temp)
# real x, A, Beta, Temp
# real temp1, temp2, val

# real exp()

# begin

#  temp1 = A * ( (x / 20000. ) ** (-Beta ) )
#  temp2 = ( -PLANCK * CCM ) / ( K * Temp * x )

#  val = temp1 * exp ( temp2 )
#  return (val)
# end

real procedure disk(x, A, Beta, Lamb)
real x, A, Beta, Lamb
real temp1, val

real exp()

begin

  temp1 = A * (x / 20000.) ** (-Beta ) 

  val = temp1 * exp( - Lamb / x)

  return (val)
end

real procedure ccmext(w, ebv, rv)
real 	w, ebv, rv

real	x
real 	a, b, y, y2, y3, extcor, aext, xp


begin

x = 10000. / w
a = 0.0
b = 0.0

# Infrared Wavelengths
if ( x >= 0.3 && x <= 1.1 ) {
	xp = x**1.61
	a = 0.574 * xp
	b = -0.527 * xp
}


#Optical wavelength
if (x > 1.1 && x <= 3.3) {
	y = x - 1.82
	a = 1.0 + 0.17699 * y - 0.50477 * y**2 - 0.02427 * y**3 +
	          0.72085 * y**4 + 0.01979 * y**5 - 0.77530 * y**6 +
		  0.32999 * y**7
	b = 0.0 + 1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 -
		  5.38434 * y**4 - 0.62551 * y**5 + 5.30260 * y**6 -
		  2.09002 * y**7
}

#Near-UV wavelengths
if ( x > 3.3 && x <= 8.0 ) {
	a = 0.
	b = 0.
	if ( x >= 5.9 ) {
		y = x - 5.9
		y2 = y**2
		y3 = y2 * y
		a = -0.04473 * y2 - 0.009779 * y3
		b = 0.21300 * y2 + .120700 * y3
	}
	a = a + 1.752 - 0.316 * x - 0.104 / (0.341 + (x-4.67)**2)
	b = b - 3.090 + 1.825 * x + 1.206 / (0.263 + (x-4.62)**2)
}

# Far-UV wavelengths
if ( x > 8.0 && x <= 20. ) {
		y = x - 8.0
		y2 = y**2
		y3 = y * y2
		a = -1.073 - 0.628 * y + 0.137 * y2 - 0.070 * y3
		b = 13.670 + 4.257 * y - 0.420 * y2 + 0.374 * y3
}
	
aext = rv * a + b
extcor = 10.**(-0.4 * ebv * aext)

return (extcor)

end
