# ophyd -- hydrogenic cross section
# bound-free plus free-free cross-section of a hydrogenic ion.
# formulae are from seaton (stimulated emission is omitted.!!!)
#
# input:
#	wave	r4 wavelength (angstroms)
#	temp	r4 kelvin temperature
#	z	r4 atomic number (1 = hi,2 = heii, etc.)
# output:
#	ophyd	r4 cross section per ion (cm^2/ion)
#--
# 1989 may kdh @ stsci .. modify comments
# jul 1989 Dave Bazell -- SPP version

real procedure ophyd( wave, temp, z )

real wave		# wavelenth (angstroms)
real temp		# temperature (Kelvin)
real z			# atomic number

real zz			# z*z
real beta		# ~1/kT
real x			# internal variable
real q			# internal variable
real qq			# q*q
real g1			# internal variable
real g2			# internal variable
real sum		# summed value
real gbf		# internal variable
real xg			# internal variable
real gff		# internal variable
real a			# internal variable
int n			# summation index 
int nn			# n*n
int n1			# lower summation limit
int n2			# upper summation limit 

real sig0
real k
data sig0/7.907e-18/
data k/3/

begin
# Check for bad initial values and return if there are any.
	if (wave <= 0. || temp <= 0. || z <= 0) return(0)
	zz = z*z
	beta = 157890.*zz/temp
	x = wave*zz/911.267

	q = x**(-0.3333333)
	qq = q*q
	g1 = 0.3458*q
	g2 = 0.03307*qq

# sum continua of first few levels explicitly
	n1 = sqrt(x) + 1
	n2 = n1 + k - 1
	sum = 0.
	do n = n1,n2
	{
	   nn = n*n
	   xg = x/nn - 0.5
	   gbf = 1. + g1*xg - g2*(xg*xg + 1.25)
	   sum = sum + exp(beta/nn)*gbf/(nn*n)
	}

# use continuum approximation for remaining levels
	xg = x/beta
	gff = 1. + g1*(xg + 0.5) - 2.*g2*(xg*(xg + 0.5) + 0.75)
	a = n2 + 0.5
	sum = sum + (exp(beta/(a*a)) - 1. + gff)/(beta + beta)

	return (sig0*x*x*x*sum/zz*exp(-beta))
end
