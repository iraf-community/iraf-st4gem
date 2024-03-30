# THERMAL -- setup for loop over wavelengths in THERMBACK  task

#* HISTORY *
#* A.Sivaramakrishnan      dd-Mmm-yy       original
#* V. Laidler              18-Nov-04	Add support for output units; clean up final conversion.

#include "libsynphot.h"
define  MAXLIST     100
define  MAXPARAM    3
define  SZ_KEYWRD       32              # Length of keyword
define  GRF_COMPID      1               # index to optical component col in graph tab
define  GRF_KEYWRD      2               # index to keyword col in graph tab
define  GRF_THMLID      5               # index to thermal component col in grap
# above defines are also in libsynphot.h under st4gem$lib/synphot
# but convention precludes including it 



procedure pflist(flist, n)
char flist[SZ_FNAME, ARB]
int n
int i
begin
	do i = 1, n {
		call  printf("\t\t%d\t%s\n")
		call pargi(i)
		call pargstr(flist[1, i])
		call flush(STDOUT)
	}
end



procedure thermal(omode, tmode, ogrftb, ocmptb, tcmptb, detcat,
                  form, wave, nwave, output, trate, verbose)

char   omode[ARB]   # i incoming  observation mode
char   tmode[ARB]   # i incoming  'special requests on temperatures' mode
char  ogrftb[ARB]   # i graph table for omode (optical synphot style)
char  ocmptb[ARB]   # i component lookup table file name
char  tcmptb[ARB]   # i component lookup table file nam
char  detcat[ARB]   # i table containing detector pixel scales
char    form[ARB]   # i form/units for output
real    wave[ARB]   # i array of wavelengths for the calculation
int  nwave          # i number of wavelengths in wave array
char  output[ARB]   # i name of output table
	
real trate          # o detected thermal photons per second per steradian
                    # in the specified instrument mode
bool verbose        # i verbosity control flag

#--
real   asumr()    # function call


# local variables:

#string  form     "counts" # detected photons,  
                          # HgCdTE (eg NICMOS) gain 1 between e and photons
                          # CCDs: gain adjustable by read electronics
                          # other detector physics... ???
#string  filecol  "FILENAME"
string  obscol   "OBSMODE"
#string  thlcol   "THLMODE"
int	    oncomp   # number of optical components found in lookup table
int	    tncomp   # number of thermal components found in lookup table

real inclA       # inclusive area of pupil, must have been put into the
                 #thruput_er stack earlier by eg. putclr("area", real)....

pointer omodenkd # 'naked' obs mode without 'band' token: 
pointer tmodenkd
                 # I don;t understand (anands@stsci.edu)
pointer tbtopn()
pointer octptr   # table pointer for obs comp table
pointer tctptr   # table pointer for thl comp table

pointer octfn    # file name from template for obs comp table
pointer tctfn    # file name from template for thl comp table

pointer onparam, oparlist, oflist  # prolly only oflist used?
pointer tnparam, tparlist, tflist  # prolly only tflist used?

pointer thermspec  #The output (accumulated)spectrum before conversion to photons
pointer fullpass #Fake bandpass to give to rate function

pointer modenum, modelist
int mxlist
int nmode


pointer sp

#string  compname  "COMPNAME"
#string  keyname   "KEYWORD"
#string  tcompname "THCOMPNAME"
#string  tkeyname  "KEYWORD"


double scale, apscale
int apx, apy
begin

mxlist = 5


    call smark(sp)

    #all salloc (octptr,    SZ_FNAME,               TY_CHAR)
    #all salloc (tctptr,    SZ_FNAME ,              TY_CHAR)
    call salloc (octfn,     SZ_FNAME ,              TY_CHAR)
    call salloc (tctfn,     SZ_FNAME ,              TY_CHAR)

    call salloc (omodenkd,  SZ_LINE,                TY_CHAR)
    call salloc (tmodenkd,  SZ_LINE,                TY_CHAR)


    call salloc (onparam,   MAXPARAM,               TY_INT)
    call salloc (tnparam,   MAXPARAM,               TY_INT)

    call salloc (oparlist,  MAXPARAM*MAXLIST,       TY_REAL)
    call salloc (tparlist,  MAXPARAM*MAXLIST,       TY_REAL)

    call salloc (oflist,    MAXLIST*(SZ_FNAME+1),   TY_CHAR)
    call salloc (tflist,    MAXLIST*(SZ_FNAME+1),   TY_CHAR)

    call salloc (modenum, mxlist, TY_INT)
    call salloc (modelist, mxlist*(SZ_KEYWRD+1), TY_CHAR)

    call salloc (thermspec, nwave, TY_REAL)
    call amovkr (0.0, Memr[thermspec], nwave)

    call salloc (fullpass, nwave, TY_REAL)
    call amovkr (1.0, Memr[fullpass], nwave)

    call get_hstarea(inclA)  

    # Produce a list of table names from the obs/thl mode: these are
    # component file names for each of the optical/pupil
    # components in the optical train. We do this once for the optical
    # and once for the thermal components.
 
    call lastfile (ocmptb, Memc[octfn], SZ_FNAME)                                  
    call lastfile (tcmptb, Memc[tctfn], SZ_FNAME)                                  

	octptr = tbtopn (Memc[octfn], READ_ONLY, 0)
	tctptr = tbtopn (Memc[tctfn], READ_ONLY, 0)

    if (verbose) {
		call printf("\n\t obscol = %s\n") 
		call pargstr(obscol)
		call flush(STDOUT)

		call printf("\tomode    = %s\n") 
		call pargstr(omode)
		call flush(STDOUT)

   }

   call getnaked (tmode, Memc[tmodenkd], SZ_LINE)

    call getnaked (omode, Memc[omodenkd], SZ_LINE)
   if (verbose) {
#		call printf("\tomodenkd = %s\n") 
#		call pargstr(Memc[omodenkd])
#		call flush(STDOUT)
#
#		call printf("\ttmodenkd = %s\n") 
#		call pargstr(Memc[tmodenkd])
#		call flush(STDOUT)
#
#		call printf("ogrftb = %s\n\t") 
#		call pargstr(ogrftb)
#		call flush(STDOUT)
#
#		call printf("compname = %s  ") 
#		call pargstr(compname)
#		call flush(STDOUT)
#
#		call printf("keyname = %s  ") 
#		call pargstr(keyname)
#		call flush(STDOUT)
#
#		call printf("ocmptb = %s\n") 
#		call pargstr(ocmptb)
#		call flush(STDOUT)
#
#		call printf("oncomp = %d\n") 
#		call pargstr(oncomp)
#		call flush(STDOUT)


 
####Search the graph table for the optical components
	call printf("\n\t Optical components\n")
	call flush(STDOUT)
    }
    call searchgraf (verbose, ogrftb, GRF_COMPID,
			ocmptb, Memc[omodenkd],
                        MAXLIST, MAXPARAM, SZ_FNAME,  
                        oncomp, Memi[onparam], Memr[oparlist],
                        Memc[oflist]) 

    if (verbose) {
		call printf("\tobs oncomp = %d\n")
		call pargi(oncomp)
		call flush(STDOUT)

#		call printf("\tpflist(Memc[oflist], oncomp)\n")
#		call flush(STDOUT)
#		call pflist(Memc[oflist], oncomp)


		call printf("ocmptb = %s\n") 
		call pargstr(ocmptb)
		call flush(STDOUT)

#		call printf("ogrftb = %s\n\t") 
#		call pargstr(ogrftb)
#		call flush(STDOUT)


####Search the graph table for the thermal components
	call printf("\n\t Thermal components\n")
	call flush(STDOUT)
    }
    call searchgraf (verbose, ogrftb, GRF_THMLID, 
			tcmptb, Memc[omodenkd],
            MAXLIST, MAXPARAM, SZ_FNAME,  
            tncomp, Memi[tnparam], Memr[tparlist],
            Memc[tflist]) 

    if (verbose) {
		call printf("\tobs tncomp = %d\n")
		call pargi(tncomp)
		call flush(STDOUT)

#		call printf("\tpflist(Memc[tflist], tncomp)\n")
#		call flush(STDOUT)
#		call pflist(Memc[tflist], tncomp)


		call printf("tcmptb = %s\n\n") 
		call pargstr(tcmptb)
		call flush(STDOUT)

    }
#......................................................................................
# See about getting the special thermal mode keywords.
	call salloc (modelist, mxlist*(SZ_KEYWRD+1), TY_CHAR)#

#	# Break mode string into list of keywords and parameters

	
	call breakmode (Memc[omodenkd], MAXLIST, MAXPARAM, SZ_KEYWRD, nmode, 
			Memi[tnparam], Memr[tparlist], Memc[modelist])
		
#.......................................................................................
# Get pixel scale, for use in computing etendue.

    call getscale ( detcat, omode, apscale, apx, apy) 
    scale = apscale * 3600.0 #convert degrees to arcsec
    if (verbose) {
    call eprintf(" Detector lookup table, mode: %s %s \n")
    call pargstr(detcat)
    call pargstr(omode)
    call eprintf("\t Apscale is %g (degrees per pixel)\n")
    call pargd(apscale)
    call eprintf("\t Pixel scale is %g (arcsec per pixel)\n")
    call pargd(scale)
    call eprintf("\t Array size %d %d \n \n")
    call pargi(apx)
    call pargi(apy)
    }

# Loop over the components
    call tcptloop(oncomp, Memi[onparam], Memr[oparlist], Memc[oflist],  #Input optical
               tncomp, Memi[tnparam], Memr[tparlist], Memc[tflist],    #Input thermal
	       Memc[modelist], nmode,                                  #Input thermal
               wave, nwave, scale, inclA, verbose,		       #Input spectrum
               Memr[thermspec])					#output spectrum     

#-------------------------------------------------------------------------------
# Lengthy digression about units: Nov 04
#...............................................................................
# Working notes used when adding support for output units.
# 
# Tcptloop returns  thermspec, in the units: 
#   photons per sec per arcsec^2 per cm^2 per wavelength.

# This is very much like the definition of photlam:
#  photlam = photons per sec per cm^2 per angstrom
#
# It is safe to call phottoany in order to convert to
# the desired units; the per arcsec^2 carries through. This makes sense
# in principle, and experiments to play with this gave the expected results.
# 
# **** NOTE that phottoany(), when converting to counts, ALSO integrates
#      over the wavelength and the area. Thus we do not have to do that
#      integral if the requested form is counts.
# -- Note also, that because we want to return the rate in the requested
#    units, we shouldn't do that integral over wavelenght, or multiply
#    by the HSTarea, in any other case either. 
#
#............................................................................................
# Original implementation, now superseded by above logic:
#............................................................................................
# So we must:
#   - integrate over wavelength => photons per sec per arcsec^2 per cm^2
#   - multiply by hst area => photons per sec per arcsec^2
#   - multiply by pixel scale of detector => photons per sec per pixel
#
#    trate = syntegral( nwave, wave, Memr[thermspec] )  * inclA * scale**2 #NIC2 pixel size
#............................................................................................
# Convert to the desired units. Returns "<any> per square arcsec".

    call phottoany(form, nwave, wave, Memr[thermspec])

# Write thermal background spectrum to output file.

    call putspec (output, form, ogrftb, ocmptb, tcmptb, 
                  omode, nwave, wave, Memr[thermspec] )

# Sum up the spectrum, and give the answer per pixel (instead of  per square arcsec)
#------------------------------------------------------------------------

    trate = asumr(Memr[thermspec],nwave) * scale**2



    #Write scalar results to STDOUT

    call printf("trate is %8.4g flux units (%s) per pixel \n")
    call pargr(trate)
    call pargstr(form)
    call flush(STDOUT)

   #Write scalar results to output parameter too
    call clputr("thermflux",trate)

     call sfree(sp)

end
