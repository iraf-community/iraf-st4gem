#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

# TCPTLOOP -- Compute the thermal background photon rate -- wavelength loop

#* HISTORY *
#* A.Sivaramakrishnan      Sep 2001       Created
#* V.Laidler		   Jan 2002	  Completed

#include "libsynphot.h"

define  MAXLIST     100
define  MAXPARAM    3
define  SZ_KEYWRD       32              # Length of keyword

# above defines are also in libsynphot.h under st4gem$lib/synphot
# but convention is to include them explicitly.
  
procedure tcptloop(oncomp, onparam, oparlist, oflist,  #Input optical
                   tncomp, tnparam, tparlist, tflist,  #Input thermal
	           modelist, nmode,                                  #Input thermal 
                   wave, nwave, scale, inclA, verbose,		#Input spectrum, verbosity control
                    thermspec)                              #Output result

int    oncomp                         # i: number of synclassic cpts in train
int    onparam[ARB]                   # i: nr of parameters per o_keyword (classic)
real   oparlist[MAXPARAM, ARB]        # i: parameters of each (classic) keyword
char   oflist[SZ_FNAME, MAXLIST]      # i: list of (classic) component file names
#
int    tncomp                         # i: number of synthermal cpts in train
int    tnparam[ARB]                   # i: nr of parameters per t_keyword (thermal)
real   tparlist[MAXPARAM, ARB]        # i: parameters of each (thermal) keyword
char   tflist[SZ_FNAME, MAXLIST]      # i: list of (thermal) component file names
char   modelist[SZ_KEYWRD, MAXLIST]   # i: list of modes matching to tnparam and tparlist
int    nmode			      # i: dimension of tnparam, tparlist, and modelist
#
real   wave[ARB]                      # i: wavelengths --- domain of calculation
int    nwave                          # i: cardinality of domain
double scale                          # i: pixel scale of detector
real   inclA      		      # i: inclusive area of pupil, must have been put into the

bool   verbose	                      # i: controls verbosity
#
real	thermspec[ARB]			# o: output thermal spectrum


#
		
#--

real   syntegral()   # function call
bool   strne(), streq()
bool   getkeyval()
int    ahivi()



pointer sp	 # stack memory management

# Global scalars


int  icomp,i	#loop counter
int  oindex     # index into the optical component tables. Note ASSUMPTION:
		#      OPTICAL COMPONENTS OCCUR IN THE SAME ORDER AS THERMAL COMPONENTS!
real   trate    # detected thermal photons s-1 sr-1

bool  checktemp # do some components have overridden temperatures

# Global vectors


pointer unity_array # Sometimes we need one of these
pointer thruput    # A temporary variable to hold the optical thruput of transmissive components


# Per-component scalars
char thlmode[SZ_FNAME]       # Component name
real beamFillingFraction  # Fraction of beam filled by component
real newtemp, temp	 # temperature, gotten from refdata or (later) overridden
char thtype[8]	  # legal values CLEAR, OPAQUE, THRU, NUMERIC
	string CLEAR   "CLEAR"
	string OPAQUE  "OPAQUE"
	string THRU    "THRU"
	string NUMERIC "NUMERIC"



# Per-component vectors
pointer emiss     # emissivity on wavetab
pointer BB        # black body function (handily computed at T on wavetab)
pointer cptThermSpec  # contribution from a component

begin

    call smark(sp)

#call eprintf("inside tcptloop\n")

    # Initialize global vectors

    call salloc (thruput, nwave, TY_REAL)
    call amovkr (0.0, Memr[thruput], nwave)

    call salloc (unity_array, nwave, TY_REAL)
    call amovkr (1.0, Memr[unity_array], nwave)

    # Initialize per-component vectors
    call salloc (emiss, nwave, TY_REAL)
    call amovkr (0.0, Memr[emiss], nwave)

    call salloc (BB, nwave, TY_REAL)
    call amovkr (0.0, Memr[BB], nwave)

    call salloc (cptThermSpec, nwave, TY_REAL)
    call amovkr (0.0, Memr[cptThermSpec], nwave)
 
    oindex = 0          
    #.................................................................................
    #FRILL - Do this later
    # Determine which default component temperatures to override: arrays are
    # of size tncomp (have to figure out how to map comp name in the tmode
    # to a component number in the thermal properties component file list --
    # call getTlist(tflist, tmode, templist, ntemps)
    #             # tflist           # i: thermal components' file list
    #            # string tmode        # i: thermal mode string
    #            # tncomp               # 1: number of overridden temps
    #            # real templist[ARB]  # o: temp of overridden component 
    #.................................................................................          

#...................................................................................	  
    if (verbose) {	
	call eprintf("\n")
	for ( i = 1; i<= nmode; i=i+1 ) {
	  call eprintf(" Keyword %d has %d tnparams; tparlist %g ; modelist %s \n")
	  call pargi(i)
	  call pargi(tnparam[i])
	  call pargr(tparlist[1,i])
	  call pargstr(modelist[1,i])
        }
	call eprintf("\n")
     }
#...................................................................................	   	
  
    # Prepare the temperature override dictionary
    if (ahivi(tnparam,nmode)>0) 
    {
	checktemp = true
	do i = 1, nmode 
	  {	  call strlwr(modelist[1,i]) }
    }

#...................................................................................	   	
    if (verbose) {
	   call eprintf(" nwave %d ; wavelength range %g %g \n \n")
	   call pargi(nwave)
	   call pargr(wave[1])
	   call pargr(wave[nwave])
    }
    # Open all obs_comp files and get array of component names into obs_compname[]

    # Loop over components in thermal component file list (icomp loop variable):
    do icomp = 1, tncomp {

    # Get the thermal data for this component: Visible area, temperature, and mode (type)

      call thcompkey(tflist[1,icomp], thlmode, beamFillingFraction, temp, thtype)  
#      if (verbose) {
#	call eprintf("Loop iteration %d...\n")
#	call pargi(icomp)
#
#	call eprintf(" Hdr: component, beamFillingFraction, temp, thtype: %s %g %g %s \n")
#	call pargstr(thlmode)
#	call pargr(beamFillingFraction)
#	call pargr(temp)
#	call pargstr(thtype)
#
#	call eprintf("tnparam %g\n")
#	call pargi(tnparam[icomp])
#      }

    #Check to be sure we need to do anything.
	if ( strne(thtype, CLEAR) ) {  


# TRANSMISSIVE SECTION: OPAQUE COMPONENTS DON'T HAVE THIS........................................
#  .............Numeric components also reduce the beam so they can be treated this way too.

        # If the component is transmissive, look up the optical data...
          if (streq(thtype, THRU) || streq(thtype, NUMERIC) ) {
		oindex = oindex + 1  #Check assumption here? Where are the names?
                # get component thruput table data  using obs_compname[] array 
#             if (verbose) { 
#               call eprintf(" %s is a transmissive component \n")
#	       call pargstr(thlmode)
#               call eprintf(" Filename %s \n")
#               call pargstr(oflist[1,oindex])
#             }	

               call evalfilt(oflist[1,oindex], onparam,
                   oparlist[1,oindex], nwave, wave, YES, 
		   Memr[thruput], Memr[unity_array])   # (Think I got this right)

        #... & reduce the beam.
                call amulr(thermspec, Memr[thruput], thermspec, nwave)
           } #end if


#THERMAL SECTION: ALL PHYSICAL COMPONENTS HAVE THIS..................................................
	   if ( strne(thtype, NUMERIC) ) {
         #Frill...........................................................
	 # Obtain default T in DEFTEMP keyword 
         # if (component name in the T-overridden table) 
         #        # overwrite T with user-requested T
         #.................................................................
	     if (checktemp) {

		if (getkeyval(thlmode, newtemp, modelist, tparlist, nmode) ) {
		   call eprintf(" Overriding %s temperature to %g \n")
		   call pargstr(thlmode)
		   call pargr(newtemp)
		   temp = newtemp
		}
	     }


        # Read emissivity table(icomp) and interpolate onto wavetab
             call evalemiss( tflist[1,icomp], tnparam, tparlist[1,icomp], 
 			   nwave, wave, YES, 
                           Memr[emiss], Memr[unity_array])

	 # Compute BB function on wavetab

#	   call bbfunc(temp, nwave, wave, Memr[BB])

# BBFUNC returns units that we don't understand - so BBSI should become a permanent part
# of the library. BBSI needs renaming: it really returns photlam per square arcsec.

	     call bbsi(temp, nwave, wave, Memr[BB])


         # Component spectrum = bff * BB * emissivity:
#	   call printf("\n\t BB 1, N: %g %g ...\n\t")
#	   call pargr(Memr[BB])
#	   call pargr(Memr[BB+nwave-1])
#	   call flush(STDOUT)

	     call amulr(Memr[BB], Memr[emiss], Memr[emiss], nwave)
	     call amulkr(Memr[emiss], beamFillingFraction, Memr[cptThermSpec], nwave)

# 
#
#	   call printf("\n\t BB*emiss*bff 1, N: %g %g ...\n\t")
#	   call pargr(Memr[cptThermSpec])
#	   call pargr(Memr[cptThermSpec+nwave-1])
#	   call flush(STDOUT)

           
	   
# ADD CONTRIBUTION OF THIS COMPONENT
           # Finally, add the component contribution to the beam.
	     call aaddr( thermspec, Memr[cptThermSpec], thermspec, nwave) 

	   } #end check for not numeric

        } #end check for clear component

# This calculation does NOT belong here, but is included so that we can print out
# the results as each component is added.
     if (verbose) {

      trate = syntegral( nwave, wave, thermspec )  * inclA * scale**2 # pixel size
      call printf("Component %d : %s  (%s, area %g, temp= %g ) \n")
      call pargi(icomp)
      call pargstr(thlmode)
      call pargstr(thtype)
      call pargr(beamFillingFraction)
      call pargr(temp)
      call flush(STDOUT)

      call printf("\tCounts per second per pixel: %g \n \n ")
      call pargr(trate)
      call flush(STDOUT)
     }
    } #end loop over all components
#.......................................................................................
     if (verbose) {
	   call eprintf(" nwave %d \n")
	   call pargi(nwave)
     }

    #clean up when finished
    call sfree(sp)           

end
