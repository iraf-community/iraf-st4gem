# Thermback.par
#  obsmode, s, a, , , ,                            Instrument observation mode(s)
#  output, s, h, none, , ,                         Output table name
#  tcmptb, s, h, "mtab$*_tmt.fits", , ,            Instrument thermal component table
#  detcat, s, h, "mtab$detectors.dat", , ,         
#  wavetab, s, h, " ", , ,                         Wavelength table name
#  refdata, pset, h, , , ,                         Reference data
#  verbose, b, h, yes, , ,                         Print results to STDOUT ?
#  mode, s, h, 'a', , , 
#  form,s,h,counts,photlam|counts|flam|fnu|photnu|jy|mjy|\
#                  abmag|stmag|vegamag|obmag,,     Form for magnitude

# THERMBACK -- Compute the thermal background count rate of the HST
#
#* HISTORY *
#* A.Sivaramakrishnan  nn-Mmm-01   original
#* V.Laidler           18-Nov-04   add support for output units and write
#*                                    thermal background rate to output param

procedure thermback

#-- task parameters

pointer obsmode     # observation mode (optical synphot style)
pointer thlmode     # thermal mode: special requests on temperatures of 
                    # pupil components

pointer ogrftb      # graph table for 'optical' optics train

pointer ocmptb      # component lookup table for 'optical' optics train
pointer tcmptb      # ibid for optical elts and pupil components
                    # for thermal background contributions

pointer detcat      # table containing pixel scales of variuos detectors

pointer wavetab     # table of wavelengths at which calculations are performed
pointer	output	    # name of output table

pointer form        # units for output table & answer

bool    verbose     # print results to STDOUT ?

#-- utility stuff

int	tbtacc()
string	blank        ""

real    hstarea     # telescope area in cm^2
real    trate       # nr of thermal background photons detected per sec per steradian.

pointer pcode

pointer wave
int     nwave

pointer sp
real    clgetr()
bool	clgetb()

begin

    # Allocate dynamic memory for strings
    #
    call smark (sp)


    call salloc (ogrftb,      SZ_FNAME, TY_CHAR)
    call salloc (ocmptb,      SZ_FNAME, TY_CHAR)
    call salloc (tcmptb,      SZ_FNAME, TY_CHAR)
    call salloc (obsmode,     SZ_FNAME, TY_CHAR)
    call salloc (thlmode,     SZ_FNAME, TY_CHAR)
    call salloc (detcat,      SZ_FNAME, TY_CHAR)
    call salloc (output,      SZ_FNAME, TY_CHAR)
    call salloc (form,     SZ_FNAME, TY_CHAR)

    call salloc (wavetab,     SZ_FNAME, TY_CHAR)
    call salloc (pcode,       SZ_LINE,  TY_INT)
 

    # Read task parameters
    verbose =  clgetb ("verbose")
    
    hstarea = clgetr ("area")  # get from refdata pset
    call put_hstarea (hstarea) # push this into where?

    call clgstr ( "grtbl",        Memc[ogrftb],        SZ_FNAME)
    if (verbose) {
	call printf("grtbl parameter is: %s\n")
	call pargstr(Memc[ogrftb])
	call flush(STDOUT)
    }

    call clgstr ("cmptbl",        Memc[ocmptb],        SZ_FNAME)
    call clgstr ("tcmptb",        Memc[tcmptb],        SZ_FNAME)
    call clgstr ("form",        Memc[form],        SZ_FNAME)

    call clgstr ("obsmode",       Memc[obsmode],       SZ_FNAME)
#    call clgstr ("thlmode",       Memc[thlmode],       SZ_FNAME)
    call clgstr ("detcat",        Memc[detcat],       SZ_FNAME)
  
    call clgnone ("output", Memc[output], SZ_FNAME)

    call clgstr ("wavetab",       Memc[wavetab],       SZ_FNAME)
    #Memc[wavetab] = EOS

    call inisyntab

    call expcompile (Memc[obsmode], Memi[pcode], SZ_LINE)


     # Read or compute the wavelength set.

     if (tbtacc (Memc[wavetab]) == NO) {
	    call getwavelen (blank, Memc[ogrftb], Memc[ocmptb], 
			     Memi[pcode], 1, SZ_LINE, wave, nwave)
     } else {
	    call getwavelen (Memc[wavetab], Memc[ogrftb], Memc[ocmptb], 
			     Memi[pcode], 1, SZ_LINE, wave, nwave)
     }

    if (verbose){
	call eprintf("wavetab: %s \n")
	call pargstr(Memc[wavetab])
	call eprintf("getwavelen: nwave = %d\n")
	call pargi(nwave)
    }

    call thermal(Memc[obsmode], Memc[thlmode], Memc[ogrftb],
                              Memc[ocmptb], Memc[tcmptb], Memc[detcat],
			      Memc[form],
                              Memr[wave], nwave, Memc[output], trate, 
                              verbose)

    call clssyntab
    call sfree(sp)


end
