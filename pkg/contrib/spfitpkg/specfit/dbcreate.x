###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	dbcreate
#
#  Description:	DBCREATE is an Iraf task that helps you create specfit database
#		files
#
#  Arguments:	New database name
#
#  Returns:	Output is to user specified database
#
#  Notes:	Information shared in common blocks defined in "specfit.com".
#
#  History:	August 1994	Created grimes
#		June 1995	grimes 	modified to add comment ability
#
###########################################################################


include "specfit.h"


procedure t_dbcreate()






real 	pars[PARPERC],lower[PARPERC]
real  	upper[PARPERC], st[PARPERC]
real 	tol[PARPERC]
int 	fix_type[PARPERC]

common /new_comp_parts/ pars,lower,upper,st,tol,fix_type


char 	compreal[SZ_CNAME,NCOMP], dbcomments[SZ_LINE,NCOMP]
int 	parreal[NCOMP]

int	done, correct
	
int 	i,k, count,j, position, type, junk

char	prompt[SZ_LINE], tempstr[SZ_FNAME], clstring[SZ_LINE]
char	datanam[SZ_FNAME],datafile[SZ_FNAME]
char 	what[3]

pointer dt


pointer	dtmap1()
int 	clgeti(), clgstr(), strcmp(), strlen()
#bool	clgetb()


include "specfit.com"

begin
	#Initialize components information
	call sfinitcomp()
	
	position = 1

	junk = clgstr("namedata",datanam,SZ_FNAME-1)

	what[1] = datanam[1]
	what[2] = datanam[2]
	what[3] = EOS
	# Check for sf prefix
	if ( strcmp("sf",what) == 0 ) {
		call sprintf(datafile,SZ_FNAME-1,"%s")
			call pargstr(datanam)
		do i = 1, SZ_FNAME {
			datanam[i] = EOS
		}

		do i = 3, strlen(datafile) {
			datanam[i-2] = datafile[i]
		}

		call printf("WARNING: sf header already present\n")
		call printf("Will use database name %s instead\n")
			call pargstr(datanam)

	} else {
		call sprintf(datafile,SZ_FNAME-1,"sf%s")
			call pargstr(datanam)
	}

	# Open file and start adding info
	dt = dtmap1("./",datafile,APPEND)

	call dtptime(dt)
	call dtput(dt,"begin\t%s\n")
		call pargstr(datanam)
	call dtput(dt,"\ttask\tspecfit\n")

count = 0
done = NO
#Main loop to add components
while ( done == NO ) {

	junk = clgstr( "comp_name" , compreal[1,count+1] , SZ_LINE + 1 )
	#Quit
	if ( strcmp(compreal[1,count+1],"quit") == 0 ) {
		break
	}
	#Print out list of component types
	if ( strcmp(compreal[1,count+1],"?") == 0 ) {
		call printf("The known component types are:\n")
		do i = 1, NKEYS {
			call printf("\t%s\n")
			  call pargstr(compkeys[1,i])
		}
		next
	}

	count = count + 1

	type = 500
	do i = 1, NKEYS {
		if ( strcmp(compreal[1,count],compkeys[1,i]) == 0 ) {
			type = i
			parreal[count] = ncpar[i]
			break
		}
	}
	#If unknown component type
	if ( type == 500 ) {
		call printf("Unknown component type %s\n")
			call pargstr(compreal[1,count])
		correct = clgeti("parambool")
		if ( correct == YES ) {
			parreal[count] = clgeti("parameters")
		} else {
			count = count - 1
			next
		}
	}

	call strcpy("",dbcomments[1,count],SZ_LINE)
	call clpstr("comment","")
	junk = clgstr("comment",dbcomments[1,count],SZ_LINE)

 switch (type) { #Go through all the components and get their values

	case 1: #linear
		call get_rparam("Flux at 1000 Angstroms",1)

		call get_rparam("Continuum slope(change in flux/A)",2)

	case 2: #powerlaw
		call get_rparam("Flux at 1000 Angstroms",1)
	
		call get_rparam("Power Law Index Alpha for f=(x/1000)^(-alpha)",2)

	case 3: #black body
		call get_rparam("Flux at 5500 Angstroms (Flambda)",1)

		call get_rparam("Temperature (Kelvin)",2)


	case 4: #gaussian
		call get_rparam("Enter the Gaussian flux",1)

		call get_rparam("Enter the centroid",2)

		call get_rparam("Enter the FWHM",3)
		
		call get_rparam("Enter the skew (1 is symmetrical)",4)




	case 5: #logarith
		call get_rparam("Flux, or area under the line",1)

		call get_rparam("Enter the centroid",2)

		call get_rparam("FWHM of the line in km/s",3)
		
		call get_rparam("Skew (1=symmetric)",4)

	case 6: #labsorp
		call get_rparam("Equivalent width of the line",1)

		call get_rparam("Enter the centroid",2)

		call get_rparam("FWHM in km/s",3)

	case 7: #eabsorp 
		call get_rparam("optical depth at the edge",1)
 	
		call get_rparam("wavelength of the edge",2)

		
	case 8: #recomb 
		call get_rparam("Flux at the edge",1)

		call get_rparam("Wavelength of the edge",2)
	
		call get_rparam("Electron Temperature (K)",3)

		call get_rparam("FWHM",4)

	case 9: #extcor
		call get_rparam("Mean Galactic Extinction E(B-V)",1)


	case 10: #usercont
		call get_rparam("Normalization (1 if model fluxes match data)",1)

		call get_rparam("Linear shift in wavelength",2)

		call get_rparam("Redshift",3)

		call get_rparam("Value of the key parameter",4)



	case 11, 12: #userline,userabs
		call get_rparam("Normalization (1.0 file input values)",1)

		call get_rparam("Linear shift in wavelength",2)

		call get_rparam("Redshift",3)

		call get_rparam("Value of the key parameter",4)


	
	case 13: #logabs
		call get_rparam("Optical Depth at Line Center",1)

		call get_rparam("Enter the centroid",2)

                call get_rparam("FWHM of the line in km/s",3)


	case 14: #lorentz
		call get_rparam("Flux, or area under the line",1)

		call get_rparam("Enter the centroid",2)

		call get_rparam("FWHM of the line in km/s",3)

		call get_rparam("Alpha -- explonent of centroid",4)


	case 15: #dampabs

		call get_rparam("Column density times the oscillator strength",1)

		call get_rparam("Enter the centroid",2)

		call get_rparam("3 Value",3)
	
        
	case 16: #bpl
                call get_rparam("Flux at the Break Wavelength",1)

                call get_rparam("The Break Wavelength",2)

                call get_rparam("Power Law Index Above the Break",3)

                call get_rparam("Power Law Index Below the Break",4)

	case 17: #ffree
		call get_rparam("Normalization at 5500A",1)

		call get_rparam("Temperature in Kelvins",2)

	case 18: #tauabs

		call get_rparam("Optical depth at line center",1)

		call get_rparam("Enter the centroid",2)

		call get_rparam("FWHM in km/s",3)

	case 19: #extdrude

		call get_rparam("Mean Galactic Extinction E(B-V)",1)

		call get_rparam("w0",2)
		
		call get_rparam("Gamma",3)

		call get_rparam("C1",4)

		call get_rparam("C2",5)
		
		call get_rparam("C3",6)

		call get_rparam("C4",7)

	case 20: #disk

		call get_rparam("Normalized Flux",1)

		call get_rparam("Beta Value",2)

		call get_rparam("Temperature in Kelvin",3)


	case 21:	#ccmext

		call get_rparam("EBV Value E(B-V)",1)

		call get_rparam("RV",2)

#If Adding your own fitting function and also putting it in interactive mode
#	case 22: WHat number of component type is it
#
#
#		To enter a parameter just use this
#			
#		call get_rparam("A Prompt String",x)
#
#		This will automatically let the user input all the relevant 
#			information.  The Prompt String is just a short 
#			description of the parameter that will be printed 
#			to the user and x is again the parameters number.  It 
#			will automatically ask the user for the parameter ( 
#			using the prompt string ), a lower limit, upper limit, 
#			step size, tolerance, or fix/free status
#
#		If you have questions look at the above cases to see examples
#

	default:
		if ( type < NKEYS ) {
			do  k = 1, ncpar[type] {
				call sprintf(prompt,SZ_LINE,"Enter Parameter #%d\n")
					call pargi(k)
				call get_rparam(prompt,k)

			}
		} else {
			do k = 1, parreal[count] {
				call sprintf(prompt,SZ_LINE,"Enter Parameter #%d")
					call pargi(k)
				call get_rparam(prompt,k)
			}
		}
	}

	do i = 0, parreal[count]-1 {
		par0[position+i] = pars[i+1]
		blim[position+i] = lower[i+1]
		tlim[position+i] = upper[i+1]
		step[position+i] = st[i+1]
		ptol[position+i] = tol[i+1]
		ifix[position+i] = fix_type[i+1]				
	}


	position = position + parreal[count]



}

	call dtput(dt, "components\t%d\n")
		call pargi(count)
	for ( i = 1; i <= count; i = i + 1) {
		call dtput(dt, "\t\t%s\n")
			call pargstr(compreal[1,i])
	}


# Write out parameters for each component
	k = 0
	for ( i = 1; i <= count; i = i + 1) {
		call sprintf (tempstr, SZ_FNAME, "%s%d")
		    call pargstr(compreal[1,i])
		    call pargi(i)
		if ( strcmp(dbcomments[1,i],"")!=0 ) {
			call dtput(dt, "\t\t%s\t%d #%s\n")
		    		call pargstr(tempstr)
		    		call pargi(parreal[i])
				call pargstr(dbcomments[1,i])
		} else {
			call dtput(dt, "\t\t%s\t%d\n")
				call pargstr(tempstr)
				call pargi(parreal[i])
		}
		for ( j = 1; j <= parreal[i]; j = j + 1) {
			k = k + 1
			call dtput(dt, "\t\t\t%g %g %g %g %g %d\n")
			    call pargr(par0[ k ])
			    call pargr(blim[ k ])
			    call pargr(tlim[ k ])
			    call pargr(step[ k ])
			    call pargr(ptol[ k ])
			    call pargi(ifix[ k ])
		}
	}

	call dtunmap (dt)

#Call utility dbcheck to see if any problems

	if ( count == 0 ) {
		call printf("No components inputted: insertion aborted\n")
	} else {
		call printf("\n\nCalling utility dbcheck\n\n")

		call sprintf(clstring,SZ_LINE,"dbcheck %s")
			call pargstr(datanam)
		call clcmdw(clstring)	
	}
	
end





