#############################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#                        
#  Synopsis:    procedure sfcomponent_add(nfree,fpar)
#		procedure sfcomponent_del(nfree,fpar)
#		procedure sfcomponent_list()
#                 
#                
#  Description: These are used only in interactive mode.  They allow the user
#		to insert new components, delete existing components, or list 
#		all of the existing components.
#
#
#  Arguments:   
#
#
#  Returns:     
#
#
#  History:     June	1994	J. Grimes
#		August 4 94	J Grimes, Added extcor, default interactive
#					parameter algorithmn, and fixed minor
#					bug about number of free parameters
#		January 95	grimes added ccmext
#		2/16/95		grimes fixed linked params 
#
#
################################################################   

include "specfit.h"
include	<gio.h>
include <gset.h>


int procedure sfcomponent_add(nfree,fpar,gfd,gt)

int 	nfree
real 	fpar[ARB]
pointer	gfd,gt

int 	comp_num, type,junk
char 	comp_name[SZ_LINE]

 
real 	pars[PARPERC],lower[PARPERC]
real  	upper[PARPERC], st[PARPERC]
real 	tol[PARPERC]
int 	fix_type[PARPERC]

common /new_comp_parts/ pars,lower,upper,st,tol,fix_type
	
int 	success
int 	i,j,k
real 	chisq
char 	strval[4]
char	prompt[SZ_LINE]

include "specfit.com"

int 	clgeti(), clgstr()

begin
#Get name of new component, its type number, and its position
 
	junk = clgstr( "comp_name" , comp_name , SZ_LINE + 1 )

	call cmatch(comp_name,type,0) #Find identifying number of that type
	if ( type == -1 ) {
		call printf("Component insertion aborted\n")
		return NO
	}

	comp_num = clgeti("comp_number")
	call printf("Return to the Text Screen.\n")

 success = YES


 call gdeactivate(gfd,0)

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
		

		call printf("Exit back to the graphics screen to input the centroid value.\n")
		call greactivate(gfd,AW_CLEAR)		
		call gclear(gfd)
		call replot(gfd,gt,Memr[lambdas],Memr[spectrum],npts)
		call overplot(gfd,gt,Memr[lambdas],Memr[fitsp],npts)
		call printf("Place the cursor at the centroid & press return\n")
		junk = 4
		call clgcur("cursor",pars[2],junk,junk,junk,strval,junk)
		call printf("Return to Text Screen.\n")
		call gdeactivate(gfd,0)
		call get_rest_of_vals(2)

		call get_rparam("Enter the FWHM",3)
		
		call get_rparam("Enter the skew (1 is symmetrical)",4)

	case 5: #logarith
		call get_rparam("Flux, or area under the line",1)


		call printf("Exit back to the graphics screen to input the centroid value.\n")
		call greactivate(gfd,AW_CLEAR)		
		call gclear(gfd)
		call replot(gfd,gt,Memr[lambdas],Memr[spectrum],npts)
		call overplot(gfd,gt,Memr[lambdas],Memr[fitsp],npts)
		call printf("Place the cursor at the centroid & press return\n")
		junk = 4
		call clgcur("cursor",pars[2],junk,junk,junk,strval,junk)
		call printf("Return to Text Screen.\n")
		call gdeactivate(gfd,0)
		call get_rest_of_vals(2)





		call get_rparam("FWHM of the line in km/s",3)
		
		call get_rparam("Skew (1=symmetric)",4)

	case 6: #labsorp
		call get_rparam("Equivalent width of the line",1)
		
		call printf("Exit back to the graphics screen to input the centroid value.\n")
		call greactivate(gfd,AW_CLEAR)		
		call gclear(gfd)
		call replot(gfd,gt,Memr[lambdas],Memr[spectrum],npts)
		call overplot(gfd,gt,Memr[lambdas],Memr[fitsp],npts)
		call printf("Place the cursor at the centroid & press return\n")
		junk = 4
		call clgcur("cursor",pars[2],junk,junk,junk,strval,junk)
		call printf("Return to Text Screen.\n")
		call gdeactivate(gfd,0)
		call get_rest_of_vals(2)


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


		call printf("Exit back to the graphics screen to input the centroid value.\n")
		call greactivate(gfd,AW_CLEAR)		
		call gclear(gfd)
		call replot(gfd,gt,Memr[lambdas],Memr[spectrum],npts)
		call overplot(gfd,gt,Memr[lambdas],Memr[fitsp],npts)
		call printf("Place the cursor at the centroid & press return\n")
		junk = 4
		call clgcur("cursor",pars[2],junk,junk,junk,strval,junk)
		call printf("Return to Text Screen.\n")
		call gdeactivate(gfd,0)
		call get_rest_of_vals(2)

                call get_rparam("FWHM of the line in km/s",3)


	case 14: #lorentz
		call get_rparam("Flux, or area under the line",1)
		
		call printf("Exit back to the graphics screen to input the centroid value.\n")
		call greactivate(gfd,AW_CLEAR)		
		call gclear(gfd)
		call replot(gfd,gt,Memr[lambdas],Memr[spectrum],npts)
		call overplot(gfd,gt,Memr[lambdas],Memr[fitsp],npts)
		call printf("Place the cursor at the centroid & press return\n")
		junk = 4
		call clgcur("cursor",pars[2],junk,junk,junk,strval,junk)
		call printf("Return to Text Screen.\n")
		call gdeactivate(gfd,0)
		call get_rest_of_vals(2)

		call get_rparam("FWHM of the line in km/s",3)

		call get_rparam("Alpha -- explonent of centroid",4)


	case 15: #dampabs

		call get_rparam("Column density times the oscillator strength",1)
		
		call printf("Exit back to the graphics screen to input the centroid value.\n")
		call greactivate(gfd,AW_CLEAR)		
		call gclear(gfd)
		call replot(gfd,gt,Memr[lambdas],Memr[spectrum],npts)
		call overplot(gfd,gt,Memr[lambdas],Memr[fitsp],npts)
		call printf("Place the cursor at the centroid & press return\n")
		junk = 4
		call clgcur("cursor",pars[2],junk,junk,junk,strval,junk)
		call printf("Return to Text Screen.\n")
		call gdeactivate(gfd,0)
		call get_rest_of_vals(2)

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

		
		call printf("Exit back to the graphics screen to input the centroid value.\n")
		call greactivate(gfd,AW_CLEAR)		
		call gclear(gfd)
		call replot(gfd,gt,Memr[lambdas],Memr[spectrum],npts)
		call overplot(gfd,gt,Memr[lambdas],Memr[fitsp],npts)
		call printf("Place the cursor at the centroid & press return\n")
		junk = 4
		call clgcur("cursor",pars[2],junk,junk,junk,strval,junk)
		call printf("Return to Text Screen.\n")
		call gdeactivate(gfd,0)
		call get_rest_of_vals(2)
		
		call get_rparam("FWHM in km/s",3)

	case 19: #extdrude

		call get_rparam("Mean Galactic Extinction E(B-V)",1)

		call get_rparam("w0",2)
		
		call get_rparam("Gamma",3)

		call get_rparam("C1",4)

		call get_rparam("C2",5)
		
		call get_rparam("C3",6)

		call get_rparam("C4",7)

	case 20:		#disk

		call get_rparam("Normalized Flux",1)

		call get_rparam("Beta Value",2)
		
		call get_rparam("Temperature in Kelvin",3)


	case 21:  #ccmext

		call get_rparam("EBV Value E(B-V)",1)
	
		call get_rparam("RV",2) 


#If Adding your own fitting function and also putting it in interactive mode
#	case 22: WHat number of component type is it
#
#	For every parameter in your component you will have to let
#			the user input a value.  If you want the user to enter
#			the value as a cursor position use this. All you have  
#			to do is change x to the number of the parameter.  The
#			function get_rest_of_vals will let the user enter all
#			the parameters information ( i.e. lower and upper bounds
#			,tolerance, step size, and fix/free status )
#		
#		call printf("Exit back to the graphics screen to input the centroid value.\n")
#		call greactivate(gfd,AW_CLEAR)		
#		call gclear(gfd)
#		call replot(gfd,gt,Memr[lambdas],Memr[spectrum],npts)
#		call overplot(gfd,gt,Memr[lambdas],Memr[fitsp],npts)
#		call printf("Place the cursor at the centroid & press return\n")
#		junk = 4
#		call clgcur("cursor",pars[2],junk,junk,junk,strval,junk)
#		call printf("Return to Text Screen.\n")
#		call gdeactivate(gfd,0)
#		call get_rest_of_vals(2)
#
#
#
#
#	If the parameter is to be entered by the keyboard all you have 
#			to do is call
#			
#		call get_rparam("A Prompt String",x)
#
#		This will automatically let the user input all the relevant 
#			information.  The Prompt String is just a short 
#			description of the parameter that will be printed 
#			to the user and x is again the parameters number
#
#		If you have questions look at the above cases to see examples
#		
#	Also if you use the dbcreate task you might want to look at the file
#		dbcreate.x and do much the same as above except that its
#		impossible to do graphics input.
#



	default:
		if ( type <= NKEYS ) {
			do  k = 1, ncpar[type] {
				call sprintf(prompt,SZ_LINE,"Enter Parameter #%d\n")
					call pargi(k)
				call get_rparam(prompt,k)
			}
		} else {
			call eprintf("Component insertion aborted, illegal component\n")  
			success = NO
		}
	}

 call printf("End of Component Insertion\n")
 call greactivate(gfd,AW_CLEAR)

 if (success==YES) { 
	#If the value for a component have been entered sucessfully
	#Update all the global values



	ncomp = ncomp + 1

	parptr[1,ncomp] = parptr[ncpar[comtyp[ncomp-1]],ncomp-1] + 1

 	if (comp_num > ncomp) {
		comp_num = ncomp 
	}

	#Move old components back in the array from new comps position

	for (j=parptr[ncpar[comtyp[ncomp-1]],ncomp-1]+ncpar[type];j>=parptr[1,comp_num]+ncpar[type];j=j-1)
		{
		par0[j] = par0[j-ncpar[type]]
		blim[j] = blim[j-ncpar[type]]
		tlim[j] = tlim[j-ncpar[type]]
		step[j] = step[j-ncpar[type]]
		ptol[j] = ptol[j-ncpar[type]]
		ifix[j] = ifix[j-ncpar[type]]
		call strcpy(comments[1,j-ncpar[type]],comments[1,j],SZ_LINE)
		}

 	npar = npar + ncpar[type]


	#Put in new comps values
 	k = 0
 	for (j=parptr[1,comp_num];j<=parptr[1,comp_num]+ncpar[type]-1;j=j+1) 
		{
		k = k + 1
		par0[j] = pars[k]
		blim[j] = lower[k]
		tlim[j] = upper[k]
		step[j] = st[k]
		ptol[j] = tol[k]
		ifix[j] = fix_type[k]
		call strcpy("",comments[1,j],SZ_LINE)
		}

	#Also update other globals
 	for (i=ncomp;i>comp_num;i = i - 1) {
		comtyp[i] = comtyp[i-1]
		inplot[i] = inplot[i-1]
	  }
 	comtyp[comp_num] = type
	inplot[comp_num] = YES

	#Update the parptr array so that everything points to the correct
		# value
	k = 0
	for (i=1;i<=ncomp;i=i+1) 
		{
		for (j=1;j<=ncpar[comtyp[i]];j=j+1) 
			{
			k = k + 1
			parptr[j,i] = k
			}
		}

	for (i=1;i<=ncomp;i=i+1) {
		for (j=1;j<=ncpar[comtyp[i]];j=j+1) {
			if (ifix[parptr[j,i]]!=0 && ifix[parptr[j,i]]!=-1 && i!=comp_num) {
 				if (ifix[parptr[j,i]]>=comp_num)
					ifix[parptr[j,i]]=ifix[parptr[j,i]]+1
				}
			}
		}


	#Update fpar
	call setpar(nfree,fpar)
	call chispec(nfree, fpar, chisq)
 
 	}
	return SUCCESS
end


procedure get_rparam(pstring,num) #Used to input real parameter vals 
char pstring[ARB]
int num

real pars[PARPERC],lower[PARPERC]
real  upper[PARPERC], st[PARPERC]
real tol[PARPERC]
int fix_type[PARPERC]


common /new_comp_parts/ pars,lower,upper,st,tol,fix_type

real clgetr()

begin
	call printf("%s\t\t")
	  call pargstr(pstring)
	pars[num] = clgetr("new_par_value")
	call get_rest_of_vals(num)
end




procedure get_rest_of_vals(num) #Used in the above procedure to automate some of
				#the input procedures
int num


real pars[PARPERC],lower[PARPERC]
real  upper[PARPERC], st[PARPERC]
real tol[PARPERC]
int fix_type[PARPERC]

common /new_comp_parts/ pars,lower,upper,st,tol,fix_type

real low,up,step

real clgetr()
int clgeti()

begin
 #Get the default multipliers for determining the new values from a parameter
 low = clgetr("low_mult")
 up = clgetr("high_mult")
 step = clgetr("step_mult")

 call clputr("lower_li",low * pars[num])
 lower[num] = clgetr("lower_li")
 call clputr("upper_li",up * pars[num])
 upper[num] = clgetr("upper_li")
 call clputr("step_siz",step * pars[num])
 st[num] = clgetr("step_siz")
 tol[num] = clgetr("param_to")
 call clputi("fix_or_f",0)
 fix_type[num] = clgeti("fix_or_f")
end







procedure sfcomponent_del(nfree,fpar) #Delete an existing component
int 	nfree
real 	fpar[ARB]

include "specfit.com"

int 	num_comp,k,j,i

real 	chisq

int answer


int 	clgeti()
#bool 	clgetb()

begin
	#Get the components number
	call printf("Which component shall be deleted?\n")
	num_comp = clgeti("comp_num") #Find out which component to be deleted

	if (num_comp > ncomp )
		num_comp = ncomp


	call printf("Are you sure you want to delete %s number %d ")
		call pargstr(compkeys[1,comtyp[num_comp]])
		call pargi(num_comp)
	
	answer = clgeti("choice")
	if ( answer == NO ) {
		call printf("Component deletion aborted\n")
		return
	}


#Update all global values

	npar = npar - ncpar[comtyp[num_comp]]

	#Move the param and associated values into the correct places

	for (i=parptr[1,num_comp];i<=parptr[ncpar[comtyp[ncomp]],ncomp]-ncpar[comtyp[num_comp]];i=i+1) 
	     {
                par0[i] = par0[i+ncpar[comtyp[num_comp]]]
                blim[i] = blim[i+ncpar[comtyp[num_comp]]]
                tlim[i] = tlim[i+ncpar[comtyp[num_comp]]]  
                step[i] = step[i+ncpar[comtyp[num_comp]]] 
                ptol[i] = ptol[i+ncpar[comtyp[num_comp]]] 
                ifix[i] = ifix[i+ncpar[comtyp[num_comp]]]
		call strcpy(comments[1,i+ncpar[comtyp[num_comp]]],comments[1,i],SZ_LINE) 
              }


	ncomp = ncomp - 1


        for (i=num_comp;i<=ncomp ;i = i + 1) {
                comtyp[i] = comtyp[i+1]
          	inplot[i] = inplot[i+1]
		}
	comtyp[ncomp+1] = 0
	inplot[ncomp+1] = NO


	#Update the pointers the values
	k = 0
        for (i=1;i<=ncomp;i=i+1) {
                for (j=1;j<=ncpar[comtyp[i]];j=j+1) {
                        k = k + 1
                        parptr[j,i] = k
                        }
                }

	for (i=1;i<=ncomp;i=i+1) {
		for (j=1;j<=ncpar[comtyp[i]];j=j+1) {
			if (ifix[parptr[j,i]]!=0 && ifix[parptr[j,i]]!=-1 && i!=num_comp) {
				if (ifix[parptr[j,i]]>=num_comp) {
					ifix[parptr[j,i]]=ifix[parptr[j,i]]-1
				}
			}
		}
	}

	#Update fpar and get the new chisq
	call setpar(nfree,fpar)
	call chispec(nfree, fpar, chisq)

end








procedure sfcomponent_list()

include "specfit.com"

int i

begin

	for (i=1;i<=ncomp;i=i+1) {
		call printf("%d \t\t%s\n")
			call pargi(i)
			call pargstr(compkeys[1,comtyp[i]])
		}



end










