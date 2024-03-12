###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	procedure sfdogrid(nfree, fpar, chisq)
#
#  Description:	SFDOGRID allows one to fit a grid of values
#		The remaining values are fit with a fitting algorithm
#		and returns chisquare and the updated parameters
#		When complete SFDOGRID puts the best values of the parameters
#		into the database and outputs the whole grid to the log files.
#
#  Arguments:	int	nfree		- Number of free params
#		real	fpar[ARB]	- Array containing the free params
#
#  Returns:	real	chisq		- Chisquare for the current fit
#
#  Notes:	Shares data in "specfit.com"
#
#  History:	Dec. 1991   Knox Long
#		6/95		grimes	modified to meld it with specfit
#
###########################################################################

include	"specfit.h"

procedure sfdogrid(nfree, fpar, chisq,fit_type)
int	nfree
real	fpar[ARB]
real	chisq
int	fit_type


real	chisq_best

extern chispec

include	"specfit.com"

int	grid_num	#No of parameterrs to be gridded
int	grid_comp[2],grid_par[2]  #component no, parameter number of that component
real	grid_par0[2], grid_blim[2],grid_tlim[2],grid_step[2]
real	grid_par0_best[2]
int	ix,iy

int	clgeti()
real	clgetr()

begin
# call printf("Entered sfdogrid.\n")

chisq_best=1.e21

grid_num=clgeti("grid_num")
ix = 0

if(grid_num==1){  #Then it's only a one dimensional grid
	grid_comp[1]=clgeti("grid_comp1")
	grid_par[1]=clgeti("grid_par1")
	call printf("Component %d Parameter %d Current Value is %g\n")
		call pargi(grid_comp[1])
		call pargi(grid_par[1])
		call pargr(par0[ parptr[grid_par[1],grid_comp[1]]])
	grid_blim[1]=clgetr("grid_blim1")
	grid_tlim[1]=clgetr("grid_tlim1")
	grid_step[1]=clgetr("grid_step1")

	call printf("Your grid will range from %g to %g with stepsize %g.\n")
		call pargr(grid_blim[1])
		call pargr(grid_tlim[1])
		call pargr(grid_step[1])

# Now do the grid

	call log_grid_init(grid_num,grid_comp,grid_par,grid_blim,grid_tlim,grid_step)

	for(grid_par0[1]=grid_blim[1];grid_par0[1]<=grid_tlim[1];grid_par0[1]=grid_par0[1]+grid_step[1]){
		par0[ parptr[grid_par[1],grid_comp[1]]]=grid_par0[1]
		
		call sfdofit(nfree,fpar,chisq,fit_type)
		#call sfdosimplex(nfree,fpar,chisq)
		
		call printf("For Parameter Value %g Best Chisq is %g.\n\n")
			call pargr(grid_par0[1])
			call pargr(chisq)
		ix = ix + 1
		call log_grid(grid_num,grid_par0,chisq,ix,0)
		if(chisq < chisq_best){
			chisq_best=chisq
			grid_par0_best[1] = grid_par0[1]
			}
		}

	par0[ parptr[grid_par[1],grid_comp[1]]]=grid_par0_best[1] #Redo best fit
	call sfdofit(nfree,fpar,chisq,fit_type)
	#call sfdosimplex(nfree,fpar,chisq)
	}

else {  #Then it's a two dimensional grid

	grid_comp[1]=clgeti("grid_comp1")
	grid_par[1]=clgeti("grid_par1")
	call printf("1:  Component %d Parameter %d Current Value is %g\n")
		call pargi(grid_comp[1])
		call pargi(grid_par[1])
		call pargr(par0[ parptr[grid_par[1],grid_comp[1]]])
	grid_blim[1]=clgetr("grid_blim1")
	grid_tlim[1]=clgetr("grid_tlim1")
	grid_step[1]=clgetr("grid_step1")
	call printf("1:  Your grid will range from %g to %g with stepsize %g.\n")
		call pargr(grid_blim[1])
		call pargr(grid_tlim[1])
		call pargr(grid_step[1])

	#second Parameter Values
	grid_comp[2]=clgeti("grid_comp2")
	grid_par[2]=clgeti("grid_par2")
	call printf("2:  Component %d Parameter %d Current Value is %g\n")
		call pargi(grid_comp[2])
		call pargi(grid_par[2])
		call pargr(par0[ parptr[grid_par[2],grid_comp[2]]])
	grid_blim[2]=clgetr("grid_blim2")
	grid_tlim[2]=clgetr("grid_tlim2")
	grid_step[2]=clgetr("grid_step2")
	call printf("2:  Your grid will range from %g to %g with stepsize %g.\n")
		call pargr(grid_blim[2])
		call pargr(grid_tlim[2])
		call pargr(grid_step[2])

# Now do the grid

	call log_grid_init(grid_num,grid_comp,grid_par,grid_blim,grid_tlim,grid_step)

	for(grid_par0[1]=grid_blim[1];grid_par0[1]<=grid_tlim[1];
				grid_par0[1]=grid_par0[1]+grid_step[1]){

		ix = ix + 1
		iy = 0

		par0[ parptr[grid_par[1],grid_comp[1]]]=grid_par0[1]
		for (grid_par0[2] = grid_blim[2];grid_par0[2] <= grid_tlim[2];
			grid_par0[2]=grid_par0[2]+grid_step[2]){

			iy = iy + 1

			par0[ parptr[grid_par[2],grid_comp[2]]]=grid_par0[2]
			call sfdofit(nfree,fpar,chisq,fit_type)
			#call sfdosimplex(nfree,fpar,chisq)

			call printf("For Parameter Values 1: %g and 2: %g Best Chisq is %g.\n\n")
				call pargr(grid_par0[1])
				call pargr(grid_par0[2])
				call pargr(chisq)

			call log_grid(grid_num,grid_par0,chisq,ix,iy)
			if (chisq < chisq_best) {
				chisq_best=chisq
				grid_par0_best[1]=grid_par0[1]
				grid_par0_best[2]=grid_par0[2]
			}
		}
	}

	call log_grid_bottom()
	par0[ parptr[grid_par[1],grid_comp[1]]]=grid_par0_best[1] #Redo best fit
	par0[ parptr[grid_par[2],grid_comp[2]]]=grid_par0_best[2] #Redo best fit
	call sfdofit(nfree,fpar,chisq,fit_type)
	#call sfdosimplex(nfree,fpar,chisq)

}

#call printf("Left sfdogrid. \n")
end


procedure log_grid(grid_num,grid_par0,chisq,ix,iy)

int	grid_num
real 	grid_par0[2],chisq
int 	ix, iy
int	j

include "specfit.com"

begin
for (j = 1; j <= nlogfd; j = j + 1) {
	
	if (grid_num == 1) {
		call fprintf(logfd[j],"   %3d\t\t\t\t%14.8g\t\t%15.7g\n")
			call pargi(ix)
			call pargr(grid_par0[1])
			call pargr(chisq)
	} else {
		call fprintf(logfd[j],"   %3d %3d\t\t%14.8g %14.8g\t%15.7g\n")
			call pargi(ix)
			call pargi(iy)
			call pargr(grid_par0[1])
			call pargr(grid_par0[2])
			call pargr(chisq)
	}
}
end

procedure log_grid_bottom()
int	j
include "specfit.com"

begin
	for (j = 1; j <= nlogfd; j = j + 1) {	
		call fprintf(logfd[j],"End Grid Values\n\n\n")
	}
end



procedure log_grid_init(grid_num,grid_comp,grid_par,grid_blim,grid_tlim,grid_step)
int	grid_num			#No of parameters to be gridded
int	grid_comp[2],grid_par[2]	#component num parameter number of that component
real	grid_blim[2],grid_tlim[2],grid_step[2]

include	"specfit.com"
int i,j
begin
	for (j = 1; j <= nlogfd; j=j+1){
		for(i = 1; i <= grid_num; i = i + 1){
			call fprintf(logfd[j],"GRID_INIT Component %d Parameter %d Minimum %g Maximum %g StepSize %g\n")
				call pargi(grid_comp[i])
				call pargi(grid_par[i])
				call pargr(grid_blim[i])
				call pargr(grid_tlim[i])
				call pargr(grid_step[i])
		}
		call fprintf(logfd[j],"Grid Position\t\t\t   Parameter(s)\t\t\t Chisq\n")
	}
end



procedure sfdofit(nfree, fpar, chisq, fit_type)
int	nfree
real	fpar[ARB]
real	chisq
int	fit_type


include	"specfit.com"

begin
	if (fit_type == 1) {
		call sfdomarquadt(nfree,fpar,chisq)
	} else if (fit_type == 2) {
		call sfdosimplex(nfree,fpar,chisq)
	} else if (fit_type == 3) {
		call sfalternate(nfree,fpar,chisq)
	} else if (fit_type == 4) {
		call sfdonumrec(nfree,fpar,chisq)
	} else {
		call eprintf("ERROR: Check procedure sfdofit in sfdogrid.x\nUnknown fit type!\n")
	}
end
