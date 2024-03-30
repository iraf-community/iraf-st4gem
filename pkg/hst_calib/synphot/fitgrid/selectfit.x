include	"../plspec/plspec.h"

# SELECTFIT -- Select the nth ranked fit from among a list of ranked fits

procedure selectfit( speclist, ebmvs, models, rank, bestmod)

char	speclist[SZ_LINE,MAXSPEC]	# i: List of spectra fit
real	ebmvs[ARB]			# i: List of extinctions
int	models[ARB]			# i: List of model indices
int	rank				# i: Rank of model to select
char	bestmod[ARB]			# o: selected model

int	ic

begin

	# Find the index of the appropriately ranked model
	ic = models[rank]

	# Select the desired model
	call strcpy( speclist[1,ic], bestmod, SZ_LINE )

	# Add extinction if value is non-zero
	if ( ebmvs[ic] != 0. )
	   call applyebmv( bestmod, bestmod, ebmvs[ic], SZ_LINE )
	   
end
