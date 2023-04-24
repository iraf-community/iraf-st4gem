# NL_POWER -- Raises array to a power.

procedure nl_power (x, z, ndata, gamma, amp, ref)

real	x[ARB]		# i: Indep. variable array.
real	z[ARB]		# o: Computed power-law function.
int	ndata		# i: Number of data points.
real	gamma		# i: Power-law index.
real	amp		# i: Power-law flux at reference x.
double	ref		# i: Reference x.

#--
int	i

begin
	do i = 1, ndata {
	    if (x[i] <= 0.) 
	        z[i] = 0.
	    else 
	        z[i] = amp * (x[i]/ref) ** gamma
	}
end
                                
