define	MIN_REAL 1.0e-37

# TOINDEF -- Set numbers less than MIN_REAL and negative numbers to INDEF

procedure toindef( data, ndat )

real	data[ARB]	#io: Data array
int	ndat		#i : number of data points
#--

int	ic

begin

	do ic = 1, ndat {
	   if ( data[ic] <= MIN_REAL )
	      data[ic] = INDEFR
	}
end
