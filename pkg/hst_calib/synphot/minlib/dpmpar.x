include <mach.h>

#* HISTORY *
#* B.Simon	05-Aug-94	Original

# DPMPAR -- Return machine dependent constants to fit routine

double procedure dpmpar (what)

int	what		# i: which constant to return
#--
double	value

begin
	switch (what) {
        case 1:	# machine precision
            value = EPSILOND
        case 2:	# smallest number
            value = 1.0 / MAX_DOUBLE
        case 3:	# largest number
            value = MAX_DOUBLE
        }

        return (value)
end
