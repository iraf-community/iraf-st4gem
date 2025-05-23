# Copyright restrictions apply - see st4gem$copyright.stsdas

#-------------------------------------------------------------------29 Jan 96---
.help wt_avg.x Jan96 nebular/fivel
.ih
NAME
wt_avg - Calculate weighted average of N values in an array, excluding INDEFs. 
.endhelp
#-------------------------------------------------------------------------------
#  WT_AVG -	Calculates weighted average of first N_AVG values in an array, 
# 		excluding INDEFs. Note that WT array should not contain INDEFs.
#		Returned value is INDEFR if all values are INDEF.  
#
#  Revision history:
#	24-Jan-96 by RAShaw	Initial implementation

real procedure wt_avg (v, wt, n_avg)

#  Calling arguments:
real	v[ARB]	 	# I: values to be averaged
real	wt[ARB]		# I: weights of values
int	n_avg		# I: no. of values to average
real	avg		# O: weighted average of !INDEFR values

#  Local variables:
int	i		# generic
int	n_good		# no. non-INDEF values in array
real	sum_val		# sum of !INDEFR values
real	sum_wts		# sum of weights of !INDEFR values

begin
	sum_val = 0.
	sum_wts = 0.
	n_good  = 0
	do i = 1, n_avg {
	    if ( !IS_INDEFR(v[i]) ) {
	    	n_good = n_good + 1
		sum_val = sum_val + v[i] * wt[i]
		sum_wts = sum_wts + wt[i]
	    }
	}

	if (n_good > 0)
	    avg = sum_val / sum_wts
	else 
	    avg = INDEFR

	return (avg)
end


