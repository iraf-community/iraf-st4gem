# MERGEWAVE -- Merges two wavelength sets
#
#input:
#	n1	= number of wavelengths of wavelength set 1
#	w1(n1)	= wavelength set 1
#	n2	= number of wavelengths of wavelength set 2
#	w2(n2)	= wavelength set 2
#	maxout	= max number of output wavelengths
# output
#	maxout	= number of output wavelengths
#	wout(maxout)	= output merged wavelength set
# nov 1987 Keith Horne @ stsci
# jul 1989 spp version by Dave Bazell

procedure mergewave( n1, w1, n2, w2, maxout, wout )

int	n1		# Number of wavelengths in w1
int	n2		# Number of wavelengths in w2
int	maxout		# Input: dim of wout;  Output: Number of wavel. in wout
int	i1		# index for w1
int	i2		# index for w2
int	nout		# index for wout 
real	w1[ARB]		# First array of input wavelengths
real	w2[ARB]		# Second array of input wavelengths
real	wout[ARB]	# Output array of merged wavelengths

begin

# trap invalid inputs
	if( n1 <= 0 || n2 <= 0 || maxout <= 0 )
	{
	   call printf("** Invalid wavelength set dimensions.\n")
	   call printf("** n1=%d, n2=%d, nout=%d\n")
	      call pargi(n1)
	      call pargi(n2)
	      call pargi(nout)
	   call printf("** Mergewave Aborted\n")
	   maxout = 0
	   return
	}

	i1 = 1			#index for array w1
	i2 = 1			#index for array w2
	nout = 1		#index for array wout

	while ( nout <= maxout)
	{

# select the smallest of the two possible wavelengths
      	   if( i1 <= n1 )
	   {
      	      if( i2 <= n2 )
	      {
      	         if( w1[i1] < w2[i2] )
	         {
	            wout[nout] = w1[i1]
	            i1 = i1 + 1
	         }
     	         else if( w1[i1] == w2[i2] )
	         {
	            wout[nout] = w1[i1]
	            i1 = i1 + 1
	            i2 = i2 + 1
	         }
                 else if( w1[i1] > w2[i2] )
	         {
	            wout[nout] = w2[i2]
	            i2 = i2 + 1
	         }
	      }
# load wavelengths from first set when second set is finished
              else
	      {
                 wout[nout] = w1[i1]
	         i1 = i1 + 1
	      }
	   }

# load wavelengths from second set if first set finishes first
           else if( i2 <= n2 )
	   {
	      wout[nout] = w2[i2]
	      i2 = i2 + 1
	   }

# both input wavesets are finished so return normally
           else
	   {
	      maxout = nout - 1
	      return
           }

# check for invalid wavelength point and delete if bad
	   if( wout[nout] <= 0. ) nout = nout - 1

# increment output counter and go back for next point
	   nout = nout + 1
	}

# Output array is full but input arrays are not depleted: signal truncation 
# warning and return
        if ( i1 <= n1 || i2 <= n2 )
	{
	  call printf("** Mergewave warning: output wavelength set overflow\n") 
	  call printf("** at nout= %d, wout= %d \n")
	     call pargi(nout)
	     call pargi(wout[nout])
	  maxout = nout
	  return
	}
end
