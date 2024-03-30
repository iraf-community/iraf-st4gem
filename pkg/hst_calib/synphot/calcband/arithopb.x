# ARITHOPB -- Perform arithmetic operations on passbands

procedure arithopb(word, pending, nwv, wv, sp1, sp)

char 	word[ARB]		# i:  operation parsed
char	pending[ARB]		# io: Pending arithmetic operation
int	nwv			# i: number of wavelengths
real	wv[ARB]			# i: Array of wavelengths
real	sp1[ARB]		# i: array of spectral values
real 	sp[ARB]			# io: array of spectral values

int 	i
int 	strsearch(), strlen()

begin

	# If there is a pending operation
	if( strlen( pending ) > 0 ) {

	   # Perform the pending spectrum arithmetic
	   if( strsearch(pending,"*") > 0 )
	      call bmulr( sp, sp1, sp, nwv )

	   else if( strsearch(pending,"/") > 0 )
	      do i=1,nwv {
	         if (sp[i] != 0. && !IS_INDEFR (sp[i]) && !IS_INDEFR (sp[i]))
	            sp[i] = sp1[i]/sp[i]
	         else
	            sp[i] = INDEFR
	      }

	   call strcpy("",pending,SZ_LINE)		# Op no longer pending
	}
	
	# Setup new pending operation and stow current spectrum (if any) 
	call strcpy( word, pending, SZ_LINE)
        if( strsearch(pending,"*") > 0 || strsearch(pending,"/") > 0 ) {
	   do i=1,nwv
	      sp1[i] = sp[i]
        } else
	   call strcpy("",pending,SZ_LINE)
	
end
