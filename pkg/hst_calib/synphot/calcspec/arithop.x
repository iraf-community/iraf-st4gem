# ARITHOP -- Perform arithmetic operations 

procedure arithop(word, pending, form1, form, nwv, wv, sp1, sp)

char 	word[ARB]		# i:  operation parsed
char	pending[ARB]		# io: Pending arithmetic operation
char 	form1[ARB]		# io: Form (units) of operand
char	form[ARB]		# io: Form (units) of other operand
int	nwv			# i: number of wavelengths
real	wv[ARB]			# i: Array of wavelengths
real	sp1[ARB]		# i: array of spectral values
real 	sp[ARB]			# io: array of spectral values

int 	i
int 	strsearch(), strlen()
bool	status

begin

	# If there is a pending operation, change spectra to same form of flux
	if( strlen( pending ) > 0 ) {
	   if( strsearch(form1, "mag") > 0 ) {
	      call specform( nwv, wv, sp1, form1, sp1, "flam", status )
	      call strcpy( "flam", form1, SZ_LINE)
	   }
	
	   if( strsearch(form,"mag") > 0 ) {
	      call specform( nwv, wv, sp, form, sp, form1, status )
	      call strcpy(form1,form,SZ_LINE)
	   }

	   # Change spectra to the same flux units
	   if( strsearch(form,form1) < 0 ) 
	      call specform( nwv, wv, sp1, form1, sp1, form, status )

	   # Perform the pending spectrum arithmetic
	   if( strsearch(pending,"+") > 0 )
              do i=1,nwv
	         sp[i] = sp[i] + sp1[i]
	   else if( strsearch(pending,"-") > 0 )
	      do i=1,nwv
	         sp[i] = sp1[i] - sp[i]

	   call strcpy("",pending,SZ_LINE)		# Op no longer pending
	}
	
	# Setup new pending operation and stow current spectrum (if any) 
	call strcpy( word, pending, SZ_LINE)
        if( strsearch(pending,"+") > 0 || strsearch(pending,"-") > 0 ) {
	   do i=1,nwv
	      sp1[i] = sp[i]
	   call strcpy(form,form1,SZ_LINE)
        } else
	   call strcpy("",pending,SZ_LINE)
	
end
