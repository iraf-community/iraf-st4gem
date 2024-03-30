define	MAXCH 32

# REPORTFIT -- Report the bestfit files

procedure reportfit( bchi2, bebmv, bscale, bmodel, speclist, nfree, nkeep)

real	bchi2[ARB]		# i: Chisqrs of best fits
real	bebmv[ARB]		# i: Extinctions of best fits
real	bscale[ARB]		# i: Scale factors of best fits
int	bmodel[ARB]		# i: Indices of models
char	speclist[SZ_LINE,ARB]	# i: Fit file names
int	nfree			# i: 
int	nkeep			# i: Number of bestfit spectra

int	ic
real	chi2rel
char	mult[MAXCH], ebmv[MAXCH], specname[SZ_FNAME]

string	ebfmt	" * ebmv(%g) "
string	mulfmt	" * %g "

begin
	call printf("\n")
	call printf("The best-fit models are as follows: N= %d\n")
	   call pargi( nfree )
	call printf("   rank    X^2/N      X^2/best\n")

	do ic = 1, nkeep {

	   # Copy current spectrum name to print string
	   call strcpy( speclist[1,bmodel[ic]], specname, SZ_LINE )

	   # If E(B-V) != 0 append extinction to print string
	   if ( bebmv[ic] != 0. ) {
	      call sprintf( ebmv, MAXCH, ebfmt )
	         call pargr( bebmv[ic] )
	      call strcat( ebmv, specname, SZ_LINE)
	   }

	   # Append scale factor to print string
	   call sprintf( mult, MAXCH, mulfmt )
	      call pargr( bscale[ic] )
	   call strcat( mult, specname, SZ_LINE )

	   # Calculate relative chisqr
	   if ( bchi2[ic] > 0. )
	      chi2rel = bchi2[ic] / bchi2[1]

	   call printf("%5d%11.4g%11.4g   %s\n")
	      call pargi( ic )
	      call pargr( bchi2[ic] )
	      call pargr( chi2rel )
	      call pargstr( specname )
	}
end
