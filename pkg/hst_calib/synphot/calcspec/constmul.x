include	"mac.h"

# CONSTMUL -- Multiply a spectrum by a scalar constant

procedure constmul( script, form, nwv, wv, iw, sp )

char	script[ARB]	# i: Script containing commands and arguments
char	form[ARB]	# i: Form of spectrum i.e. units
int	nwv		# i: number of wavelengths
real	wv[ARB]		# i: array of wavelengths
int	iw		# io: index of last word parsed
real	sp[ARB]		# i: raw spectrum;  o: spectrum * constant

char	cword[SZ_FNAME]
int	nchar, ip
int	strsearch(), ctowrd(), ctor()
real 	constant, addmag
bool	status

begin	

	# Parse the next word in script to get value of constant
	nchar = ctowrd( script, iw, cword, SZ_FNAME)
	ip = 1
	nchar = ctor(cword,ip,constant)

	# If constant <= 0 and form is magnitudes change form out of magnitudes
	#   if not mags then multiply by const

        if( constant <= 0. ) {

	   if( strsearch(form,"stmag") > 0 ) {
	      call specform( nwv, wv, sp, form, sp, "flam", status )
	      call strcpy("flam",form,SIZE_STR)
	   }

	   else if( strsearch(form,"mag") > 0 ) {
	      call specform( nwv, wv, sp, form, sp, "fnu" , status )
	      call strcpy("fnu", form, SIZE_STR)
	   }
	   call bmulkr( sp, constant, sp, nwv )


	} else {

	   if( strsearch(form,"mag") > 0 ) {
	      addmag = -2.5 * alog10( constant )
	      call baddkr( sp, addmag, sp, nwv )
	   } else
	      call bmulkr( sp, constant, sp, nwv )
	}
end
