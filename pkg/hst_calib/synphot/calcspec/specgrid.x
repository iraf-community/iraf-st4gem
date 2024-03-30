include	"mac.h"
include "../plspec/plspec.h"
 
# SPECGRID -- Spectrum Grid

procedure specgrid( script, nwave, wave, form, iw, spec)

char	script[ARB]	# i: Command script
int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: array of wavelengths
char	form[ARB]	# i: form (units) of spectrum
int	iw		# io: Position in script
real	spec[ARB]	# o: array of interpolated spectral points

char	fname[SZ_FNAME], spname[SZ_LINE,MAXSPEC], cseq[SIZE_STR]
char	form2[SIZE_STR], name[SZ_FNAME]

int	ip, i, i1, i2, nchar
int	ctowrd(), ctor(), strsearch()
real	part, omp, seq
bool	status
pointer	spec2

begin

	nchar = ctowrd( script, iw, fname, SZ_FNAME)
	nchar = ctowrd( script, iw, cseq, SZ_FNAME)
	ip = 1
	nchar = ctor( cseq, ip, seq)	

	# Evaluate first spectrum
	i1 = - max( 1, int( seq ) )
	call strcpy("@", name, SZ_FNAME)
	call strcat(fname, name, SZ_FNAME)
	call loadlist( name, i1, spname, MAXSPEC)
	call strlwr( spname[1,1] )

	#  Check if file is ascii
	if ( strsearch( spname[1,1], ".dat" ) > 0 )
	   call ascspec( nwave, wave, spec, spname[1,1], form, status )
	# Otherwise assume a its a table
	else
	   call evalspec( nwave, wave, spec, spname[1,1], "FLUX", form, status)

	# Evaluate second spectrum
	i2 = int(seq) + 1
	i2 = - i2

	call loadlist(name, i2, spname, MAXSPEC)
	call strlwr( spname[1,1] )

	if( i2 != i1 ) {
	   call malloc( spec2, nwave, TY_REAL)

	   #  Check if file is ascii
	   if ( strsearch(spname[1,1],".dat") > 0 )
	      call ascspec( nwave, wave, Memr[spec2], spname[1,1], form2,
	                    status )
	   # Otherwise assume a its a table
	   else
	      call evalspec( nwave, wave, Memr[spec2], spname[1,1], 
	                  "FLUX", form2, status )
	   part = ( seq - i1 ) / float( i2 - i1 ) 
	   if( part == 1. ) {
	      do i=1,nwave
	         spec[i] = Memr[spec2+i-1]
	      call strcpy( form2, form, SZ_FNAME)
	   } else if( part != 0. ) {
	   
	      # Convert spectra to magnitudes before interpolation
	      if( strsearch( form, "mag") > 0 || strsearch(form,"MAG") > 0)
	         call specform( nwave, wave, Memr[spec2], form2, Memr[spec2], 
	                        form, status )
	      else if( strsearch(form2,"mag") > 0 || strsearch(form2,"MAG")>0){
	         call specform( nwave, wave, spec, form, spec, form2, status )
	         call strcpy(form2, form, SIZE_STR)

	      } else {
	         call specform( nwave, wave, spec, form, spec, "abmag", status )
	         call specform( nwave, wave, Memr[spec2], form2, Memr[spec2], 
	                        "abmag", status )
	         call strcpy("abmag", form, SIZE_STR)

   	      }

	      # Linear interpolation
	      omp = 1.-part
	      do i=1,nwave
	         if ( !IS_INDEFR (spec[i]) && Memr[spec2+i-1] != INDEFR )
	            spec[i] = spec[i] * omp + Memr[spec2+i-1] * part
	         else
	            spec[i] = INDEFR
	   }
	}
end
