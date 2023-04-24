include "../fourier.h"

# ft_fft -- fast Fourier transform
# The input and output are pairs of arrays containing the real & imaginary
# parts.  The input pair are packed into a complex array, the Fourier
# transform is taken (either forward or inverse), and the result is unpacked
# into the output pair.
#
# C. D. Biemesderfer, STScI, 2 Dec 87
# Phil Hodge, 18-Aug-1988  Just do the transform.

procedure ft_fft (work, inreal, inimag, outreal, outimag, fwd)

pointer work		# i: structure for work space
real	inreal[ARB]	# i: input data (real part)
real	inimag[ARB]	# i: input data (imaginary part)
real	outreal[ARB]	# o: output data (real part)
real	outimag[ARB]	# o: output data (imaginary part)
bool	fwd		# i: forward transform?
#--
int	npts

begin
	npts = NPTS(work)

	if (fwd) {

	    # Pack input data into complex.
	    call apkxr (inreal, inimag, Memx[XWORK(work)], npts)

	    # Do complex forward transform.
	    call cfftf (npts, Memx[XWORK(work)], Memr[TRIGTAB(work)])

	    # Unpack complex results.
	    call aupxr (Memx[XWORK(work)], outreal, outimag, npts)

	} else {

	    # Pack input data into complex.
	    call apkxr (inreal, inimag, Memx[XWORK(work)], npts)

	    # Do complex inverse transform.
	    call cfftb (npts, Memx[XWORK(work)], Memr[TRIGTAB(work)])

	    # Unpack complex results.
	    call aupxr (Memx[XWORK(work)], outreal, outimag, npts)
	}
end
