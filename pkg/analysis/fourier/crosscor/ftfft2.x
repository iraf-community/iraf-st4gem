include "../fourier.h"

# This file contains ft_fft_2 and ft_xcorr.
#
# Phil Hodge,  3-Aug-1988  Subroutines created.
# Phil Hodge, 20-Jul-1993  In ft_xcorr, take complex conjugate of the
#			first array instead of the second; delete ft_conv.

# ft_fft_2 -- Fourier transform of two input arrays
# This routine takes two input sets (real & imaginary parts of each set)
# and packs them into a pair of complex arrays.  The complex arrays are
# longer than the input sets and are padded on the right with zero to
# a size of NPTS(work).  The Fourier transform is taken for both of the
# complex arrays, and the results are unpacked into the output arrays.
#
# If either npts1 or npts2 is zero, the corresponding output array
# is simply filled with zero.

procedure ft_fft_2 (work, npts1, npts2,
		inreal1, inimag1, inreal2, inimag2,
		outreal1, outimag1, outreal2, outimag2)

pointer work			# i: structure for work space
int	npts1			# i: length of first input array pair
int	npts2			# i: length of second input array pair
real	inreal1[ARB]		# i: Input data (real part)
real	inimag1[ARB]		# i: Input data (imaginary part)
real	inreal2[ARB]		# i: Input data (real part)
real	inimag2[ARB]		# i: Input data (imaginary part)
real	outreal1[ARB]		# o: transform of first array (real part)
real	outimag1[ARB]		# o: transform of first array (imaginary part)
real	outreal2[ARB]		# o: transform of second array (real part)
real	outimag2[ARB]		# o: transform of second array (imaginary part)
#--
int	allpts			# length of output array pair

begin
	allpts = NPTS(work)

	# Pack input data into complex, pad with zero, do complex forward
	# Fourier transform, and unpack into output arrays.

	if (npts1 > 0) {
	    call apkxr (inreal1, inimag1, Memx[XWORK(work)], npts1)
	    if (npts1 < allpts)
		call aclrx (Memx[XWORK(work)+npts1], allpts-npts1)
	    call cfftf (allpts, Memx[XWORK(work)], Memr[TRIGTAB(work)])
	    call aupxr (Memx[XWORK(work)], outreal1, outimag1, allpts)
	} else {
	    call aclrr (outreal1, allpts)
	    call aclrr (outimag1, allpts)
	}

	if (npts2 > 0) {
	    call apkxr (inreal2, inimag2, Memx[XWORK2(work)], npts2)
	    if (npts2 < allpts)
		call aclrx (Memx[XWORK2(work)+npts2], allpts-npts2)
	    call cfftf (allpts, Memx[XWORK2(work)], Memr[TRIGTAB(work)])
	    call aupxr (Memx[XWORK2(work)], outreal2, outimag2, allpts)
	} else {
	    call aclrr (outreal2, allpts)
	    call aclrr (outimag2, allpts)
	}
end

# ft_xcorr -- cross correlation of two input arrays
# This routine takes two input sets (real & imaginary parts of each set)
# and multiplies the complex conjugate of the first by the second,
# putting the result in the same location as the second set.

procedure ft_xcorr (real1, imag1, real2, imag2, allpts)

real	real1[ARB]		# i: first array (real part)
real	imag1[ARB]		# i: first array (imaginary part)
real	real2[ARB]		# io: second array & output (real part)
real	imag2[ARB]		# io: second array & output (imaginary part)
int	allpts			# i: size of arrays
#--
real	temp			# scratch for real part of product
int	k			# loop index

begin
	do k = 1, allpts {
	    temp     =  real1[k] * real2[k] + imag1[k] * imag2[k]
	    imag2[k] =  real1[k] * imag2[k] - imag1[k] * real2[k]
	    real2[k] =  temp
	}
end
