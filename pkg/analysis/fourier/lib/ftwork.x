include "../fourier.h"

# ft_work_init -- create scratch; initialize
# Create scratch arrays for NCAR routines (complex array to contain
# the input data; work space including a table of sines & cosines).
# If phase shifting according to coordinate info was requested,
# create scratch for cosine & sine table.  If output is to be centered,
# create an array that is half as long as the data array.
#
# Phil Hodge, 15-Jul-1988  Subroutine created.

procedure ft_work_init (work, coord_shift, center)

pointer work		# i: pointer to scratch and flags
bool	coord_shift	# i: create space for phase shifting?
bool	center		# i: create space for centering output?
#--
int	npts		# size of array to be transformed
int	hnpts		# npts / 2 + 1
int	n_shift		# npts - shift, for centering (or decentering)
int	wsiz

begin
	npts = NPTS(work)
	hnpts = npts / 2 + 1
	n_shift = npts - abs (SHIFT(work))

	# Allocate (one or two) complex-type arrays for NCAR routine.
	call malloc (XWORK(work), npts, TY_COMPLEX)
	if (N_FILES(work) > 1)
	    call malloc (XWORK2(work), npts, TY_COMPLEX)

	# Allocate scratch space and initialize NCAR trig table.
	wsiz = 4*npts + 15
	call calloc (TRIGTAB(work), wsiz, TY_REAL)
	call cffti (npts, Memr[TRIGTAB(work)])

	# Create scratch for phase shifting according to coordinate info.
	if (coord_shift) {
	    call malloc (COSTAB(work), hnpts, TY_REAL)
	    call malloc (SINTAB(work), hnpts, TY_REAL)
	}

	# Allocate space for centering output.
	if (center)
	    call malloc (C_COPY(work), n_shift, TY_REAL)
end


# ft_work_free -- deallocate memory for NCAR routines and shifting

procedure ft_work_free (work, coord_shift, center)

pointer work		# i: pointer to scratch and flags
bool	coord_shift	# i: deallocate space used for phase shifting?
bool	center		# i: deallocate space used for centering output?
#--

begin
	# Deallocate NCAR scratch space.
	call mfree (TRIGTAB(work), TY_REAL)
	call mfree (XWORK(work), TY_COMPLEX)
	if (N_FILES(work) > 1)
	    call mfree (XWORK2(work), TY_COMPLEX)

	# Deallocate space for phase shift.
	if (coord_shift) {
	    call mfree (COSTAB(work), TY_REAL)
	    call mfree (SINTAB(work), TY_REAL)
	}

	# Deallocate space for centering output.
	if (center)
	    call mfree (C_COPY(work), TY_REAL)
end
