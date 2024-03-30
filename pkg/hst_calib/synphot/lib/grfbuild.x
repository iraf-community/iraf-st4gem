include	"../grflist.h"

# GRF_BUILD -- Build the graph list structlure from the graph arrays
#
# B.Simon	12-Aug-88	First Code

procedure grf_build (ngraph, nkeyword, ncomp, innode, outnode, 
		     compname, keyword, inidx, gp)

int	ngraph		# i: Number of rows in graph
int	nkeyword	# i: Length of keyword string
int	ncomp		# i: Length of component name string
pointer	innode		# i: Array of innode values
pointer	outnode		# i: Array of outnode values
pointer	compname	# i: Array of component names
pointer	keyword		# i: Array of mnemonic keywords
int	inidx		# i: Index array sorted on innode and compname
pointer	gp		# o: Pointer to graph list structure
#--
int	temp, idx, jdx, kdx

begin
	# Create and initialize graph data structure

	gp = NULL
	call malloc (gp, LEN_GRFSTRUCT, TY_STRUCT)

	GRF_SIZE(gp) = ngraph
	GRF_CMPLEN(GP) = ncomp + 1
	GRF_KEYLEN(gp) = nkeyword + 1
	GRF_INARRAY(gp) = innode
	GRF_OUTARRAY(gp) = outnode
	GRF_CMPARRAY(gp) = compname
	GRF_KEYARRAY(gp) = keyword

	# Use the row at the end of the data structure for temporary storage

	temp = ngraph + 1

	# Loop over all rows, moving them into sorted order

	do idx = 1, ngraph {

	    # The inidx array forms one or more cycles. Move the first
	    # row in the cycle to the temporary location. Repeatedly
	    # move the remaining rows in the cycle until the final
	    # location of the first row is found. Move the first row
	    # from its temporary location to its final location. Update
	    # the inidx array to indicate which rows have been moved.

	    if (Memi[inidx+idx-1] != idx) {
		call grf_copy (gp, idx, temp)

		jdx = idx
		while (Memi[inidx+jdx-1] != idx) {
		    kdx = Memi[inidx+jdx-1]
		    call grf_copy (gp, kdx, jdx)
		    Memi[inidx+jdx-1] = jdx
		    jdx = kdx
		}

		call grf_copy (gp, temp, jdx)
		Memi[inidx+jdx-1] = jdx
	    }

	}

end

# GRF_COPY -- Copy one row in a graph list to another

procedure grf_copy (gp, idx, jdx)

pointer	gp		# i: Pointer to graph list structure
int	idx		# i: Row to copy from
int	jdx		# i: Row to copy to
#--

begin
	GRF_INNODE(gp,jdx) = GRF_INNODE(gp,idx)
	GRF_OUTNODE(gp,jdx) = GRF_OUTNODE(gp,idx)

	call strcpy (GRF_COMPNAME(gp,idx), GRF_COMPNAME(gp,jdx),
		     GRF_CMPLEN(gp)-1)
	call strcpy (GRF_KEYNAME(gp,idx), GRF_KEYNAME(gp,jdx),
		     GRF_KEYLEN(gp)-1)

end
