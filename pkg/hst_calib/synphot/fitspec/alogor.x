#* HISTORY *
#* B. Simon	29-Jul-94	original

# ALOGOR -- Logical or of two boolean arrays

procedure alogor (a, b, c, n)

bool	a[ARB]		# i: first array
bool	b[ARB]		# i: second array
bool	c[ARB]		# o: output array
int	n		# i: length of arrays
#--
int	i

begin
	do i = 1, n
	    c[i] = a[i] || b[i]
end
