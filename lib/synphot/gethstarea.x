include	<synphot.h>

# GET_HSTAREA -- Read hstarea from parameter file

procedure get_hstarea (area)

real	area		# u: HST primary mirror area
#--
real	hstarea
data	hstarea	 / HSTAREA /

begin
	area = hstarea

	entry put_hstarea (area)

	hstarea = area
end
