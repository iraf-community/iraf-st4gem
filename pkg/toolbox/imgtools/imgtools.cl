procedure imgtools()
string	mode="al"

begin
	package imgtools

	task addmasks,
	     imcalc,
	     iminsert,
	     improject,
	     rd2xy,
	     xy2rd,
	     xyztable,
	     xyztoim              = "imgtools$x_imgtools.e"

	cl()
end
