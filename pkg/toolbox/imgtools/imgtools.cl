procedure imgtools()
string	mode="al"

begin
        set  mstools = "imgtools$mstools/"
        set gcombine = "imgtools$gcombine/"

	package imgtools

	task addmasks,
	     boxinterp,
	     countfiles,
	     gcopy,
	     gstatistics,
	     imcalc,
             imfill,
	     iminsert,
	     improject,
	     listarea,
	     mkgauss,
	     moveheader,
	     pickfile,
	     pixedit,
	     pixlocate,
	     rbinary,
	     rd2xy,
	     stack,
	     xy2rd,
	     xyztable,
	     xyztoim    = "imgtools$x_imgtools.e"

        task    gcombine    = "gcombine$x_gcombine.e"
        task 	mstools.pkg = "mstools$mstools.cl"


# Psets for gstat:
	task	gstpar = "imgtools$gstpar.par"

	cl()

end
