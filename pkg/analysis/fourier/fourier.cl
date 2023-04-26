# FOURIER.CL -- Script to set up tasks in the IRAF fourier package

procedure fourier ()

string	mode="al"
#------------------------------------------------------------------------------
begin

package fourier

task    crosscor,
	forward,
	taperedge,
	inverse	        = "fourier$x_fourier.e"
cl()
end
