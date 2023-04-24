procedure fitting()
string	mode="al"

begin
	package fitting
	task gfit1d  = "fitting$x_fitting.e"
	task samplepars  = "fitting$samplepars.par"
	task errorpars   = "fitting$errorpars.par"
	cl()
end
