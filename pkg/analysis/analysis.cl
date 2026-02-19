# ANALYSIS.CL -- The STScI data analysis suite of packages

procedure analysis()
string	mode="al"

begin
	set dither     = "analysis$dither/"
	set fitting    = "analysis$fitting/"
	set fourier    = "analysis$fourier/"
	set isophote   = "analysis$isophote/"
        set nebular    = "analysis$nebular/"
        set at_data    = "nebular$atomic_data/"

	package analysis

	task dither.pkg 	= "dither$dither.cl"
	task fitting.pkg 	= "fitting$fitting.cl"
	task fourier.pkg 	= "fourier$fourier.cl"
	task isophote.pkg 	= "isophote$isophote.cl"
        task nebular.pkg        = "nebular$nebular.cl"

	cl()
end
