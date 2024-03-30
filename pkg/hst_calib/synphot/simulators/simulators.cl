procedure simulators()
string	mode="al"

begin

	package simulators


	task	simbackgd,
		simimg,
		simnoise,
		simspec	= "simulators$x_simulators.e"

	# pset definitions

	task	refdata	 =	"simulators$refdata.par"
	task	simmodp  =	"simulators$simmodp.par"
	task	simbackp =	"simulators$simbackp.par"
	task	simcatp  =	"simulators$simcatp.par"

	# Write the banner message

	if (access ("simulators$simulators.msg")) {
		type ("simulators$simulators.msg")
	}

	clbye()
end
