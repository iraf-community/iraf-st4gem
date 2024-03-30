procedure synphot()

string	version = "6 June 2002"
string	mode="al"

begin
	package synphot


	task	bandpar,
		calcspec,
		calcband,
		calcphot,
		countrate,
		genwave,
		fitband,
		fitgrid,
		fitspec,
		grafcheck,
		graflist,
		grafpath,
		imspec,
		mkthru,
		obsmode,
		plspec,
		pltrans,
		plratio,
		plband,
		showfiles,
                thermback = "synphot$x_synphot.e"


	# pset for all tasks

	task	refdata	 =	"synphot$refdata.par"

	# subpackage for simulators

        set     simulators = "synphot$simulators/"
        task    simulators.pkg = "simulators$simulators.cl"

	clbye()
end
