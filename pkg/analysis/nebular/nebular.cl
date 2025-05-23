#  NEBULAR package -- Nebular diagnostics routines
#
#	Last update: 9-Jul-97

procedure nebular()

string	mode		= "al"

begin
	package nebular

	task	abund,
		ionic,
		ntcontour,
		ntplot,
		redcorr,
		temden,
		zones		= "nebular$x_nebular.e"

	task	$at_data	= "nebular$at_data.cl"
	task	$nlevel		= "nebular$nlevel.cl"

	# Psets - Line flux column names:
	task	diagcols	= "nebular$diagcols.par"
	task	fluxcols	= "nebular$fluxcols.par"

	task	faluminum	= "nebular$faluminum.par"
	task	fargon		= "nebular$fargon.par"
	task	fcalcium	= "nebular$fcalcium.par"
	task	fcarbon		= "nebular$fcarbon.par"
	task	fchlorine	= "nebular$fchlorine.par"
	task	fmagnesium	= "nebular$fmagnesium.par"
	task	fneon		= "nebular$fneon.par"
	task	fnitrogen	= "nebular$fnitrogen.par"
	task	foxygen		= "nebular$foxygen.par"
	task	fpotassium	= "nebular$fpotassium.par"
	task	fsilicon	= "nebular$fsilicon.par"
	task	fsodium		= "nebular$fsodium.par"
	task	fsulfur		= "nebular$fsulfur.par"

	hidetask faluminum 
	hidetask fargon 
	hidetask fcalcium 
	hidetask fcarbon 
	hidetask fchlorine 
	hidetask fmagnesium
	hidetask fneon 
	hidetask fnitrogen 
	hidetask foxygen 
	hidetask fpotassium 
	hidetask fsilicon 
	hidetask fsodium 
	hidetask fsulfur

	cl()
end
