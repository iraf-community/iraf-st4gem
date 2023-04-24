# TOOLBOX.CL -- General Tools Package
#

procedure toolbox()
string	mode="al"

begin
	set headers     = "toolbox$headers/"
	set imgtools    = "toolbox$imgtools/"
	set tools       = "toolbox$tools/"

	package toolbox

	task headers.pkg = "headers$headers.cl"
	task imgtools.pkg = "imgtools$imgtools.cl"
	task tools.pkg = "tools$tools.cl"

	cl()
end
