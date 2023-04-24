#} GRAPHICS.CL - Package script for GRAPHICS Package
procedure graphics()
string	mode="al"

begin
	set stplot = "graphics$stplot/"

	package graphics
	task stplot.pkg = "stplot$stplot.cl"
	clbye()
end
