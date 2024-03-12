# SPFITPKG.CL - Package script task for SPFITPKG.
# 16-Aug-95 by GAKriss, Johns Hopkins Univeristy

procedure spfitpkg()

string	mode		= "al"

begin
	set specfit="spfitpkg$specfit/"
	package spfitpkg

	task	specfit,
		dbcreate,
		dbcheck		= "spfitpkg$x_spfitpkg.e"

	cl()
end
