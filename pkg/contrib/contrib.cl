procedure contrib()
string	mode="al"

begin
	set spfitpkg     = "contrib$spfitpkg/"

	package contrib

	task spfitpkg.pkg	= "spfitpkg$spfitpkg.cl"

	# Write the welcome message
	type contrib$contrib.msg

	clbye
end
