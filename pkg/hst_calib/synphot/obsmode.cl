# OBSMODE.CL -- List information on the observation mode string

procedure obsmode ()

string	mode = "a"

begin
	page ("synphot$doc/obsmode.hlp",first=1,map+,clear+)
end
