#* HISTORY *
#* B. Simon	05-Aug-94	original

# VARNAME -- Construct a fit variable name

procedure varname (ivar, name, maxch)

int	ivar		# i: variable number
char	name[ARB]	# o: variable name
int	maxch		# i: maximum length of name
#--
string	badnumber  "varname: illgeal index to fit variable"
string	vnames     "vone,vtwo,vthree,vfour,vfive,vsix,vseven,veight,vnine"

int	word_find()

begin
	if (word_find (ivar, vnames, name, maxch) <= 0)
	    call printerr_int (badnumber, ivar)

end
