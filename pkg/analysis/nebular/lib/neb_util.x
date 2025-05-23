# Copyright restrictions apply - see st4gem$copyright.stsdas

include	"../atom.h"
include	"../neberr.h"

#--------------------------------------------------------------------25 Aug 08--
.help neb_util.x Aug08 nebular
.ih
NAME
.nf
     get_atom_typ - Get atomic number from string containing name
      get_ion_num - Create the magic "ionum" and check for validity
      get_ion_str - Fetch string appropriate for a given diagnostic ratio
     get_atom_str - Generate atom/ion string. 
get_atomic_symbol - Return atomic symbol given an atomic number. 
   get_ion_symbol - Generate ion symbol given an ion number. 
.fi
.endhelp
#-------------------------------------------------------------------------------
#  GET_ATOM_TYP -- Interpret string to determine enumerated type of atom.  
#		Input string must be in lower case.  Unrecognized types 
#		reported as zero.  

int procedure get_atom_typ (atom_string)

#  Arguments:
char	atom_string[ARB]	# I: string containing atomic name 
int	atom_type		# O: enumerated unit type; zero if unknown

#  Function used:
bool	streq()			# strings equal?

begin
	call strlwr (atom_string)

	# Set atno based upon dictionary entry...
	if (streq (atom_string, "hydrogen") )
	    atom_type = HYDROGEN

	else if (streq (atom_string, "helium") )
	    atom_type = HELIUM

	else if (streq (atom_string, "carbon") )
	    atom_type = CARBON

	else if ( streq (atom_string, "nitrogen") )
	    atom_type = NITROGEN

	else if ( streq (atom_string, "oxygen") )
	    atom_type = OXYGEN

	else if ( streq (atom_string, "neon") )
	    atom_type = NEON

	else if ( streq (atom_string, "flourine") )
	    atom_type = FLOURINE

	else if ( streq (atom_string, "sodium") )
	    atom_type = SODIUM

	else if ( streq (atom_string, "magnesium") )
	    atom_type = MAGNESIUM

	else if ( streq (atom_string, "aluminum") )
	    atom_type = ALUMINUM

	else if ( streq (atom_string, "silicon") )
	    atom_type = SILICON

	else if ( streq (atom_string, "sulfur") )
	    atom_type = SULFUR

	else if ( streq (atom_string, "chlorine") )
	    atom_type = CHLORINE

	else if ( streq (atom_string, "argon") )
	    atom_type = ARGON

	else if ( streq (atom_string, "potassium") )
	    atom_type = POTASSIUM

	else if ( streq (atom_string, "calcium") )
	    atom_type = CALCIUM

	else if ( streq (atom_string, "iron") )
	    atom_type = IRON

	else if ( streq (atom_string, "nickel") )
	    atom_type = NICKEL

	else if ( streq (atom_string, "selenium") )
	    atom_type = SELENIUM

	else if ( streq (atom_string, "krypton") )
	    atom_type = KRYPTON

	# ...or give up.
	else
	    atom_type = UNKNOWN

	return (atom_type)

end


#-------------------------------------------------------------------------------
#  GET_ION_NUM -- Create the magic "ionum" and check for validity.  

int procedure get_ion_num (atomic_no, ion, temden_code)

#  Arguments:
int	atomic_no	# I: atomic number
int	ion		# I: ion ([S II] => 1, etc.)
int	temden_code	# I: |T_e|N_e => |1|2
int	ion_num		# O: enumerated unit type; zero if unknown

#  Functions used:
#bool	ion_valid()	# valid ion_num specified?

begin
	ion_num = atomic_no * 1000 + (ion) * 10 + temden_code

	# Ensure a valid ion/spectrum was specified for calculation. 
#	if (!ion_valid (atomic_no, ion)) {
#	    ion_num = 0
#	    call fl_err (ION_INVALID)
#	} 

	return (ion_num)
end


#-------------------------------------------------------------------------------
#  GET_ION_STR -- Fetch string appropriate for a given diagnostic ratio.  

int procedure get_ion_str (at, diag, ion_string, max_chars)

#  Calling arguments:
pointer	at			# I: atomic data structure
int	diag			# I: diagnostic type (density|temperature)
char	ion_string[ARB]		# O: string containing units description 
int	max_chars		# I: size of string

#  Local variable:
int	ionum			# magic ion code

#  Function used:
int	get_ion_num()		# create the magic "ionum" 
int	strlen()		# length of a string

include	"../at.h"

begin
	ionum = get_ion_num (AT_ATOM(at), AT_ION(at), diag)

	# Build ion string based upon atomic number, spectrum, and tem_den switch
	switch (ionum) {
	case 6002:
	    call sprintf (ion_string, max_chars, 
		"C I: I(9823+9849)/I(8728)")
	case 6011:
	    call sprintf (ion_string, max_chars, 
		"C II: I(2326.1)/I(2327.6)")
	case 6021:
	    call sprintf (ion_string, max_chars, 
		"C III: I(1907)/I(1909)")
	case 7001:
	    call sprintf (ion_string, max_chars, 
		"N I: I(5198)/I(5200)")
	case 7002:
	    call sprintf (ion_string, max_chars, 
		"N I: I(5198+5200)/I(10397+10407)")
	case 7012:
	    call sprintf (ion_string, max_chars, 
		"N II: I(6548+6583)/I(5755)")
	case 7021:
	    call sprintf (ion_string, max_chars, 
		"N III: I(1749.7)/I(1752.2)")
	case 7031:
	    call sprintf (ion_string, max_chars, 
		"N IV: I(1483)/I(1487)")
	case 8002:
	    call sprintf (ion_string, max_chars, 
		"O I: I(6300+6363)/I(5577)")
	case 8011:
	    call sprintf (ion_string, max_chars, 
		"O II: I(3726)/I(3729)")
	case 8012:
	    call sprintf (ion_string, max_chars, 
		"O II: I(3726+3729)/I(7320+7330)")
	case 8022:
	    call sprintf (ion_string, max_chars, 
		"O III: I(4959+5007)/I(4363)")
	case 8031:
	    call sprintf (ion_string, max_chars, 
		"O IV: I(1401.2)/I(1404.8)")
	case 8041:
	    call sprintf (ion_string, max_chars, 
		"O V: I(1214)/I(1218)")
	case 10011:
	    call sprintf (ion_string, max_chars, 
		"Ne II: I(128216)")
	case 10022:
	    call sprintf (ion_string, max_chars, 
		"Ne III: I(3869+3969)/I(3342)")
	case 10031:
	    call sprintf (ion_string, max_chars, 
		"Ne IV: I(2423)/I(2425)")
	case 10032:
	    call sprintf (ion_string, max_chars, 
		"Ne IV: I(2422.4+2425.0)/I(1601.6+1601.4)")
	case 10042:
	    call sprintf (ion_string, max_chars, 
		"Ne V: I(3426+3346)/I(2975)")
	case 11032:
	    call sprintf (ion_string, max_chars, 
		"Na IV: I(3242+3362)/I(2805)")
	case 11052:
	    call sprintf (ion_string, max_chars, 
		"Na VI: I(2871+2970)/I(2569)")
	case 12042:
	    call sprintf (ion_string, max_chars, 
		"Mg V: I(2783+2928)/I(2418)")
	case 12062:
	    call sprintf (ion_string, max_chars, 
		"Mg VII: I(2506+2626)/I(2262)")
	case 13011:
	    call sprintf (ion_string, max_chars, 
		"Al II: I(2661)/I(2670)")
	case 13012:
	    call sprintf (ion_string, max_chars, 
		"Al II: I(2661+2670)/I(1671)")
	case 14011:
	    call sprintf (ion_string, max_chars, 
		"Si II: I(2335)/I(2345)")
	case 14021:
	    call sprintf (ion_string, max_chars, 
		"Si III: I(1883)/I(1892)")
	case 14022:
	    call sprintf (ion_string, max_chars, 
		"Si III: I(1883+1892)/I(1206)")
	case 16011:
	    call sprintf (ion_string, max_chars, 
		"S II: I(6716)/I(6731)")
	case 16012:
	    call sprintf (ion_string, max_chars, 
		"S II: I(6725)/I(4072)")
	case 16022:
	    call sprintf (ion_string, max_chars, 
		"S III: I(9069+9532)/I(6312)")
	case 16031:
	    call sprintf (ion_string, max_chars, 
		"S IV: I(1406)/I(1417)")
	case 17021:
	    call sprintf (ion_string, max_chars, 
		"Cl III: I(5517)/I(5537)")
	case 17022:
	    call sprintf (ion_string, max_chars, 
		"Cl III: I(5517+5537)/I(3353+3343)")
	case 17032:
	    call sprintf (ion_string, max_chars, 
		"Cl IV: I(7530+8045)/I(5323)")
	case 18011:
	    call sprintf (ion_string, max_chars, 
		"Ar II: I(69890)")
	case 18022:
	    call sprintf (ion_string, max_chars, 
		"Ar III: I(7136+7751)/I(5192)")
	case 18031:
	    call sprintf (ion_string, max_chars, 
		"Ar IV: I(4711)/I(4740)")
	case 18032:
	    call sprintf (ion_string, max_chars, 
		"Ar IV: I(4711+4740)/I(2854+2868)")
	case 18042:
	    call sprintf (ion_string, max_chars, 
		"Ar V: I(6435+7006)/I(4626)")
	case 19032:
	    call sprintf (ion_string, max_chars, 
		"K IV: I(6102+6796)/I(4511)")
	case 19041:
	    call sprintf (ion_string, max_chars, 
		"K V: I(6223)/I(6349)")
	case 19042:
	    call sprintf (ion_string, max_chars, 
		"K V: I(4123+4163)/I(2515+2495)")
	case 20042:
	    call sprintf (ion_string, max_chars, 
		"Ca V: I(5309+6087)/I(3996)")
	default:
	    call sprintf (ion_string, max_chars, "")
	}

	return (strlen(ion_string))
end


#-------------------------------------------------------------------------------
#  GET_ATOM_STR -	Generate atom/ion string. 

int procedure get_atom_str (atom, ion, atom_str, sz_str)

#  Calling arguments:
int	atom			# I: atomic number
int	ion			# I: ion
char	atom_str[sz_str]	# O: atom/ion string
int	sz_str			# I: size of input string

#  Declarations:
char	at_str[2]		# Atomic symbol
char	ion_str[4]		# Ion symbol
int	strlen()		# length of a string

begin
	call get_atomic_symbol (atom, at_str, 2)

	call get_ion_symbol (ion, ion_str, 4)

	call sprintf (atom_str, sz_str, "%s %s")
	    call pargstr (at_str)
	    call pargstr (ion_str)

	return ( strlen(atom_str) )
end


#-------------------------------------------------------------------------------
#  GET_ATOMIC_SYMBOL -	Generate atomic symbol given an atomic number. 

procedure get_atomic_symbol (atom, symbol, sz_str)

#  Calling arguments:
int	atom		# I: atomic number
char	symbol[sz_str]	# O: atomic symbol
int	sz_str		# I: size of input string

begin
	switch (atom) {
	case HYDROGEN:
	    call strcpy ("H", symbol, sz_str)

	case HELIUM:
	    call strcpy ("He", symbol, sz_str)

	case CARBON:
	    call strcpy ("C", symbol, sz_str)

	case NITROGEN:
	    call strcpy ("N", symbol, sz_str)

	case OXYGEN:
	    call strcpy ("O", symbol, sz_str)

	case NEON:
	    call strcpy ("Ne", symbol, sz_str)

	case SODIUM:
	    call strcpy ("Na", symbol, sz_str)

	case MAGNESIUM:
	    call strcpy ("Mg", symbol, sz_str)

	case ALUMINUM:
	    call strcpy ("Al", symbol, sz_str)

	case SILICON:
	    call strcpy ("Si", symbol, sz_str)

	case SULFUR:
	    call strcpy ("S", symbol, sz_str)

	case CHLORINE:
	    call strcpy ("Cl", symbol, sz_str)

	case ARGON:
	    call strcpy ("Ar", symbol, sz_str)

	case POTASSIUM:
	    call strcpy ("K", symbol, sz_str)

	case CALCIUM:
	    call strcpy ("Ca", symbol, sz_str)

	case IRON:
	    call strcpy ("Fe", symbol, sz_str)

	case NICKEL:
	    call strcpy ("Ni", symbol, sz_str)

	case SELENIUM:
	    call strcpy ("Se", symbol, sz_str)

	case KRYPTON:
	    call strcpy ("Kr", symbol, sz_str)

	default:
	    call error (INVLD_ATOM_ION, "Unrecognized atomic symbol")
	}
end


#-------------------------------------------------------------------------------
#  GET_ION_SYMBOL -	Generate ion symbol given an ion number. 

procedure get_ion_symbol (ion, symbol, sz_str)

#  Calling arguments:
int	ion		# I: atomic number
char	symbol[sz_str]	# O: atomic symbol
int	sz_str		# I: size of input string

begin
	switch (ion+1) {
	case 1:
	    call strcpy ("I", symbol, sz_str)

	case 2:
	    call strcpy ("II", symbol, sz_str)

	case 3:
	    call strcpy ("III", symbol, sz_str)

	case 4:
	    call strcpy ("IV", symbol, sz_str)

	case 5:
	    call strcpy ("V", symbol, sz_str)

	case 6:
	    call strcpy ("VI", symbol, sz_str)

	case 7:
	    call strcpy ("VII", symbol, sz_str)

	case 8:
	    call strcpy ("VIII", symbol, sz_str)

	case 9:
	    call strcpy ("IX", symbol, sz_str)

	case 10:
	    call strcpy ("X", symbol, sz_str)

	default:
	    call error (INVLD_ATOM_ION, "Unrecognized ion symbol")
	}
end


