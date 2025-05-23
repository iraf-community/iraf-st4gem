# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<lexnum.h>
include	<evvexpr.h>

#--------------------------------------------------------------------6 Feb 98--
.help params Aug97 nebular/lib
.ih
NAME
params -- Parameter/value handling
.ih
DESCRIPTION
The following library handles the setting of arbitrary "parameters" or
variables.

There are two forms of input/output.  The first is the direct setting
of a parameter or retrieval of a parameter.  The second is the use of
strings for setting/inquiring.  Specifically, a string with the format:

.nf
	param value param value ...

or alternatively with the restricted format:

	param=value param=value ...
.fi

can be used to set parameter values, or to retrieve parameter values.
Note that parameter/value pairs are either separated by spaces or by 
equals signs, depending upon how the "restrict_syntax" value is set in 
the call to pr_alloc.  If a string value contains a space, the value 
should be quoted using double quotes.  To a place double quote within 
a value, escape it by preceeding the character with a '\'.  A 
backslash can be included by specifying "\\".

Note that all parameter values are stored as strings.  Conversion
to/from strings are performed for all the type-dependent calls.  The
parameter names are case-insensitive.

.ls OBJECT MEMORY STRUCTURE
The following is the memory structure of the parameter object.  The
primary variabls are listed first.
.ls PNAMES
The string containing the list of parameters.  The position of a
parameter name in the list determines where the value of that
parameter is in the VALUES array.
.le
.ls VALUE
The value of a parameter.
.le
.ls NPARS
The number of parameters defined.
.le

The following variables/macros help manage the memory object.
.ls PN_PTR
Pointer to the PNAMES string.
.le
.ls VA_PTR
Pointer to the V_PTR array.
.le
.ls V_PTR
The array of pointer to the values.
.le
.ls PN_SZ
Current maximum length of the PNAMES string
.le
.ls V_SZ
Current maximum number of values that can be stored.
.le
.ls PN_MAX_SZ
Current maximum length of a parameter name.
.le
.ls GROW
The amount the V_PTR array grows.
.le
.ls CURRENT
Index of the "current" parameter.  This is the index of the parameter
returned by the last pr_nextX call.
.le
.ls RESTRICT_SYNTAX
Is the syntax of the param/value pairs to be restricted?  If so, the 
pairs must be separated with an equals sign.  
.le

.ih OBJECT ROUTINE SUMMARY
.ls Public Routines
The set of routines that make up the public interface of the parameter 
object are as follows:

.nf
     pr_alloc -- Create a parameter object
      pr_free -- Free a parameter object

      pr_sset -- Set parameter values from a parameter string
       pr_get -- Return parameter values in a string
      pr_getp -- Return parameter values as a string pointer

pr_sset[idsb] -- Set a specific parameter
 pr_get[idsb] -- Return value of a specific parameter
       pr_gsp -- Return pointer to value string

       pr_rew -- Rewind to beginning of list
pr_next[idsb] -- Return value of the "next" parameter
    pr_nextsp -- Return pointer to value string of "next" parameter

      pr_copy -- Copy a params object
.fi
.le

.ls Private Routines
.nf
     pr_ssetg -- Convert string parameter value to specified type
     pr_c2bin -- Convert parameter value of specified type given index
     pr_index -- Return parameter index given its name
    pr_initop -- 
      pr_vtoc -- 
    pr_npname -- Get the next parameter name
      pr_test -- Test of parameter object functionality
.fi
.le

.endhelp
#---------------------------------------------------------------------------
# The params object memory structure
define	PN_PTR			Memi[$1+0]
define	PN_SZ			Memi[$1+1]
define	VA_PTR			Memi[$1+2]
define	V_SZ			Memi[$1+3]
define	NPARS			Memi[$1+4]
define	PN_MAX_SZ		Memi[$1+5]
define	CURRENT			Memi[$1+6]
define	RESTRICT_SYNTAX		Memb[$1+7]
define	SZ_PARS			8

define	PNAMES			Memc[PN_PTR($1)]
define	V_PTR			Memi[VA_PTR($1)+$2-1]
define	GROW			10
#---------------------------------------------------------------------------
# End of the params object memory structure.
#---------------------------------------------------------------------------
# Generic memory management
define	Pname			Memc[pname]
define	Sx			Memc[sx]
define	Value			Memc[value]
define	Vin			Memc[vin+$1-1]
define	Vout			Memc[vout+$1-1]
#
# The dreaded gotos
define	def_param_		10
#---------------------------------------------------------------------------
pointer procedure pr_alloc (restrict_syntax)

pointer	o			# The params object.

# Declarations
bool	restrict_syntax		#I: Use '=' to separate parameter/value pairs?

errchk	malloc

begin
	call malloc (o, SZ_PARS, TY_STRUCT)
	
	call malloc (PN_PTR(o), 1, TY_CHAR)
	call strcpy ("", PNAMES(o), 1)
	PN_SZ(o) = 0
	PN_MAX_SZ(o) = 0

	VA_PTR(o) = NULL
	V_SZ(o) = 0

	NPARS(o) = 0
	CURRENT(o) = 0
	RESTRICT_SYNTAX(o) = restrict_syntax

	return (o)
end
#---------------------------------------------------------------------------
# End of pr_alloc
#---------------------------------------------------------------------------
procedure pr_free (o)

pointer	o			# IO: The params object, NULL on return.

# Declarations
int	i			# Generic.

errchk	evvfree, mfree

begin
	# If there is no object, return.
	if (o == NULL)
	    return

	# Free the values
	do i = 1, NPARS(o)
	    call evvfree (V_PTR(o,i))

	# Free the primary arrays.
	call mfree (VA_PTR(o), TY_POINTER)
	call mfree (PN_PTR(o), TY_CHAR)

	# Free the object.
	call mfree (o, TY_STRUCT)
end
#---------------------------------------------------------------------------
# End of pr_free
#---------------------------------------------------------------------------
procedure pr_sset (o, s)

pointer	o			# I:  The params object.
char	s[ARB]			# I:  Parameter string.

# Declarations
bool	add			# 'true' to character to build value.
int	b			# Current character in build value.
int	c			# Current character.
int	i			# Generic.
bool	p_found			# 'true' if a parameter name has been found.
pointer	pname			# Parameter name.
bool	quote			# 'true' if in a quoted string.
int	strlen()		# String length.
pointer	vout			# Build string.

errchk	malloc, mfree, pr_ssets, realloc

begin
	c = 1
	b = 0
	p_found = false
	quote = false
	add = false
	i = strlen (s)
	call malloc (pname, i, TY_CHAR)
	call malloc (vout, i, TY_CHAR)
	
	# Go through the parameter string, find each parameter/value
	# pair, then go and set that parameter to the value.
	while (s[c] != EOS) {
	    switch (s[c]) {
            case '\\':

                # An escape.  Add the next character to whatever is being
                # built.
                c = c + 1
                add = true
                if (s[c] == EOS)
                    break
                
	    case '"':

                # Quoting a value.
                # If a parameter has not been found, there is a serious
                # problem with the input.
                if (!p_found)
                    call error (1, "syntax error: value without parameter")

                # If initial quote, reset beginning of value to next character.
                if (!quote) {
                    b = 0
                    add = false
                    quote = true
                }

                # Else, if ending quote, set the parameter value.
                else
                    goto def_param_
                
	    case '=':

		# If the more restrictive syntax is imposed, a parameter name 
		# has been found.  Otherwise, treat it as any other character.
		if (RESTRICT_SYNTAX(o)) {
                    call strcpy (Vout(1), Pname, b)
		    p_found = true
                    b       = 0
		    add     = false
		}

	    case ' ':

		# The global delimiter if not quoted.
		add = true
		if (!quote) {

		    # If a parameter has already been found, then the end
		    # of the value has been detected.  Set the parameter.
		    if (p_found) {

def_param_
			# Found end of parameter value, set the parameter.
			Vout(b+1) = EOS
			call pr_ssetg (o, Pname, Vout(1), INDEFI)
			b = 0
			p_found = false
			quote = false
			add = false
		    }

		    # Else, check to see if any parameter has been found.
		    # If not skip the space.
		    else if (b == 0) {
			add = false
		    }

		    # Else, the end of the parameter has been found.
		    else {
                        call strcpy (Vout(1), Pname, b)
                        b = 0
                        p_found = true
                        add = false
		    }
		}

            default:

                # Add the character to the build string.
                add = true
	    }

            # Add the character to the build string.
            if (add) {
                b = b + 1
                Vout(b) = s[c]
            }
            
	    # Increment the character pointer.
	    c = c + 1
	}

	# Make sure the final parameter is set.
	if (p_found) {
            Vout(b+1) = EOS
	    call pr_ssetg (o, Pname, Vout(1), INDEFI)
        }

	# That's all folks.
	call mfree (pname, TY_CHAR)
	call mfree (vout, TY_CHAR)
end
#---------------------------------------------------------------------------
# End of pr_sset
#---------------------------------------------------------------------------
pointer procedure pr_getp (o)

pointer	o			# I:  Params object

# Declarations
bool	first			# First parameter.
int	i, j			# Generic.
pointer	pname			# Parameter name.
bool	pr_nextsp()		# Get next parameter as a string pointer.
bool	quote			# True if the value should be quoted.
pointer	s			# String to return.
pointer	sb			# String buffer.
pointer	sb_open()		# Create a string buffer.
pointer	sb_string()		# Return the string buffer.
int	stridx()		# Find character in string.
int	strlen()		# Length of string.
pointer	value			# Parameter value.
pointer	vin			# Pointer to stored value.
pointer	vout			# Version of value with escapes.

errchk	malloc, mfree, sb_cat, sb_close, sb_open, sb_string

begin
	# Create the output string buffer.
	sb = sb_open()
	
	# Loop through all the parameters, building up the string.
	first = true
	call pr_rew (o)
	pname = NULL
	value = NULL
	while (pr_nextsp (o, pname, value)) {

	    # Append a space if there is more than 1 parameter.
	    if (!first)
		call sb_cat (sb, " ")
	    
	    # See if there are any spaces in the value.  If so, quote
	    # the whole value.
	    quote = (stridx (" ", Value) > 0 ||
		     strlen (Value) <= 0)

	    # Build the string.  Put the parameter name in the string.
	    call sb_cat (sb, Pname)

	    # For the restricted syntax, insert an equals sign between 
	    # parameter/value pairs.
	    if (RESTRICT_SYNTAX(o))
	    	call sb_cat (sb, "=")
	    else
	    	call sb_cat (sb, " ")

	    # Place the value in the string.
	    if (quote)
		call sb_cat (sb, "\"")

            # If there are no quotes in the value, just concatentate it.
            # If there are, then escape the quotes.
            if (stridx ("\"", Value) > 0) {
                call malloc (vout, 2 * strlen (Value), TY_CHAR)
                vin = value
                i = 1
                j = 1
                while (Vin(i) != EOS) {
                    if (Vin(i) == '"') {
                        Vout(j) = '\\'
                        j = j + 1
                    }
                    Vout(j) = Vin(i)
                    j = j + 1
                    i = i + 1
                }
                Vout(j) = EOS
                call sb_cat(sb, Vout(1))
                call mfree (vout, TY_CHAR)
            } else
                call sb_cat (sb, Value)

            # If quoting, close quote.
	    if (quote)
		call sb_cat (sb, "\"")

	    first = false
	}

	# That's all folks.
	s = sb_string (sb)
	call sb_close (sb)
	call mfree (pname, TY_CHAR)
	call mfree (value, TY_CHAR)
	return (s)
end
#---------------------------------------------------------------------------
# End of pr_getp
#---------------------------------------------------------------------------
procedure pr_get (o, out, max)

pointer	o			# I:  Params object.
char	out[ARB]		# O:  String containing parameters
int	max			# O:  Maximum length of string.

# Declarations
pointer	pr_getp()		# Get parameter values.
pointer	value			# Value string.

errchk	mfree, pr_getp

begin
	value = pr_getp (o)
	call strcpy (Value, out, max)
	call mfree (value, TY_CHAR)
end
#---------------------------------------------------------------------------
# End of pr_get
#---------------------------------------------------------------------------
procedure pr_ssets (o, pname, value)

pointer	o			# I:  Params object.
char	pname[ARB]		# I:  Parameter name.
char	value[ARB]		# I:  Parameter value.

begin
	call pr_ssetg (o, pname, value, TY_CHAR)
end
#---------------------------------------------------------------------------
# End of pr_ssets
#---------------------------------------------------------------------------
procedure pr_sseti (o, name, value)

pointer	o			# I:  Params object.
char	name[ARB]		# I:  Parameter name.
int 	value			# I:  Parameter value.

# Declarations
int	index			# parameter index.
int	pr_index()		# Get parameter index.

begin
	index = pr_index (o, name, true)
	call pr_initop (V_PTR(o,index), TY_INT)
	O_VALI(V_PTR(o,index)) = value
end
#---------------------------------------------------------------------------
# End of pr_sseti
#---------------------------------------------------------------------------
procedure pr_ssetb (o, name, value)

pointer	o			# I:  Params object.
char	name[ARB]		# I:  Parameter name.
bool 	value			# I:  Parameter value.

# Declarations
int	btoi()			# Convert boolean to integer.
int	index			# parameter index.
int	pr_index()		# Get parameter index.

begin
	index = pr_index (o, name, true)
	call pr_initop (V_PTR(o,index), TY_BOOL)
	O_VALI(V_PTR(o,index)) = btoi (value)
end
#---------------------------------------------------------------------------
# End of pr_ssetb
#---------------------------------------------------------------------------
procedure pr_ssetd (o, name, value)

pointer	o			# I:  Params object.
char	name[ARB]		# I:  Parameter name.
double	value			# I:  Parameter value.

# Declarations
int	index			# parameter index.
int	pr_index()		# Get parameter index.

begin
	index = pr_index (o, name, true)
	call pr_initop (V_PTR(o,index), TY_DOUBLE)
	O_VALD(V_PTR(o,index)) = value
end
#---------------------------------------------------------------------------
# End of pr_ssetd
#---------------------------------------------------------------------------
pointer procedure pr_gsp (o, name)

pointer	o			# I:  Params object.
char	name[ARB]		# I:  Parameter to get value of.

pointer	value			# O:  Value of parameter as string pointer.

# Declarations
int	l			# Length of string.
int	p			# Parameter index.
int	pr_index()		# Get parameter index.
pointer	pr_vtoc()		# Get value as a string.
pointer	sp			# Stack pointer.
int	strlen()		# Length of string.
pointer	sx			# Generic.

errchk	malloc
errchk	salloc, sfree, smark

begin
	call smark (sp)

	l = strlen (name)
	call salloc (sx, SZ_LINE+l, TY_CHAR)
	
	# Find the parameter.
	p = pr_index (o, name, false)
	if (IS_INDEFI(p)) {
	    call sprintf (Sx, SZ_LINE+strlen (name), "parameter %s not defined")
	    call pargs (name)
	    call error (1, Sx)
	}

	# Return the value.
	value = pr_vtoc (V_PTR(o,p))
	
	# That's all folks.
	call sfree (sp)
	return (value)
end
#---------------------------------------------------------------------------
# End of pr_gsp
#---------------------------------------------------------------------------
procedure pr_gets (o, name, v, max)

pointer	o			# I:  Params object.
char	name[ARB]		# I:  Parameter name.
char	v[max]			# O:  Value of parameter.
int	max			# I:  Maximum length of output string.

# Declarations
pointer	pr_gsp()		# Get pointer to parameter value.
pointer	value			# Pointer to parameter value.

errchk	mfree, pr_gsp

begin
	value = pr_gsp (o, name)
	call strcpy (Value, v, max)
	call mfree (value, TY_CHAR)
end
#---------------------------------------------------------------------------
# End of pr_gets
#---------------------------------------------------------------------------
int procedure pr_geti (o, name)

pointer	o			# I:  Params object.
char	name[ARB]		# I:  Parameter name.

int	v			# O:  Value.

# Declarations
int	ctoi()			# String to integer.
int	i			# Generic.
int	index			# Index of parameter.
pointer	op			# Operand value.
int	pr_index()		# Get parameter index.
char	sx[SZ_LINE]		# Generic.

begin
	index = pr_index (o, name, false)
	if (IS_INDEFI (index)) {
	    call sprintf (sx, SZ_LINE, "parameter '%s' not found")
	    call pargstr (name)
	    call error (1, sx)
	}
	op = V_PTR(o,index)
	switch (O_TYPE(op)) {
	case TY_CHAR:
	    i = 1
	    i = ctoi (O_VALC(op), i, v)

	case TY_SHORT:
	    v = O_VALS(op)

	case TY_INT, TY_BOOL:
	    v = O_VALI(op)

	case TY_LONG:
	    v = O_VALL(op)

	case TY_REAL:
	    v = O_VALR(op)

	case TY_DOUBLE:
	    v = O_VALD(op)

	case TY_POINTER:
	    v = O_VALP(op)

	default:
	    call error (1, "pr_geti: unkown type")
	}	
	    
	return (v)
end
#---------------------------------------------------------------------------
# End of pr_geti
#---------------------------------------------------------------------------
bool procedure pr_getb (o, name)

pointer	o			# I:  Params object.
char	name[ARB]		# I:  Parameter name.

bool	v			# O:  Value.

# Declarations
int	ctoi()			# String to integer
int	i, j			# Generic.
int	index			# Index of parameter.
bool	itob()			# Convert integer to boolean.
pointer	op			# Operand value.
int	pr_index()		# Get parameter index.
char	sx[SZ_LINE]		# Generic.

begin
	index = pr_index (o, name, false)
	if (IS_INDEFI (index)) {
	    call sprintf (sx, SZ_LINE, "parameter '%s' not found")
	    call pargstr (name)
	    call error (1, sx)
	}
	op = V_PTR(o,index)
	switch (O_TYPE(op)) {
	case TY_CHAR:
	    i = 1
	    i = ctoi (O_VALC(op), i, j)
	    v = itob (j)

	case TY_SHORT:
	    v = O_VALS(op) > 1

	case TY_INT, TY_BOOL:
	    v = itob (O_VALI(op))

	case TY_LONG:
	    v = O_VALL(op) > 1

	case TY_REAL:
	    v = O_VALR(op) > 1

	case TY_DOUBLE:
	    v = O_VALD(op) > 1

	case TY_POINTER:
	    v = O_VALP(op) > 1

	default:
	    call error (1, "pr_getd: unkown type")
	}	
	    
	return (v)
end
#---------------------------------------------------------------------------
# End of pr_getb
#---------------------------------------------------------------------------
double procedure pr_getd (o, name)

pointer	o			# I:  Params object.
char	name[ARB]		# I:  Parameter name.

double	v			# O:  Value.

# Declarations
int	ctod()			# String to double.
int	i			# Generic.
int	index			# Index of parameter.
pointer	op			# Operand value.
int	pr_index()		# Get parameter index.
char	sx[SZ_LINE]		# Generic.

begin
	index = pr_index (o, name, false)
	if (IS_INDEFI (index)) {
	    call sprintf (sx, SZ_LINE, "parameter '%s' not found")
	    call pargstr (name)
	    call error (1, sx)
	}
	op = V_PTR(o,index)
	switch (O_TYPE(op)) {
	case TY_CHAR:
	    i = 1
	    i = ctod (O_VALC(op), i, v)

	case TY_SHORT:
	    v = O_VALS(op)

	case TY_INT, TY_BOOL:
	    v = O_VALI(op)

	case TY_LONG:
	    v = O_VALL(op)

	case TY_REAL:
	    v = O_VALR(op)

	case TY_DOUBLE:
	    v = O_VALD(op)

	case TY_POINTER:
	    v = O_VALP(op)

	default:
	    call error (1, "pr_getd: unkown type")
	}	
	    
	return (v)
end
#---------------------------------------------------------------------------
# End of pr_getd
#---------------------------------------------------------------------------
bool procedure pr_nextsp (o, pname, value)

pointer	o			# I:  Params object.
pointer	pname			# IO: Pointer to string containing parameter.
pointer	value			# IO: Pointer to string containing value.

# Declarations.
pointer	pr_gsp()		# Get value as string pointer.
bool	pr_npname()		# Get next parameter name.

begin
	if (!pr_npname (o, pname))
	    return (false)
	
	if (value != NULL)
	    call mfree (value, TY_CHAR)
	value = pr_gsp (o, Pname)

	return (true)
end
#---------------------------------------------------------------------------
# End of pr_nextsp
#---------------------------------------------------------------------------
bool procedure pr_nexts (o, p, pmax, v, vmax)

pointer	o			# I:  Params object.
char	p[pmax]			# O:  Parameter name.
int	pmax			# I:  Maximum length of pname.
char	v[vmax]			# O:  Parameter value.
int	vmax			# I:  Maximum length of value.

# Declarations.
pointer	pname			# Parameter name.
bool	pr_nextsp()		# Get next parameter as string.
pointer	value			# Parameter value as string.

errchk	mfree, pr_nextsp

begin
	# Get next paramter as string.  If no more parameters, return.
	pname = NULL
	value = NULL
	if (!pr_nextsp (o, pname, value))
	    return (false)

	call strcpy (Pname, p, pmax)
	call strcpy (Value, v, vmax)
	
	call mfree (value, TY_CHAR)
	call mfree (pname, TY_CHAR)
	return (true)
end
#---------------------------------------------------------------------------
# End of pr_nexts
#---------------------------------------------------------------------------
bool procedure pr_nexti (o, p, pmax, v, valid)

pointer	o			# I:  Params object.
char	p[pmax]			# O:  Parameter name.
int	pmax			# I:  Maximum length of name.
int	v			# O:  Parameter value.
bool	valid			# O:  'true' if a valid number.

# Declarations.
pointer	pname			# Name of next parameter.
bool	pr_npname()		# Get next parameter name.
int	pr_geti()		# Get integer-valued parameter.

begin
	# Get next paramter name.
	pname = NULL
	if (!pr_npname (o, pname))
	    return (false)

	# Copy name.
	call strcpy (Pname, p, pmax)
	
	# Get the value.
	v = pr_geti (o, Pname)
	valid = true

	# That's all folks.
	call mfree (pname, TY_CHAR)
	return (true)
end
#---------------------------------------------------------------------------
# End of pr_nexti
#---------------------------------------------------------------------------
bool procedure pr_nextb (o, p, pmax, v, valid)

pointer	o			# I:  Params object.
char 	p[pmax]			# O:  Parameter name.
int	pmax			# I:  Maximum length of name.
bool	v			# O:  Parameter value.
bool	valid			# O:  'true' if a valid number.

# Declarations.
pointer	pname			# Name of next parameter.
bool	pr_npname()		# Get next parameter name.
bool	pr_getb()		# Get boolean-valued parameter.

begin
	# Get next paramter name.
	pname = NULL
	if (!pr_npname (o, pname))
	    return (false)

	# Copy name.
	call strcpy (Pname, p, pmax)
	
	# Get the value.
	v = pr_getb (o, Pname)
	valid = true

	# That's all folks.
	call mfree (pname, TY_CHAR)
	return (true)
end
#---------------------------------------------------------------------------
# End of pr_nextb
#---------------------------------------------------------------------------
bool procedure pr_nextd (o, p, pmax, v, valid)

pointer	o			# I:  Params object.
char	p[pmax]			# O:  Parameter name.
int	pmax			# I:  Maximum length of name.
double 	v			# O:  Parameter value.
bool	valid			# O:  'true' if a valid number is found.

# Declarations.
pointer	pname			# Name of next parameter.
bool	pr_npname()		# Get next parameter name.
double 	pr_getd()		# Get double-valued parameter.

begin
	# Get next paramter name.
	pname = NULL
	if (!pr_npname (o, pname))
	    return (false)

	# Copy name.
	call strcpy (Pname, p, pmax)
	
	# Get the value.
	v = pr_getd (o, Pname)
	valid = true

	# That's all folks.
	call mfree (pname, TY_CHAR)
	return (true)
end
#---------------------------------------------------------------------------
# End of pr_nextd
#---------------------------------------------------------------------------
procedure pr_rew (o)

pointer	o			# I:  Params object.

begin
	CURRENT(o) = 0
end
#---------------------------------------------------------------------------
# End of pr_rew
#---------------------------------------------------------------------------
pointer procedure pr_copy (o)

pointer	o			# I:  Params object to copy.

pointer	n			# O:  New object.

# Declarations.
int	i			# Generic.
pointer	pr_alloc()		# Create a params object.

errchk	calloc, malloc

begin
	n = pr_alloc (RESTRICT_SYNTAX(o))

	PN_SZ(n) = PN_SZ(o)
	V_SZ(n) = V_SZ(o)
	NPARS(n) = NPARS(o)
	PN_MAX_SZ(n) = PN_MAX_SZ(o)
	CURRENT(n) = CURRENT(o)

	call realloc (PN_PTR(n), PN_SZ(o), TY_CHAR)
	call strcpy (PNAMES(o), PNAMES(n), PN_SZ(o))

	call realloc (VA_PTR(n), V_SZ(n), TY_POINTER)
	do i = 1, NPARS(n) {
	    call calloc (V_PTR(n,i), LEN_OPERAND, TY_STRUCT)
	    call amovi (Memi[V_PTR(o,i)], Memi[V_PTR(n,i)], LEN_OPERAND)
	    if (O_TYPE(V_PTR(n,i)) == TY_CHAR) {
		call xvv_initop (V_PTR(n,i), O_LEN(V_PTR(o,i)), TY_CHAR)
		call strcpy (O_VALC(V_PTR(o,i)), O_VALC(V_PTR(n,i)),
			     O_LEN(V_PTR(o,i)))
	    }
	}

	return (n)
end
#---------------------------------------------------------------------------
# End of pr_copy
#---------------------------------------------------------------------------
procedure pr_ssetg (o, pname, v, type)

pointer o			# I:  The params object.
char	pname[ARB]		# I:  Parameter to set.
char	v[ARB]			# I:  Value to set.
int	type			# I:  Type to force value to.

# Declarations
int	index			# Index of the parameter.
int	pr_index()		# Get index of the parameter.

begin
	# Get index of parameter.
	index = pr_index (o, pname, true)

	# Place the value.
	call pr_c2bin (v, type, V_PTR(o,index))
	
	# That's all folks.
end
#---------------------------------------------------------------------------
# End of pr_ssetg
#---------------------------------------------------------------------------
procedure pr_c2bin (v, type, op)

char	v[ARB]			# I:  Value as string.
int	type			# I:  Type to convert to.
pointer	op			# IO: Operand value.

# Declarations
int	btoi()			# Convert boolean to integer.
int	ctod(), ctoi(), ctol(), ctor()
int	i, ix			# Generic.
int	l			# Length of value.
int	lexnum()		# Determine whether string is a number.
int	ntype			# Type of output operand.
bool	streq()			# Strings equal.
int	strlen()		# Get string length.

begin
	# If the type is unspecified, determine it.
	l = strlen (v)
	if (IS_INDEFI(type)) {
	    i = 1
	    ntype = lexnum (v, i, ix)
	    if ((ntype != LEX_DECIMAL && ntype != LEX_REAL) ||
		ix != l) {
		ntype = TY_CHAR
		if (streq (v, "yes") || streq (v, "no"))
		    ntype = TY_BOOL
	    } else
		switch (ntype) {
		case LEX_DECIMAL:
		    ntype = TY_INT
		    
		case LEX_REAL:
		    ntype = TY_DOUBLE
		}	
	} else
	    ntype = type

	# Convert the value to the specified type.
	if (ntype == TY_CHAR) {
	    call xvv_initop (op, l, TY_CHAR)
	} else
	    call xvv_initop (op, 0, ntype)
	O_FLAGS(op) = O_FREEVAL+O_FREEOP

	switch (ntype) {
	case TY_BOOL:
	    O_VALI(op) = btoi (streq (v, "yes"))

	case TY_CHAR:
	    call strcpy (v, O_VALC(op), l)

	case TY_SHORT:
	    i = 1
	    i = ctoi (v, i, ix)
	    O_VALS(op) = ix

	case TY_INT:
	    i = 1
	    i = ctoi (v, i, O_VALI(op))

	case TY_LONG:
	    i = 1
	    i = ctol (v, i, O_VALL(op))

	case TY_REAL:
	    i = 1
	    i = ctor (v, i, O_VALR(op))

	case TY_DOUBLE:
	    i = 1
	    i = ctod (v, i, O_VALD(op))

	case TY_POINTER:
	    i = 1
	    i = ctoi (v, i, O_VALP(op))

	default:
	    call error (1, "pr_c2bin: type not supported")
	}
end
#---------------------------------------------------------------------------
# End of pr_c2bin
#---------------------------------------------------------------------------
int procedure pr_index (o, name, allocate)

pointer	o			# I:  The params object.
char	name[ARB]		# I:  Parameter name to get index of.
bool	allocate		# I:  TRUE to allocate a new parameter.

int	index			# O:  Index, INDEFI if no parameter.

# Declarations
int	n_size			# Size of the parameter name.
pointer	pname			# Full name of the parameter.
pointer	sp			# Stack pointer.
int	strcic()		# Get index.
int	strlen()		# Length of string.

errchk	salloc, sfree, smark

begin
	n_size = max (PN_MAX_SZ(o), strlen(name))

	call smark (sp)
	call salloc (pname, n_size, TY_CHAR)
	
	# Check parameter list to see if the parameter already exists.
	index = strcic (name, Pname, n_size, PNAMES(o))
	
	# If the parameter doesn't exist, add it to the list.
	if (index <= 0 && allocate) {
	    NPARS(o) = NPARS(o) + 1

	    # Add the parameter name.
	    if (n_size + strlen (PNAMES(o)) + 1 > PN_SZ(o)) {
		PN_SZ(o) = PN_SZ(o) + SZ_LINE + n_size + 1
		call realloc (PN_PTR(o), PN_SZ(o), TY_CHAR)
	    }
	    call strcat (",", PNAMES(o), PN_SZ(o))
	    call strcat (name, PNAMES(o), PN_SZ(o))
	    PN_MAX_SZ(o) = n_size

	    # Add a spot to the values array if necessary.
	    if (NPARS(o) > V_SZ(o)) {
		V_SZ(o) = V_SZ(o) + GROW
		call realloc (VA_PTR(o), V_SZ(o), TY_POINTER)
	    }

	    index = NPARS(o)
	    call calloc (V_PTR(o,index), LEN_OPERAND, TY_STRUCT)
	    O_TYPE(V_PTR(o,index)) = 0
	}

	# That's all folks.
	call sfree (sp)
	return (index)
end
#---------------------------------------------------------------------------
# End of pr_index
#---------------------------------------------------------------------------
procedure pr_initop (op, type)
pointer	op			# I:  Operand to initialize.
int	type			# I:  Type to initialize.

begin
	call xvv_initop (op, 0, type)
	O_FLAGS(op) = O_FREEVAL + O_FREEOP
end
#---------------------------------------------------------------------------
# End of pr_initop
#---------------------------------------------------------------------------
pointer procedure pr_vtoc (op)

pointer	op			# I:  Operand to convert to string.

pointer	value			# O:  String version of value.

# Declarations
int	i			# generic.
int	strlen()		# Get string length.

begin
	if (O_TYPE(op) == TY_CHAR) {
	    i = strlen (O_VALC(op))
	    call malloc (value, i, TY_CHAR)
	    call strcpy (O_VALC(op), Value, i)
	} else {
	    call malloc (value, SZ_LINE, TY_CHAR)
	    switch (O_TYPE(op)) {
	    case TY_BOOL:
		call sprintf (Value, SZ_LINE, "%b")
		call pargi (O_VALI(op))

	    case TY_SHORT:
		call sprintf (Value, SZ_LINE, "%d")
		call pargs (O_VALS(op))

	    case TY_INT, TY_POINTER:
		call sprintf (Value, SZ_LINE, "%d")
		call pargi (O_VALI(op))

	    case TY_LONG:
		call sprintf (Value, SZ_LINE, "%d")
		call pargl (O_VALL(op))

	    case TY_REAL:
		call sprintf (Value, SZ_LINE, "%g")
		call pargr (O_VALR(op))

	    case TY_DOUBLE:
		call sprintf (Value, SZ_LINE, "%g")
		call pargd (O_VALD(op))
	    }	
	}

	# That's all folks,
	return (value)
end
#---------------------------------------------------------------------------
# End of pr_vtoc
#---------------------------------------------------------------------------
bool procedure pr_npname(o, pname)

pointer	o			# I:  The params object.
pointer	pname			# IO: The next parameter name.

# Declarations.
int	i			# Generic.
int	word_find()		# Find specific word in list.

errchk	realloc

begin
	# If the last parameter has been retrieved, return false.
	if (CURRENT(o) >= NPARS(o)) {
	    CURRENT(o) = 0
	    return (false)
	}

	# Return the next parameter/value pair
	CURRENT(o) = CURRENT(o) + 1

	call realloc (pname, PN_MAX_SZ(o), TY_CHAR)
	i = word_find (CURRENT(o), PNAMES(o), Pname, PN_MAX_SZ(o))

	return (true)
end
#---------------------------------------------------------------------------
# End of pr_npname
#---------------------------------------------------------------------------
procedure pr_test()

bool	bx, clgetb(), pr_getb(), pr_nexti(), streq()
int	i, strlen(), ix, clgeti(), pr_geti()
pointer	p, pr_alloc(), s, pr_getp(), pn, pr_copy()
char	sx[SZ_COMMAND], sy[SZ_COMMAND]
double	dx, clgetd(), pr_getd()

begin
	call printf ("pr_test: Allocate/deallocate test.\n")
	call printf ("\tAllocating params object.\n")
	p = pr_alloc (false)
	call printf ("\tDeallocating params object.\n")
	call pr_free (p)
	call printf ("\tRe-deallocating params object.\n")
	call pr_free (p)
	call printf ("\tTest done.\n")

	call printf ("pr_test: String operations.\n")
	p = pr_alloc (false)
	call clgstr ("p1", sx, SZ_COMMAND)
	call printf ("\tSetting params from string.\n")
	call pr_sset (p, sx)
	call printf ("\tExamining structure.\n")
	call printf ("\t\tPN_SZ = %d V_SZ = %d NPARS = %d PN_MAX_SZ = %d\n")
	call pargi (PN_SZ(p))
	call pargi (V_SZ(p))
	call pargi (NPARS(p))
	call pargi (PN_MAX_SZ(p))
	call printf ("\t\tCURRENT = %d RESTRICT_SYNTAX = %b\n")
	call pargi (CURRENT(p))
	call pargb (RESTRICT_SYNTAX(p))
	call printf ("\t\tPNAMES = '%s'\n")
	call pargstr (PNAMES(p))
	do i = 1, NPARS(p) {
	    call printf ("\t\tV[%d]: O_TYPE = %d O_LEN = %d O_VAL =")
	    call pargi (i)
	    call pargi (O_TYPE(V_PTR(p,i)))
	    call pargi (O_LEN(V_PTR(p,i)))
	    switch (O_TYPE(V_PTR(p,i))) {
	    case TY_BOOL:
		call printf ("%b\n")
		call pargi (O_VALI(V_PTR(p,i)))
	    case TY_CHAR:
		call printf ("%s\n")
		call pargstr (O_VALC(V_PTR(p,i)))
	    case TY_SHORT:
		call printf ("%d\n")
		call pargs (O_VALS(V_PTR(p,i)))
	    case TY_INT:
		call printf ("%d\n")
		call pargi (O_VALI(V_PTR(p,i)))
	    case TY_LONG:
		call printf ("%d\n")
		call pargl (O_VALL(V_PTR(p,i)))
	    case TY_REAL:
		call printf ("%g\n")
		call pargr (O_VALR(V_PTR(p,i)))
	    case TY_DOUBLE:
		call printf ("%g\n")
		call pargd (O_VALD(V_PTR(p,i)))
	    case TY_POINTER:
		call printf ("%d\n")
		call pargi (O_VALP(V_PTR(p,i)))
	    }	
	}
	s = pr_getp (p)
	call printf ("\tString = '%s'\n")
	call pargstr (Memc[s])
	i = strlen (Memc[s])
	Memc[s] = EOS
	call pr_get (p, Memc[s], i)
	call printf ("\tString2 = '%s'\n")
	call pargstr (Memc[s])
	call mfree (s, TY_CHAR)

	call clgstr ("sname", sx, SZ_COMMAND)
	call clgstr ("sval", sy, SZ_COMMAND)
	call pr_ssets (p, sx, sy)
	call pr_get (p, sy, SZ_COMMAND)
	call printf ("\n\tAdding a string parameter: '%s'\n")
	call pargstr (sy)
	call pr_gets (p, sx, sy, SZ_COMMAND)
	call printf ("\tpr_gets = '%s'\n")
	call pargstr (sy)
	
	call clgstr ("bname", sx, SZ_COMMAND)
	bx = clgetb ("bval")
	call pr_ssetb (p, sx, bx)
	call pr_get (p, sy, SZ_COMMAND)
	call printf ("\tAdding a bool parameter: '%s'\n")
	call pargstr (sy)
	bx = pr_getb (p, sx)
	call printf ("\tpr_getb = %b\n")
	call pargb (bx)
	
	call clgstr ("iname", sx, SZ_COMMAND)
	ix = clgeti ("ival")
	call pr_sseti (p, sx, ix)
	call pr_get (p, sy, SZ_COMMAND)
	call printf ("\tAdding a int parameter: '%s'\n")
	call pargstr (sy)
	ix = pr_geti (p, sx)
	call printf ("\tpr_geti = %d\n")
	call pargi (ix)
	
	call clgstr ("dname", sx, SZ_COMMAND)
	dx = clgetd ("dval")
	call pr_ssetd (p, sx, dx)
	call pr_get (p, sy, SZ_COMMAND)
	call printf ("\tAdding a double parameter: '%s'\n")
	call pargstr (sy)
	dx = pr_getd (p, sx)
	call printf ("\tpr_getd = %g\n")
	call pargd (dx)

	call printf ("\n\tTesting pr_nexti:\n")
	call pr_rew (p)
	while (pr_nexti (p, sx, SZ_COMMAND, ix, bx)) {
	    call printf ("\t\tpname ='%s' ix = %d\n")
	    call pargstr (sx)
	    call pargi (ix)
	}

	call printf ("\n\tTesting pr_copy:\n")
	call pr_get (p, sx, SZ_COMMAND)
	call printf ("\t\told = '%s'\n")
	call pargstr (sx)
	pn = pr_copy (p)
	call pr_get (p, sx, SZ_COMMAND)
	call pr_get (pn, sy, SZ_COMMAND)
	call printf ("\t\told = '%s'\n")
	call pargstr (sx)
	call printf ("\t\tnew = '%s'\n")
	call pargstr (sy)
	if (streq (sx, sy))
	    call printf ("\t\tNew copy string is equal.\n")
	else {
	    call printf ("\t\tNew copy string is different.\n")
	}
	call printf ("\n\t\tDeleting original object.\n")
	call pr_free (p)
	call pr_get (pn, sy, SZ_COMMAND)
	call printf ("\t\told = '%s'\n")
	call pargstr (sx)
	call printf ("\t\tnew = '%s'\n")
	call pargstr (sy)
	if (streq (sx, sy))
	    call printf ("\t\tNew copy string is equal.\n")
	else {
	    call printf ("\t\tNew copy string is different.\n")
	}
	p = pn
	
	call pr_free (p)
end
#---------------------------------------------------------------------------
# End of pr_test
#---------------------------------------------------------------------------
