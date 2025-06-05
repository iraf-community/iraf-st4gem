include <syserr.h>
include <iraf77.h>

define	SZ_FNT		512
define  MAX_OPNF	25		# Maximun open templates.
define  CH_DELIM 	20B		# used to flag image section

# UIMOTP -- Is the adaptation of IMTOPEN to the IRAF77 interface. The main
#	    difference lie in the treatment of multiple open. We need to keep
#	    track of those list pointers and the flags to count for a first
#	    time extraction in 'uimxtp' to handle templates with multigroup
#	    files (as in filenam[*])
#
#
# IMTOPEN -- Open an image template.  The filename template package is
# sophisticated enough to do all the necessary filename editing, etc., so all
# we need do is recast the image notation into a FNT edit operation, e.g.,
# `*.imh[*,-*]' becomes `*.hhh%%?\[\*\,-\*]%', with the ? (CH_DELIM, actually
# an unprintable ascii code) being included to make it easy to locate the
# section string in the filenames returned by FNT.  We then open the resultant
# template and perform the inverse mapping upon the filenames returned by FNT.

procedure uimotp (f77tpl, listp, istat)

%	character*(*)	f77tpl
int		istat

char	template[SZ_FNT]		# image template
int	sort, level, ip, ch, ier
pointer	sp, listp, fnt, op
define	output {Memc[op]=$1;op=op+1}
int	fntopnb(), strlen(), errcode()
int	imtbuf[MAX_OPNF]
bool	flags[MAX_OPNF]
int	count
int	dum1[MAX_OPNF]
int	dum2[MAX_OPNF]

common	/tpcom/imtbuf,flags,count,dum1,dum2
#data	count/0/   # Commented out April 26 1993 NZ. Failed to compile
		   # on an IBM/RISC machine (Common BLock thing)

begin
	count = 0    # April 26 '93. This routine gets call once per process.
		     # Hopefully it will not break anything. 
	istat = ER_OK
	call f77upk (f77tpl, template, SZ_FNT)
	call smark (sp)
	call salloc (fnt, strlen(template)*12/10 + SZ_FNT, TY_CHAR)

	# Sorting is disabled as input and output templates, derived from the
	# same database but with string editing used to modify the output list,
	# may be sorted differently as sorting is performed upon the edited
	# output list.

	sort = NO

	op = fnt
	for (ip=1;  template[ip] != EOS;  ip=ip+1) {
	    ch = template[ip]

	    if (ch == '[') {
		if (ip > 1 && template[ip-1] == '!') {
		    # ![ -- Pass a [ to FNT (character class notation).
		    Memc[op-1] = '['

		} else if (ip > 1 && template[ip-1] == '\\') {
		    # \[ -- The [ is part of the filename.  Pass it on as an
		    # escape sequence to get by the FNT.

		    output ('[')

		} else {
		    # [ -- Unescaped [.  This marks the beginning of an image
		    # section sequence.  Output `%%[...]%' and escape all
		    # pattern matching metacharacters until a comma template
		    # delimiter is encountered.  Note that a comma within []
		    # is not a template delimiter.

		    output ('%')
		    output ('%')
		    output (CH_DELIM)

		    level = 0
		    for (;  template[ip] != EOS;  ip=ip+1) {
			ch = template[ip]
			if (ch == ',') {		# ,
			    if (level <= 0)
				break			# exit loop
			    else {
				output ('\\')
				output (ch)
			    }
			} else if (ch == '[') {		# [
			    output ('\\')
			    output (ch)
			    level = level + 1
			} else if (ch == ']') {		# ]
			    output (ch)
			    level = level - 1
			} else if (ch == '*') {		# *
			    output ('\\')
			    output (ch)
			} else				# normal chars
			    output (ch)
		    }
		    output ('%')
		    ip = ip - 1
		}

	    } else if (ch == '@') {
		# List file reference.  Output the CH_DELIM code before the @
		# to prevent further translations on the image section names
		# returned from the list file, e.g., "CH_DELIM // @listfile".

		output (CH_DELIM)
		output ('/')
		output ('/')
		output (ch)

	    } else 
		output (ch)
	}

	Memc[op] = EOS
	iferr (listp = fntopnb (Memc[fnt], sort)) {
		ier = errcode()
		if (ier == SYS_FNTMAGIC)
		   istat = ER_IMTOPN
		else if (ier == SYS_FNTMAXPAT)
		   istat = ER_IMTMAXPAT
		else
		   istat = ER_IMTILLPAT
	}

	count = count + 1
	imtbuf[count] = listp
	flags[count] = true	# not reach uimxtp yet
	call sfree (sp)
	return
end
