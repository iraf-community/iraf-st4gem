# RDTABHDB -- Read a boolean keyword from the table header

bool procedure rdtabhdb (tp, keyword)

pointer	tp		# i: table descriptor
char	keyword[ARB]	# i: keyword name
#--
bool	value
pointer	htp

bool	tbhgtb()
pointer	rdtabnew()

begin
	htp = rdtabnew (tp)
	if (htp == NULL) {
	    value = tbhgtb (tp, keyword)
	} else {
	    value = tbhgtb (htp, keyword)
	    call tbtclo (htp)
	}

	return (value)
end

# RDTABHDD -- Read a double keyword from the table header

double procedure rdtabhdd (tp, keyword)

pointer	tp		# i: table descriptor
char	keyword[ARB]	# i: keyword name
#--
double	value
pointer	htp

double	tbhgtd()
pointer	rdtabnew()

begin
	htp = rdtabnew (tp)
	if (htp == NULL) {
	    value = tbhgtd (tp, keyword)
	} else {
	    value = tbhgtd (htp, keyword)
	    call tbtclo (htp)
	}

	return (value)
end

# RDTABHDI -- Read an integer keyword from the table header

int procedure rdtabhdi (tp, keyword)

pointer	tp		# i: table descriptor
char	keyword[ARB]	# i: keyword name
#--
int	value
pointer	htp

int	tbhgti()
pointer	rdtabnew()

begin
	htp = rdtabnew (tp)
	if (htp == NULL) {
	    value = tbhgti (tp, keyword)
	} else {
	    value = tbhgti (htp, keyword)
	    call tbtclo (htp)
	}

	return (value)
end

# RDTABHDR -- Read a real keyword from the table header

real procedure rdtabhdr (tp, keyword)

pointer	tp		# i: table descriptor
char	keyword[ARB]	# i: keyword name
#--
real	value
pointer	htp

real	tbhgtr()
pointer	rdtabnew()

begin
	htp = rdtabnew (tp)
	if (htp == NULL) {
	    value = tbhgtr (tp, keyword)
	} else {
	    value = tbhgtr (htp, keyword)
	    call tbtclo (htp)
	}

	return (value)
end

# RDTABHDT -- Read a text keyword from the table header

procedure rdtabhdt (tp, keyword, value, maxch)

pointer	tp		# i: table descriptor
char	keyword[ARB]	# i: keyword name
char	value[ARB]	# o: keyword value
int	maxch		# i: max value length
#--
pointer	htp

pointer	rdtabnew()

begin
	htp = rdtabnew (tp)
	if (htp == NULL) {
	    call tbhgtt (tp, keyword, value, maxch)
	} else {
	    call tbhgtt (htp, keyword, value, maxch)
	    call tbtclo (htp)
	}
end

# RDTABNEW -- Return a pointer to the primary header if this is an extension

pointer procedure rdtabnew (tp)

pointer	tp		# i: table descriptor
#--
pointer	sp, oldname, newname, htp

bool	streq()
pointer	tbtopn()

begin
	call smark (sp)
	call salloc (oldname, SZ_FNAME, TY_CHAR)
	call salloc (newname, SZ_FNAME, TY_CHAR)

	call tbtnam (tp, Memc[oldname], SZ_FNAME)
	call syntabname (tp, Memc[newname], SZ_FNAME)

	if (streq (Memc[oldname], Memc[newname])) {
	    htp = NULL

	} else {
	    call strcat ("[0]", Memc[newname], SZ_FNAME)
	    htp = tbtopn (Memc[newname], READ_ONLY, NULL)
	}

	call sfree (sp)
	return (htp)
end
