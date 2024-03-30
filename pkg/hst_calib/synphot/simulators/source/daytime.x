# DAYTIME -- Convert time string to julian date

define	SZ_FIELD	10

#* HISTORY *
#* B.Simon	14-Apr-95	original

double procedure daytime (time)

char	time[ARB]	# i: time string
#--
int	ic, year, month, day, hour, minute, sec, mid
pointer	sp, field, saved


string	mlist    "|january|february|march|april|may|june|july|august|\
september|october|november|december|"
string	meridian "|am|pm|"
string	badtime  "Cannot parse time string"

double	jdate()
int	strdic()

begin
	# Allocate temporary array

	call smark (sp)
	call salloc (field, SZ_FIELD, TY_CHAR)
	call salloc (saved, SZ_FIELD, TY_CHAR)

	ic = 1
	Memc[saved] = EOS

	# Month field is string or number

	call timefield (time, ic, month, Memc[field], SZ_FIELD)
	if (IS_INDEFI (month)) {
	    call strlwr (Memc[field])
	    month = strdic (Memc[field], Memc[field], SZ_FIELD, mlist)
	}

	if (month < 1 || month > 12)
	    call printerr_str (badtime, time)
	    
	# Day field

	call timefield (time, ic, day, Memc[field], 2)
	if (day < 1 || day > 31)
	    call printerr_str (badtime, time)

	# Add century if year is less than a hundred

	call timefield (time, ic, year, Memc[field], 4)
	if (year < 50) {
	    year = year + 2000
	} else if (year < 100) {
	    year = year + 1900
	} else if (year < 1901 || year > 2099) {
	    call printerr_str (badtime, time)
	}

	# Hour, minute, and seconds fields

	call timefield (time, ic, hour, Memc[field], 2)
	if (hour < 0 || hour > 23)
	    call printerr_str (badtime, time)

	# Minute and Second fields are optional. If a string is read, 
	# save the string as the meridian

	call timefield (time, ic, minute, Memc[field], 2)
	if (IS_INDEFI(minute)) {
	    minute = 0
	    call strcpy (Memc[field], Memc[saved], SZ_FIELD)
	}

	if (minute < 0 || minute > 59)
	    call printerr_str (badtime, time)

	call timefield (time, ic, sec, Memc[field], 2)
	if (IS_INDEFI(sec) && Memc[saved] == EOS) {
	    sec = 0
	    call strcpy (Memc[field], Memc[saved], SZ_FIELD)
	}

	if (sec < 0 || sec > 59)
	    call printerr_str (badtime, time)

	if (Memc[saved] == EOS) {
	    call timefield (time, ic, mid, Memc[field], 2)
	} else {
	    call strcpy (Memc[saved], Memc[field], SZ_FIELD)
	}
	    
	call strlwr (Memc[field])
	mid = strdic (Memc[field], Memc[field], SZ_FIELD, meridian)

	# Convert the hours field according to the meridian

	switch (mid) {
	case 0:
	    if (Memc[field] != EOS)
		call printerr_str (badtime, time)
	case 1:
	    if (hour == 12)
		hour = hour - 12
	case 2:
	    if (hour != 12)
		hour = hour + 12
	}

	call sfree (sp)
	return (jdate (year, month, day, hour, minute, sec))

end
