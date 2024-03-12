#######################################################################
#
#		Procedure Julian(year, day, hr, min,jd)
#				int year, day, hr, min
#
#		Returns double jd
#
#		Takes a year and day and converts them to Julian dates
#			Works on days after 1900 and maybe not always after 2000
#		
#		Notes
#			Conversion Algorithmn from Van Flandern & Pulkkinen
#			
#		History
#			June 29 1994	Grimes created
#
#
#
#
#
#
#######################################################################
procedure julian(year,day,hr,min,jd)
int year,day,hr,min
double jd


int 	month, date
int 	leapyear



begin


	if ( mod (year, 4) == 0 ) {
		leapyear = YES
	} else {
		leapyear = NO
	}

	

	if (day <= 31) {
		month = 1 
	} else if ( day <= 59 && leapyear == NO || day<=60 && leapyear == YES) {
		month = 2
		day = day - 31
	} else if ( day <= 90 && leapyear == NO ) {
		month = 3
		day = day - 59
	} else if ( day <= 91 && leapyear == YES ) {
		month = 3
		day = day - 60
	} else if ( day <= 120 && leapyear == NO ) {
                month = 4
                day = day - 90
        } else if ( day <= 121 && leapyear == YES ) {
                month = 4
                day = day - 91
        } else if ( day <= 151 && leapyear == NO ) {
                month = 5
                day = day - 120
        } else if ( day <= 152 && leapyear == YES ) {
                month = 5
                day = day - 121
        } else if ( day <= 181 && leapyear == NO ) {
                month = 6
                day = day - 151
        } else if ( day <= 182 && leapyear == YES ) {
                month = 6
                day = day - 152
        } else if ( day <= 212 && leapyear == NO ) {
                month = 7
                day = day - 181
        } else if ( day <= 213 && leapyear == YES ) {
                month = 7
                day = day - 182
	} else if ( day <= 243 && leapyear == NO ) {
                month = 8
                day = day - 212
        } else if ( day <= 244 && leapyear == YES ) {
                month = 8
                day = day - 213
        } else if ( day <= 273 && leapyear == NO ) {
                month = 9
                day = day - 243
        } else if ( day <= 274 && leapyear == YES ) {
                month = 9
                day = day - 244
        } else if ( day <= 304 && leapyear == NO ) {
                month = 10
                day = day - 273
        } else if ( day <= 305 && leapyear == YES ) {
                month = 10
                day = day - 274
	} else if ( day <=  334 && leapyear == NO ) {
                month = 11
                day = day - 304
        } else if ( day <= 335 && leapyear == YES ) {
                month = 11
                day = day - 305
	} else if ( leapyear == NO ) {
		month = 12
		day = day - 334
	} else {
		month = 12
		day = day - 335
	}



	date = 367*year-7*(year+(month+9)/12)/4+275*month/9+day+1721013.

	if ( year == 2000 && month>2 )
		date = date + 1

	jd = date + 0.5+00 + (hr + (min/60.))/24.
	

end






