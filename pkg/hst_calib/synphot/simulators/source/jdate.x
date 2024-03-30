#* HISTORY *
#* B.Simon	14-Apr-95	original

# JDATE -- Convert date to Julian date

double procedure jdate (year, month, day, hour, minute, sec)

int	year		# i: year in common era
int	month		# i: month number
int	day		# i: day in month
int	hour		# i: hour number
int	minute		# i: minute 
int	sec		# i: second
#--
double	date
int	yr, mo, a, b

begin
	# Julian date routine from Jean Meeus, "Astronomical Algorithms", 
	# chapter 7

	# Month and day  to day of year formula assumes year starts in March 

	yr = year
	mo = month

	if (mo < 3) {
	    yr = yr - 1
	    mo = mo + 12
	}

	# Compute century leap days

	a = yr / 100
	b = (a / 4) + 2 - a

	# Julian date formula

	date = int (365.25 * (yr + 4716)) + int (30.6001 * (mo + 1)) + 
	       day + b - 1524.5

	# Add fractional day

	date = date + (sec + 60 * (minute + 60 * hour)) / 86400.0
	return (date)
end

