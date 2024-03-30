# SCALETRANS -- Scale the data for plotting in pltrans.

procedure scaletrans( npt, xval, yval, nrow, xdat, ydat, nfit, xfit, yfit, 
	              xform, yform, xmode, ymode, xmin, xmax, ymin, ymax )

int	npt		# i: number of points in xval, yval
real	xval[ARB]	# i: x-values of data
real	yval[ARB]	# i: y-values of data
int	nrow		# i: number of rows in xdat and ydat
real	xdat[ARB]	# i: x-values of data form disk (to be fit)
real	ydat[ARB]	# i: y-values ...
int	nfit		# i: number of points in fit
real	xfit[ARB]	# i: x-values of fit
real	yfit[ARB]	# i: y-values of fit
char	xform[ARB]	# i: form of x- data
char	yform[ARB]	# i: form of y-data
char	xmode[ARB]	# i: modestring for xdata
char	ymode[ARB]	# i: modestring for ydata
real	xmin, xmax, ymin, ymax # o: computed plot limits
#--
bool	xmag, ymag
char	mode1[SZ_LINE], mode2[SZ_LINE]
int	nmode
real	xxmin, xxmax, yymin, yymax, temp

int	strsearch()

begin
	# Set the magnitude flags

	xmag = strsearch(xform,"mag") > 0 || strsearch(xform,"MAG") > 0
	ymag = strsearch(yform,"mag") > 0 || strsearch(yform,"MAG") > 0

	# Save the user limits

	xxmin = xmin
	xxmax = xmax
	yymin = ymin
	yymax = ymax

	# Calculate plot limits from data

	call ylimit (xval, npt, xmin, xmax)
	call ylimit (yval, npt, ymin, ymax)

	call ylimit (xdat, nrow, xmin, xmax)
	call ylimit (ydat, nrow, ymin, ymax)

	call ylimit (xfit, nfit, xmin, xmax)
	call ylimit (yfit, nfit, ymin, ymax)

	# Reset the plot limits to the user limits where they are not INDEF

	if (! IS_INDEFR (xxmin))
	    xmin = xxmin

	if (! IS_INDEFR (xxmax))
	    xmax = xxmax

	if (! IS_INDEFR (yymin))
	    ymin = yymin

	if (! IS_INDEFR (yymax))
	    ymax = yymax

	# Invert limits if mag scale

	if (xmag) {
	   call splitmode(xmode, nmode, mode1, mode2)
	   if (nmode == 1) {
	      temp = xmin
	      xmin = xmax
	      xmax = temp
	   }
	}

	if (ymag) {
	    temp = ymin
	    ymin = ymax
	    ymax = temp
	}
	      
end
