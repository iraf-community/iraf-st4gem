# DEFLABEL -- Create default plot labels depending on input form

procedure deflabel( form, xlabel, ylabel, title)

char	form[ARB]	# i: Form of dependent variable
char	xlabel[ARB]	# o: X axis label
char	ylabel[ARB]	# o: Y axis label
char	title[ARB]	# o: Plot title

int	strsearch()

begin

	call strlwr(form)
	call strcpy("WAVELENGTH (A)", xlabel, SZ_LINE)
	call strcpy("SYNPHOT.PLSPEC", title, SZ_LINE)

	if ( strsearch(form, "flam") > 0)
	   call strcpy("FLAM (ergs cm-2 s-1 A-1)",ylabel, SZ_LINE)
	else if (strsearch(form, "fnu") > 0)
	   call strcpy("FNU (ergs cm-2 s-1 Hz-1)", ylabel, SZ_LINE)
	else if (strsearch(form, "mjy") > 0)
	   call strcpy("FNU (mJy)", ylabel, SZ_LINE)
	else if (strsearch(form, "jy") > 0)
	   call strcpy("FNU (Jy)", ylabel, SZ_LINE)
	else if (strsearch(form, "photlam") > 0)
	   call strcpy("NLAM (PHOTONS cm-2 s-1 A-1)", ylabel, SZ_LINE)
	else if (strsearch(form, "photnu") > 0)
	   call strcpy("NNU (PHOTONS cm-2 s-1 Hz-1)", ylabel, SZ_LINE)
	else if (strsearch(form, "counts") > 0 ||
	         strsearch(form, "photpix") > 0 )
	   call strcpy("DETECTED COUNTS (Pixel-1 s-1)", ylabel, SZ_LINE)
	else if ( strsearch(form, "abmag") > 0)
	   call strcpy("ABNU (MAGS)", ylabel, SZ_LINE)
	else if (strsearch(form, "stmag") > 0)
	   call strcpy("STLAM (MAGS)", ylabel, SZ_LINE)
	else {
	   call strupr(form)
	   call strcpy(form, ylabel, SZ_LINE)
	}

end
