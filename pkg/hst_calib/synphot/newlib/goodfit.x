define	FACSQ		0.8483037	# 1/(log10(e) * 2.5 )**2 

#* HISTORY *
#* B.Simon	13-Jul-94	rewritten from stats.x

# GOODFIT -- Compute goodness of fit statistics

procedure goodfit (data, daterr, theory, ndata, form, chisq, bias, rms, nstat)

real	data[ARB]	# i: measured data
real	daterr[ARB]	# i: error in data, used for weighting
real	theory[ARB]	# i: theoretical data
int	ndata		# i: number of data points
char	form[ARB]	# i: form of data
real	chisq		# o: chi squared
real	bias		# o: bias
real	rms		# o: root mean squared
int	nstat		# o: number of data points used in statistics
#--
int	mag, idat, nbad, ngood
real	dif, weight, badchi, badsum, badbias, badrms
real	goodchi, goodsum, goodbias, goodrms

int	is_magunit()

begin
	# Initialize sums to zero

	nbad = 0
	ngood = 0

	badchi = 0.0
	badsum = 0.0
	badbias = 0.0
	badrms = 0.0

	goodchi = 0.0
	goodsum = 0.0
	goodbias = 0.0
	goodrms = 0.0

	mag = is_magunit (form)

	# Compute two sets of sums, the first for points with no errors
	# (bad points), the second for points with errors (good points)

	do idat = 1, ndata {
	    if (IS_INDEFR(data[idat]) || IS_INDEFR(theory[idat]))
		next

	    if (IS_INDEFR(daterr[idat]) || daterr[idat] <= 0.0) {
		nbad = nbad + 1
		dif = data[idat] - theory[idat]
		badchi = badchi + dif * dif

		if (mag == NO) {
		    if (data[idat] <= 0.0 || theory[idat] <= 0.0)
			next
		    dif = -2.5 * alog10 (data[idat] / theory[idat])

		} else {
		    dif = data[idat] - theory[idat]
		}

		badsum = badsum + 1.0
		badbias = badbias + dif
		badrms = badrms + dif * dif

	    } else {
		ngood = ngood + 1
		dif = (data[idat] - theory[idat]) / daterr[idat]
		goodchi = goodchi + dif * dif

		if (mag == NO) {
		    if (data[idat] <= 0.0 || theory[idat] <= 0.0)
			next
		    dif = -2.5 * alog10 (data[idat] / theory[idat])
		    weight = FACSQ * (data[idat] / daterr[idat]) ** 2

		} else {
		    dif = data[idat] - theory[idat]
		    weight = 1.0 / (daterr[idat] * daterr[idat])
		}

		goodsum = goodsum + weight
		goodbias = goodbias + weight * dif
		goodrms = goodrms + weight * dif * dif
	    }
	}

	# If any points have errors, compute the statistics from these points
	# Otherwise compute statistics from points without errors
	# If no points at all, set statistics to INDEF

	if (goodsum > 0.0) {
	    chisq = goodchi / real (ngood)
	    bias = goodbias / goodsum
	    rms = sqrt (goodrms / goodsum)
	    nstat = ngood

	} else if (badsum > 0.0) {
	    chisq = badchi / real (nbad)
	    bias = badbias / badsum
	    rms = sqrt (badrms / badsum)
	    nstat = nbad

	} else {
	    chisq = INDEFR
	    bias = INDEFR
	    rms = INDEFR
	    nstat = INDEFI
	}

end
