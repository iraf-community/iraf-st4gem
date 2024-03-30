# HST_CALIB.CL -- The STScI HST data analysis suite of packages
# Created: R.L. Williamson, 14-Jul-1993
#

procedure hst_calib()
string mode="al"

begin

	set synphot	= "hst_calib$synphot/"

	package hst_calib

	task	synphot.pkg	= "synphot$synphot.cl"

# Implicitly load synphot
synphot

	cl()

end
