#{ stplot package -- STSDAS vector graphics utilities.

package	stplot

task	igi,
	sgraph		= "stplot$x_stplot.e"
task    psikern         = "stplot$x_psikern.e"

# Plot attributes parameters pset
task    pltpar          = "stplot$pltpar.par"

# Axis attributes parameters pset
task    axispar         = "stplot$axispar.par"

# Device parameters pset (device name, append, and viewport)
task    dvpar           = "stplot$dvpar.par"

clbye()
