#{ ST4GEM.CL -- STSDAS Tasks needed by the GEMINI package.
#

cl < "st4gem$/lib/zzsetenv.def"

package	st4gem, bin = st4gembin$

task	analysis.pkg	= "analysis$analysis.cl"
task	graphics.pkg	= "graphics$graphics.cl"
task	toolbox.pkg	= "toolbox$toolbox.cl"


# Load some of the ST4GEM packages
#cl < st4gem.cl

# Print the Welcome banner
if (motd)
type st4gem$st4gem_motd
;
clbye ()
