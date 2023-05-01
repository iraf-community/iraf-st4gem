# Global MKPKG definitions for the ST4GEM packages.

$set XFLAGS	= "$(XFLAGS) -p st4gem"
$set XVFLAGS	= "$(XVFLAGS) -p st4gem"
$set LFLAGS	= "$(LFLAGS) -p st4gem"


$set USE_CDBS	    = yes
$set USE_CONTRIB    = yes	# Compile and link contrib package?
$set USE_PLAYPEN    = yes
$set USE_SUNVIEW   = yes

# Set up the NO_UNDERSCORE macro that can be used by C programs.
$ifeq (MACH, linux64)
    $include "st4gem$lib/mkpkg.sf.LNUX64"
$endif
