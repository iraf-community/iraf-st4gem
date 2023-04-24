# Global MKPKG definitions for the ST4GEM packages.

$set XFLAGS	= "$(XFLAGS) -p st4gem"
$set XVFLAGS	= "$(XVFLAGS) -p st4gem"
$set LFLAGS	= "$(LFLAGS) -p st4gem"


$set USE_CDBS	    = yes
$set USE_CONTRIB    = yes	# Compile and link contrib package?
$set USE_PLAYPEN    = yes
$set USE_SUNVIEW   = yes

# Set up the NO_UNDERSCORE macro that can be used by C programs.
$ifeq (MACH, hp700, rs6000, vms, _vax, _alpha)
    $set XFLAGS = "$(XFLAGS) -DNO_UNDERSCORE"
$else $ifeq (MACH, macosx)
    $set XFLAGS = "$(XFLAGS) -DMACOSX"
$endif
