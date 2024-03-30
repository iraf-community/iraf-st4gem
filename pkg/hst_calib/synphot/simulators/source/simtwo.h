# SIMTWO.H -- Constants used in simimg and simspec

define	RA_UNITS	"hours"      # right ascension units
define	DEC_UNITS	"degrees"    # target declination units
define	DETECT_UNITS	"seconds"    # detector units (skycoord = NO)
define	RADIUS_UNITS	"seconds"    # extended object radius units
define	OBJ_UNITS	"counts"     # units used when computing effstim
define	APER_UNITS	"seconds"    # aperture pixel size units
define	BACK_UNITS	"seconds"    # units of earthight and thermal spectra
define	ZOD_UNITS	"degrees"    # units of zodiacal light
define	ZOD_BAND	"band(v)"    # magnitude passband of zodiacal light
define	APERKEY		"APERTURE"   # aperture keyword name

# Mathematical operations supported by fileops

define	ADDOP		1
define	SUBOP		2
define	MULOP		3
define	DIVOP		4
