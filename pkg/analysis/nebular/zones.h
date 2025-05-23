include <iraf.h>
#  zones.h --	Definition of zones data structure		25-Aug-2008

#  Symbolics
define	N_ZONES		3		# Number of ionization zones
define	LOW		1		# index of low I.P. zone
define	MED		2		# index of medium I.P. zone
define	HI		3		# index of high I.P. zone

#  Structure for nebular stratification
define	Z_ATL		Memi[$1]	 # Atomic data object list
define	Ne_Ciii		Memr[P2R($1+1)]	 # Density for med ionization region
define	Ne_Oii		Memr[P2R($1+2)]	 # Density for low ionization region
#define	Ne_NEii		Memr[P2R($1+3)]	 # Density for low ionization region
define	Ne_NEiv		Memr[P2R($1+4)]	 # Density high low ionization region
define	Ne_Sii		Memr[P2R($1+5)]	 # Density for low ionization region
define	Ne_CLiii	Memr[P2R($1+6)]	 # Density med low ionization region
#define	Ne_ARii		Memr[P2R($1+7)]	 # Density for low ionization region
define	Ne_ARiv		Memr[P2R($1+8)]	 # Density for med ionization region
define	Te_Nii		Memr[P2R($1+9)]	 # Temp. for low ionization region
define	Te_Oii		Memr[P2R($1+10)] # Temp. for low ionization region
define	Te_Oiii		Memr[P2R($1+11)] # Temp. for med ionization region
define	Te_NEiii	Memr[P2R($1+12)] # Temp. for med ionization region
#define	Te_NEiv		Memr[P2R($1+13)] # Temp. for high ionization region
define	Te_NEv		Memr[P2R($1+14)] # Temp. for high ionization region
define	Te_Sii		Memr[P2R($1+15)] # Temp. for low ionization region
define	Te_Siii		Memr[P2R($1+16)] # Temp. for med ionization region
define	Te_ARiii	Memr[P2R($1+17)] # Temp. for med ionization region
define	Te_ARiv		Memr[P2R($1+18)] # Temp. for med ionization region
define	Te_ARv		Memr[P2R($1+19)] # Temp. for high ionization region

define	Ne_Low		Memr[P2R($1+20)] # Density for low ionization region
define	Ne_Med		Memr[P2R($1+21)] # Density for medium ionization region
define	Ne_Hi		Memr[P2R($1+22)] # Density for high ionization region
define	Te_Low		Memr[P2R($1+23)] # Temp. for low ionization region
define	Te_Med		Memr[P2R($1+24)] # Temp. for medium ionization region
define	Te_Hi		Memr[P2R($1+25)] # Temp. for high ionization region
define	LEN_ZONES	26		# size of ZONES structure

