# ADDMASKS.H -- Constants used by the addmasks task

include	<mach.h>

define	MAX_FLAGS	(NBITS_INT-1)

# Tokens from grammar used to define maskcompile

define	Y_WRONG		257
define	Y_LPAR		258
define	Y_RPAR		259
define	Y_IDENT		260
define	Y_IF		261
define	Y_DONE		262
define	Y_THEN		263
define	Y_ELSE		264
define	Y_PUSH		265
define	Y_CONST		266
define	Y_OR		267
define	Y_AND		268
define	Y_NOT		269
define	Y_EQ		270
define	Y_NE		271
define	Y_LT		272
define	Y_GT		273
define	Y_LE		274
define	Y_GE		275
