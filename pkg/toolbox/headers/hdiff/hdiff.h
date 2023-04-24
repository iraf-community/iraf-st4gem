# HDIFF.H -- The size and location of information in the user area

# Note that hdiff does not use the imio$db interface to access header
# keywords. If the structure of the user area should ever change, hdiff
# would have to be rewritten. The current structure is that the user area 
# is divided into fixed length records, each terminated by a newline. Each
# record has the form "keyword = value / comment", where the keyword is 
# located in the first eight characters and the value starts in or after 
# the tenth character.

define	SZ_KEYWORD	8
define	SZ_RECORD	81
define	SZ_VALUE	70

define	START_VALUE	10
