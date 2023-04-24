# ADDCOMPILE.COM -- Global variables used by addmasks compiler

pointer	line		# Buffer containing next line in expression
pointer	ch		# Pointer to next character in expression
int	ncode		# Length of code array
pointer	code		# Pointer to next available code

common	/addcom/  line, ch, ncode, code
