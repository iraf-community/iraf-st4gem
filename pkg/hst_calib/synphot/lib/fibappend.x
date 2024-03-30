include	"../graferr.h"

define	EOM	(EOS-1)

# FIB_APPENDC -- Append a character to a dynamically allocated array

procedure fib_appendc (val, buffer, index)

char	val		#  i: Character to append to array
pointer	buffer		# io: Pointer to memory buffer
int	index		# io: Array index at which to place character
#--
errchk	fib_malloc, fib_realloc

begin
	if (buffer == NULL) {
	    index = 1
	    call fib_malloc (buffer, 1, TY_CHAR)
	}

	if (Memc[buffer+index-1] == EOM)
	    call fib_realloc (buffer, index, TY_CHAR)

	Memc[buffer+index-1] = val
end
 
# FIB_APPENDI -- Append an integer to a dynamically allocated array

procedure fib_appendi (val, buffer, index)

int	val		#  i: Value to append to array
pointer	buffer		# io: Pointer to memory buffer
int	index		# io: Array index at which to place value
#--
errchk	fib_malloc, fib_realloc
                        
begin
	if (buffer == NULL) {
	    index = 1
	    call fib_malloc (buffer, 1, TY_INT)
	}

	if (Memi[buffer+index-1] == EOM)
	    call fib_realloc (buffer, index, TY_INT)

	Memi[buffer+index-1] = val
end
 
# FIB_MALLOC -- Get an oversized block of memory for a buffer

procedure fib_malloc (buffer, length, dtype)

pointer	buffer		# o: Pointer to the memory buffer
int	length		# i: Amount of memory needed
int	dtype		# i: Data type of the buffer
#--
int	fib1, fib2, block_size
pointer	ptr

string	fiberr  "Invalid data type for memory allocation by FIB_MALLOC"

errchk	malloc

begin
	# Find the first Fibonacci number that is larger than the size of
	# of the memory block requested. This is the size of the memory
	# block that will be allocated.

	fib1 = 1
	fib2 = 2
	repeat {
	    block_size = fib1 + fib2
	    fib1 = fib2
	    fib2 = block_size
	} until (block_size > length)

	call malloc (buffer, block_size, dtype)

	# Set markers which indicate end of user length and allocated length

	switch (dtype) {
	case TY_CHAR:
	    for (ptr = buffer; ptr < buffer + block_size - 1; ptr = ptr + 1)
		Memc[ptr] = EOS
	    Memc[ptr] = EOM
	case TY_SHORT:
	    for (ptr = buffer; ptr < buffer + block_size - 1; ptr = ptr + 1)
		Mems[ptr] = EOS
	    Mems[ptr] = EOM
	case TY_INT:
	    for (ptr = buffer; ptr < buffer + block_size - 1; ptr = ptr + 1)
		Memi[ptr] = EOS
	    Memi[ptr] = EOM
	case TY_LONG:
	    for (ptr = buffer; ptr < buffer + block_size - 1; ptr = ptr + 1)
		Meml[ptr] = EOS
	    Meml[ptr] = EOM
	default:
	    call error (INTERNAL, fiberr)
	}

end

# FIB_REALLOC -- Reallocate an oversized block of memory

procedure fib_realloc (buffer, length, dtype)

pointer	buffer		# io: Pointer to the memory buffer
int	length		#  i: Current length of memory block
int	dtype		#  i: Data type of the buffer
#--
int	fib1, fib2, block_size
pointer	ptr

string	fiberr  "Invalid data type for memory allocation by FIB_REALLOC"

errchk	realloc

begin
	# Find the first Fibonacci number that is larger than the size of
	# of the memory block requested. This is the size of the memory
	# block that will be allocated.

	fib1 = 1
	fib2 = 2
	repeat {
	    block_size = fib1 + fib2
	    fib1 = fib2
	    fib2 = block_size
	} until (block_size > length)

	call realloc (buffer, block_size, dtype)

	# Set markers which indicate end of user length and allocated length

	switch (dtype) {
	case TY_CHAR:
	    for (ptr = buffer + length - 1;
		 ptr < buffer + block_size - 1;
		 ptr = ptr + 1)
		    Memc[ptr] = EOS
	    Memc[ptr] = EOM
	case TY_SHORT:
	    for (ptr = buffer + length - 1;
		 ptr < buffer + block_size - 1;
		 ptr = ptr + 1)
 		    Mems[ptr] = EOS
	    Mems[ptr] = EOM
	case TY_INT:
	    for (ptr = buffer + length - 1;
		 ptr < buffer + block_size - 1;
		 ptr = ptr + 1)
		    Memi[ptr] = EOS
	    Memi[ptr] = EOM
	case TY_LONG:
	    for (ptr = buffer + length - 1;
		 ptr < buffer + block_size - 1;
		 ptr = ptr + 1)
		    Meml[ptr] = EOS
	    Meml[ptr] = EOM
	default:
	    call error (INTERNAL, fiberr)
	}

end
