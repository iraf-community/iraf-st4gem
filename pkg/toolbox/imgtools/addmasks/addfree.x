# ADDFREE -- Free the memory used to store the pseudocode structure

procedure addfree (pcode)

pointer	pcode		# i: pseudocode structure
#--

begin
	call mfree (pcode, TY_INT)
end
