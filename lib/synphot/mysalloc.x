# Debging versions of smark and sfree

procedure my_smark (sp)

pointer		sp
#--
int	ptop, itop
pointer	stack[1000]

data	ptop	/ 0 /

begin
	call smark (sp)

	ptop = ptop + 1
	stack[ptop] = sp
	return

	entry my_sfree (sp)

	for (itop = ptop; itop > 0; itop = itop - 1)
	    if (stack[itop] == sp)
		break

	if (itop == 0) {
	    call error (1, "Pointer not on stack")
	} else if (itop < ptop) {
	    call eprintf ("Pointer not at top of stack: itop = %d ptop = %d\n")
	    call pargi (itop)
	    call pargi (ptop)
	}

	ptop = itop - 1
	call sfree (sp)
end

