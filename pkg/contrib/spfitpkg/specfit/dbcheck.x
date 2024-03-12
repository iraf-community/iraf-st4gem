###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	dbcheck
#
#  Description:	DBCHECK is an IRAF task for use with specfit that makes 
#		checking your specfit database files a little easier
#
#
#  Arguments:	database file specfier
#
#  Returns:	Errors about reading the file
#
#  Notes:	Information shared in common blocks defined in "specfit.com".
#
#  History:	August	1994	Created by jpg
#
###########################################################################

include "specfit.h"
include	<error.h>
include	<pkg/dttext.h>


procedure t_dbcheck ( )



char 	instring[SZ_CNAME,NCOMP]
int 	comtype[NCOMP]

double	stepa, par0a, tlima, blima, ptola, ifixa

char	fitnames[SZ_LINE], fitname[SZ_FNAME], actual_fit[SZ_FNAME]
char	temp[SZ_LINE], tempname[SZ_LINE], what[3]

int	i, j, rec, np, istat, nc , num, k
pointer	dt, list

int	dtlocate(), dtgeti(), fscan(), clgstr(), fntgfnb(), fntlenb(), strcmp()
pointer	dtmap1(), fntopnb()

int 	found,a_indef


errchk	 dtgeti(), dtgar(), dtlocate(), dtmap(), gargi(), gargr(), gargd()

include	"specfit.com"


begin
	#Initialize a list of the component types and how many associated params
	call sfinitcomp()


	nc = clgstr("database",tempname,SZ_LINE)

	what[1] = tempname[1]
	what[2] = tempname[2]
	what[3] = EOS

	# See if they have already added th sf prefix
	if ( strcmp(what,"sf") == 0 ) {
		call sprintf(fitnames,SZ_LINE,"%s")
			call pargstr(tempname)
	} else {
		call sprintf(fitnames,SZ_LINE,"sf%s")
		  call pargstr(tempname)
	}

	# Get a list of all the files that meet the file name descriptor

	list = fntopnb(fitnames,YES)
	num = fntlenb(list)

	if ( num == 0 ) {
		call printf("No files found to check with search string %s\n")
		call pargstr(fitnames)
	}

	# Go through all the files
	do k = 1, num {
		nc = fntgfnb(list,fitname,SZ_FNAME)

		call printf("\nChecking file %s\n")
		  call pargstr(fitname)


		iferr ( dt = dtmap1("./", fitname, READ_ONLY) ) {
			call printf("\tERROR: Can't open file\n")
			next;
		}

		do i = 3, SZ_FNAME {
			actual_fit[i-2] = fitname[i]
		}

		call sprintf(temp,SZ_LINE,"begin %s")
		  call pargstr(actual_fit)

		iferr ( rec = dtlocate (dt, actual_fit) ) {
			call printf("\tERROR: Problem in file: Can't find statement - %s !\n")
			  call pargstr(temp)
			next;
		}

		iferr ( ncomp = dtgeti(dt, rec, "components") ) {
			call printf("\tERROR: Problem in file\n")
			call printf("\t\tEither unable to locate keyword 'components' or\n")
			call printf("\t\tunable to find the integer representing the number\n")
			call printf("\t\tof components that should follow the keyword.\n")

			next;
		}




# Scan the database file to collect the components and decode their names
		for ( i = 0; i < ncomp; i = i + 1) {
			istat = fscan(DT(dt))
			    call gargwrd(instring[1,i+1], SZ_CNAME)

				found = NO
				do j = 1, NKEYS {
					if ( strcmp(instring[1,i+1],compkeys[1,j]) == 0 ) {
						found = YES
						comtype[i+1] = j
					}
				}
				
				if ( found == NO ) {

				#	call printf("%s\n")
				#		call pargstr(temp)
					if ( i!=0 && strcmp(instring[1,1],instring[1,i+1]) == 0 ) {
						call printf("\tERROR:  Beginning of parameter values reached before\n")
						call printf("\t\tdone reading component types\n")
					} else { 
						call printf("\tWARNING: component type `%s` not known\n")
						call pargstr(instring[1,i+1])
						comtype[i+1] = 500
					}
				}

			call sprintf(instring[1,i+1], SZ_CNAME, "%s%d")
			    call pargstr(instring[1,i+1])
			    call pargi(i+1)
		}



		for ( i = 0; i < ncomp; i = i + 1) {

			iferr ( np = dtgeti(dt, rec, instring[1,i+1]) ) {
				call printf("\tERROR: Can't find the word `%s`\n")
					call pargstr(instring[1,i+1])
				next;
			}

			if (comtype[i+1] != 500 && np != ncpar[comtype[i+1]] ) {
				call printf("\tWARNING: usually %d parameters not %d in component type %s\n")
					call pargi(ncpar[comtype[i+1]])
					call pargi(np)
					call pargstr(compkeys[1,comtype[i+1]])
			}

			for ( j = 1; j <= np; j = j + 1) {
				par0a = INDEFD
				blima = INDEFD
				tlima = INDEFD
				stepa = INDEFD
				ptola = INDEFD
				ifixa = INDEFD

				istat = fscan (DT(dt))
					call gargd(par0a)
					call gargd(blima)
					call gargd(tlima)
					call gargd(stepa)
					call gargd(ptola)
					call gargd(ifixa)


				if ( istat != OK ) 
					call printf("\tERROR: EOF before end of read data\n")
				else {
					a_indef = NO

					if ( IS_INDEFD(par0a) ) {
						call printf("\tERROR: Can't read Parameter value: Component %d, Parameter %d\n")
					  		call pargi(i+1)
					  		call pargi(j)
						a_indef = YES
					}

					if ( IS_INDEFD(blima) ) {
						call printf("\tERROR: Can't read Lower Limit value: Component %d, Parameter %d\n")
					  		call pargi(i+1)
					  		call pargi(j)
						a_indef = YES
					}


					if ( IS_INDEFD(tlima) ) {
						call printf("\tERROR: Can't read Upper Limit value: Component %d, Parameter %d\n")
					  		call pargi(i+1)
					  		call pargi(j)
						a_indef = YES
					}



					if ( IS_INDEFD(stepa) ) {
						call printf("\tERROR: Can't read Step Size value: Component %d, Parameter %d\n")
					  		call pargi(i+1)
					  		call pargi(j)
						a_indef = YES
					}


					if ( IS_INDEFD(ptola) ) {
						call printf("\tERROR: Can't read Tolerance value: Component %d, Parameter %d\n")
					  		call pargi(i+1)
					  		call pargi(j)
						a_indef = YES
					}


					if ( IS_INDEFD(ifixa) ) {
						call printf("\tERROR: Can't read Fix/Free value: Component %d, Parameter %d\n")
					  		call pargi(i+1)
					  		call pargi(j)
						a_indef = YES
					}

					if ( a_indef == NO ) {
						if ( par0a < blima ) {
							call printf("\tWARNING: Parameter is lower than its bottom limit\n")
							call printf("\t\tComponent %d, Parameter %d\n")
					  			call pargi(i+1)
					  			call pargi(j)
						}
						if ( par0a > tlima ) {
							call printf("\tWARNING: Parameter is larger than its upper limit\n")
							call printf("\t\tComponent %d, Parameter %d\n")
					  			call pargi(i+1)
					  			call pargi(j)
						}

						if ( ptola >= 1 ) {
							call printf("\tWARNING: Tolerance is very large\n")
							call printf("\t\tComponent %d, Parameter %d\n")
					  			call pargi(i+1)
					  			call pargi(j)
						}

						if ( ifixa < -1 || ifixa > ncomp) {
							call printf("\tERROR: Fix/Free has an impossible value\n")
							call printf("\t\tComponent %d, Parameter %d\n")
					  			call pargi(i+1)
					  			call pargi(j)
						}


					
					}

				}
			}

		}
	

		call dtunmap (dt)
		call printf("Done checking file %s\n")
		  call pargstr(fitname)

	}
	call fntclsb(list)

end








