
#############################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#                        
#	Synopsis	procedure mat_inverse(nfree,alpha,beta)
#				int nfree	number of variables
#				real alpha[nfree,nfree] covariance matrix
#							to invert
#				real beta[nfree]	
#
#	Description:	Uses a Gauss-Jordan Method to invert a matrix.  
#
#
#	Arguments:  A is the Matrix to Invert
#
#
#
#	History:	Dec 1994	J Grimes 
#
#############################################################

include "specfit.h"

procedure mat_inverse(nfree,alpha, beta)
int 	nfree
real 	alpha[nfree,nfree]
real 	beta[nfree]

real 	temp, high, invert
int 	i, j , k, row, column

real 	Abs()

pointer	pivotptr,rowptr,columnptr

begin
	call malloc(pivotptr,nfree,TY_INT)
	call malloc(rowptr,nfree,TY_INT)
	call malloc(columnptr,nfree,TY_INT)
	

	do j= 1, nfree {
		Memi[pivotptr+j-1] = 0
	}

	do i = 1, nfree {
		high = 0
		do j = 1, nfree {
			if ( Memi[pivotptr+j-1] != 1 ) {
				do k = 1, nfree {
					if ( Memi[pivotptr+k-1] > 1 ) {
						call eprintf("Matrix Inversion Error - Singular Matrix (Check Input Parameters)!!\n")
					}
					else if ( Memi[pivotptr+k-1] == 0 ) {
						if ( Abs(alpha[j,k]) >= high ){
							high = Abs(alpha[j,k])
							row = j
							column = k
						}
					}
				}	 
			}
		}
		
		Memi[pivotptr+column-1] = Memi[pivotptr+column-1] + 1
		if ( row != column ) {
			temp = beta[row]
			beta[row] = beta[column]
			beta[column] = temp

			do k = 1, nfree {
				temp = alpha[row,k]
				alpha[row,k] = alpha[column,k]
				alpha[column,k] = temp
			}
		}


		Memi[columnptr+i-1] = column
		Memi[rowptr+i-1] = row
		
		if ( alpha[column,column] == 0 )
			call eprintf("Matrix Inversion Error - Singular Matrix (Check Input Parameters)!!\n")


		invert = 1.0 / alpha[column,column]
		alpha[column,column] = 1
	
		do k = 1, nfree {
			alpha[column,k] = alpha[column,k] * invert
		}

		beta[column] = beta[column] * invert
	

		do k = 1, nfree {
			if ( k != column ) {
				temp = alpha[k, column]
				beta[k] = beta[k] - beta[column] * temp
				alpha[k,column] = 0
				do j = 1, nfree {
					alpha[k,j]=alpha[k,j]-alpha[column,j]*temp
				}
			}
		}
	} 

	do j = nfree, 1, -1 {
		if ( Memi[rowptr+j-1] != Memi[columnptr+j-1] ) {
			do k = 1, nfree {
				temp = alpha[k, Memi[rowptr+j-1]]
				alpha[k,Memi[rowptr+j-1]] = alpha[k, Memi[columnptr+j-1]]
				alpha[k,Memi[columnptr+j-1]] = temp
			}
		}
		
	}

	call mfree(rowptr,TY_INT)
	call mfree(columnptr,TY_INT)
	call mfree(pivotptr,TY_INT)
end


