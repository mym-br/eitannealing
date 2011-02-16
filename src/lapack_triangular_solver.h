/*
 * lapack_triangular_solver.h
 *
 *  Created on: Jul 28, 2010
 *      Author: thiago
 */

#ifndef _LAPACK_TRIANGULAR_SOLVER_
#define _LAPACK_TRIANGULAR_SOLVER_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*
 * calls into Fortran
 */
void dtptrs_(
	char UPLO,			/*	(input) CHARACTER*1
							= 'U':  A is upper triangular;
							= 'L':  A is lower triangular.
						*/
	char TRANS, 		/*	(input) CHARACTER*1
							Specifies the form of the system of equations:
							= 'N':  A * X = B  (No transpose)
							= 'T':  A**T * X = B  (Transpose)
							= 'C':  A**H * X = B  (Conjugate transpose = Transpose)
						*/
	char DIAG, 			/*	(input) CHARACTER*1
							= 'N':  A is non-unit triangular;
							= 'U':  A is unit triangular
						*/
	int N, 				/*	(input) INTEGER
							The order of the matrix A.  N >= 0.
						*/
	int NRHS, 			/*	(input) INTEGER
							The number of right hand sides, i.e., the number of columns
							of the matrix B.  NRHS >= 0.
						*/
	const double *AP, 	/*	 (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
							 The upper or lower triangular matrix A, packed columnwise in
							 a linear array.  The j-th column of A is stored in the array
							 AP as follows:
							 if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
							 if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
						 */


	double *B,
	int LDB,
	int *INFO)


#ifdef __cplusplus
}		// extern "C"
#endif // __cplusplus


#endif	/* _LAPACK_TRIANGULAR_SOLVER_ */
