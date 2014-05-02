/* ---------------------------------------------------------------------
 *
 * -- Integer Matrix Library (IML)
 *    (C) Copyright 2004, 2006 All Rights Reserved
 *
 * -- IML routines -- Version 1.0.1 -- November, 2006
 *
 * Author         : Zhuliang Chen
 * Contributor(s) : Arne Storjohann
 * University of Waterloo -- School of Computer Science
 * Waterloo, Ontario, N2L3G1 Canada
 *
 * ---------------------------------------------------------------------
 *
 * -- Copyright notice and Licensing terms:
 *
 *  Redistribution  and  use in  source and binary forms, with or without
 *  modification, are  permitted provided  that the following  conditions
 *  are met:
 *
 * 1. Redistributions  of  source  code  must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce  the above copyright
 *    notice,  this list of conditions, and the  following disclaimer in
 *    the documentation and/or other materials provided with the distri-
 *    bution.
 * 3. The name of the University,  the IML group,  or the names of its
 *    contributors  may not be used to endorse or promote products deri-
 *    ved from this software without specific written permission.
 *
 * -- Disclaimer:
 *
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,  INDIRECT, INCIDENTAL, SPE-
 * CIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO,  PROCUREMENT  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEO-
 * RY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  (IN-
 * CLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */



/*
 * Functions:
 *
 *   - kernelLong: computing the right kernel of an integer matrix 
 *     with type mpz_t and optionally reduce with LLL
 *
 *   - nullspaceMP: computing the right nullspace of an integer matrix 
 *     with type mpz_t 
 *
 *   - nullspaceLong: computing the right nullspace of an integer matrix 
 *     with type Long
 *
 *   nonsingular system solving:
 *
 *   - nonsingSolvMM: solve nonsingular system of linear equations, where the
 *     left hand side input matrix is a signed long matrix
 *
 *   - nonsingSolvLlhsMM: solve nonsingular system of linear equations, where
 *     the left hand side input matrix is a mpz_t matrix
 *
 *   - nonsingSolvRNSMM: solve nonsingular system of linear equations, where 
 *     the left hand side input matrix is represented in RNS
 *
 *   certified system solving:
 *
 *   - certSolveLong: certified solve a system of linear equations without 
 *     reducing the solution size, where the left hand side input matrix is 
 *     represented by signed long integers
 *   
 *   - certSolveRedLong: certified solve a system of linear equations and 
 *     reduce the solution size, where the left hand side input matrix is 
 *     represented by signed long integers
 *
 *   - certSolveMP: certified solve a system of linear equations without 
 *     reducing the solution size, where the left hand side input matrix is 
 *     represented by mpz_t integers
 *
 *   - certSolveRedMP: certified solve a system of linear equations andf
 *     reduce the solution size, where the left hand side input matrix is 
 *     represented by signed mpz_t integers
 *
 *   matrix operations in Z/p:
 *  
 *   - RowEchelonTransform: compute a mod p row-echelon transform of a mod p 
 *     input matrix
 *   
 *   - mAdjoint: compute the adjoint of a mod p square matrix
 *
 *   - mBasis: compute a basis for the rowspace and/or a basis for the left 
 *     nullspace of a mod p matrix
 *
 *   - mDeterminant: compute the determinant of a square mod p matrix
 *
 *   - mInverse: certified compute the inverse of a square mod p matrix inplace
 *
 *   - mRank: compute the rank of a mod p matrix
 *
 *   - mRankProfile: compute the rank profile of a mod p matrix
 *   
 *  Note: all the matrices in this package are treated as row majored matrices.
 */


#ifndef IML_H
#define IML_H 1

#ifdef __cplusplus
#  define BEGIN_C_DECLS    extern "C" {
#  define END_C_DECLS              }
#else
#  define BEGIN_C_DECLS
#  define END_C_DECLS
#endif

BEGIN_C_DECLS

typedef unsigned long FiniteField; 
typedef double Double;

#ifndef ENUM_DEFINED_H
enum SOLU_POS { LeftSolu=101, RightSolu=102 };
enum MULT_POS { LeftMul=201, RightMul=202 };
#define ENUM_DEFINED_H
#endif

extern long nullspaceMP(const long n, 
 		        const long m, 
			const mpz_t *A, 
			mpz_t * *mp_N_pass);


extern long nullspaceLong(const long n, 
			  const long m, 
			  const long *A, 
			  mpz_t * *mp_N_pass);


extern long kernelLong(const long n, 
		       const long m, 
		       const long *A, 
		       mpz_t * *mp_N_pass,
		       const long reduce);


extern void nonsingSolvMM (const enum SOLU_POS solupos, 
			   const long n, 
			   const long m, 
			   const long *A, 
			   mpz_t *mp_B, 
			   mpz_t *mp_N, 
			   mpz_t mp_D);

extern void nonsingSolvLlhsMM_X (const enum SOLU_POS solupos, 
			       const long n, 
			       const long m, 
			       mpz_t *mp_A, 
			       mpz_t *mp_B, 
			       mpz_t *mp_N, 
			       mpz_t mp_D, double ** _Ainv, unsigned long * _X);
extern void nonsingSolvLlhsMM (const enum SOLU_POS solupos, 
			       const long n, 
			       const long m, 
			       mpz_t *mp_A, 
			       mpz_t *mp_B, 
			       mpz_t *mp_N, 
			       mpz_t mp_D);


extern void nonsingSolvRNSMM (const enum SOLU_POS solupos, 
			      const long n, 
			      const long m, 
			      const long basislen, 
			      const FiniteField *basis,
			      Double **ARNS, 
			      mpz_t *mp_B,
			      mpz_t *mp_N,
			      mpz_t mp_D);


extern long certSolveLong (const long certflag, 
			   const long n,
			   const long m,
			   const long *A, 
			   mpz_t *mp_b, 
			   mpz_t *mp_N,
			   mpz_t mp_D, 
			   mpz_t *mp_NZ,
			   mpz_t mp_DZ);


extern long certSolveRedLong (const long certflag, 
			      const long nullcol, 
			      const long n,
			      const long m,
			      const long *A,
			      mpz_t *mp_b,
			      mpz_t *mp_N,
			      mpz_t mp_D,
			      mpz_t *mp_NZ,
			      mpz_t mp_DZ);


extern long certSolveMP (const long certflag, 
			 const long n,
			 const long m,
			 mpz_t *mp_A,
			 mpz_t *mp_b,
			 mpz_t *mp_N, 
			 mpz_t mp_D,
			 mpz_t *mp_NZ,
			 mpz_t mp_DZ);


extern long certSolveRedMP (const long certflag,
			    const long nullcol,
			    const long n, 
			    const long m, 
			    mpz_t *mp_A,
			    mpz_t *mp_b, 
			    mpz_t *mp_N, 
			    mpz_t mp_D,
			    mpz_t *mp_NZ,
			    mpz_t mp_DZ);


extern void RowEchelonTransform (const FiniteField p,
				 Double *A, 
				 const long n, 
				 const long m, 
				 const long frows,
				 const long lrows, 
				 const long redflag, 
				 const long eterm, 
				 long *Q, 
				 long *rp,
				 FiniteField *d);


extern Double * mAdjoint (const FiniteField p, 
			  Double *A, 
			  const long n);


extern long mBasis (const FiniteField p, 
		    Double *A, 
		    const long n, 
		    const long m, 
		    const long basis, 
		    const long nullsp, 
		    Double **B, 
		    Double **N);


extern long mDeterminant (const FiniteField p, 
			  Double *A, 
			  const long n);


extern long mInverse (const FiniteField p,
		      Double *A,
		      const long n);


extern long mRank (const FiniteField p, 
		   Double *A,
		   const long n, 
		   const long m);


extern long * mRankProfile (const FiniteField p, 
			    Double *A,
			    const long n,
			    const long m);



END_C_DECLS

#endif /* !IML_H */
