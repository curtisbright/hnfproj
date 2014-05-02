/* ---------------------------------------------------------------------
 *
 * -- Integer Matrix Library (IML)
 *    (C) Copyright 2004, 2006 All Rights Reserved
 *
 * -- IML routines -- Version 1.0.1 -- November, 2006
 *
 * Author         : Zhuliang Chen
 * Contributor(s) : Arne Storjohann, Cory Fletcher
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
 *
 * certsolve.h includes routines certified solving system of linear equations
 * with input matrix in any shape.
 *
 * Functions:
 *   - revseq: perform operations on the permutation vector
 * 
 *   - compressBoundLong: compute the maximum magnitude of a compressed signed
 *      long submatrix
 *   
 *   - compressBoundMP: compute the maximum magnitude of a compressed mpz_t 
 *     submatrix
 *
 *   - LVecSMatMulCmp: verify the equality of a matrix-vector product and a
 *     scalar-vector product
 * 
 *   - specialminSolveLong: compute the minimal denominator solution of a full
 *     row rank system, represented by signed long integers, and corresponding
 *     certificate vector(optional). The solution size could be reduced.
 *
 *   - specialminSolveRNS: compute the minimal denominator solution of a full
 *     row rank system, represented in RNS, and corresponding certificate 
 *     vector(optional). The solution size could be reduced.
 *
 *   - certVerify: verify the correctness of a certificate vector for a system
 *     of linear equations
 *
 *   - minSolnoncompressLong: compute the minimal denominator solution of a 
 *     full row rank system without compression, represented by signed long 
 *     integers, and the corresponding certificate vector(optional). The 
 *     solution size could be reduced.
 *
 *   - minSolnoncompressRNS: compute the minimal denominator solution of a 
 *     full row rank system without compression, represented in RNS, and the
 *     corresponding certificate vector(optional). The solution size could be
 *     reduced.
 *
 *   - specialHermite: perform unimodular transformation inplace in matrix M
 *     and return an certificate matrix Cert(optional) at the same time
 *
 *   - mpz_init_array: allocate a mpz_t array length n dynamically and set all
 *     the entries in the array to be zero.
 *  
 *   - migcdex: compute solutions to the modulo N extended GCD problem
 *
 *   - kernelBasis: compute a kernel basis of a full row rank matrix
 *
 */


#ifndef CERTSOLVE_H
#define CERTSOLVE_H 1

#include "cblas.h"
#include "gmp.h"
#include "common.h"
#include "basisop.h"
#include "latreduce.h"
#include "padiclift.h"
#include "reconsolu.h"
#include "RNSop.h"
#include "iml.h"

#define _M(i,j)	   (M[((i)-1)*(k+t)+(j)-1]) 
#define _Cert(i,j) (Cert[((i)-1)*(n+k+t)+(j)-1])
#define _C(i,j)    (C[((i)-1)*(n-1)+(j)-1])

#define _P(i)    (P[(i)-1])
#define _Q(i)    (Q[(i)-1])
#define _P1(i)   (P1[(i)-1])
#define _Q1(i)   (Q1[(i)-1])
#define _c(i)    (c[(i)-1])
#define _tmpz(i) (tmpz[(i)-1])

long LVecSMatMulCmp (const enum MULT_POS mulpos, const long basislen, \
		     const long n, const long m, const FiniteField *basis, \
		     Double **ARNS, mpz_t mp_s, mpz_t *mp_V, mpz_t *mp_b);

long * revseq (const long r, const long m, const long *A);

void compressBoundLong (const long n, const long m, const long *Pt, \
			const long *A, mpz_t mp_bd);

void compressBoundMP (const long n, const long m, const long *Pt, \
		      mpz_t *mp_A, mpz_t mp_bd);


long specialminSolveLong (const long certflag, const long redflag, \
			  const long nullcol, const long n, const long m, \
			  const mpz_t mp_bdC, const long *C, mpz_t *mp_b, \
			  mpz_t *mp_N, mpz_t mp_D, mpz_t *mp_NZ, mpz_t mp_DZ);

long LongRNSbound(void);

long specialminSolveRNS (const long certflag, const long redflag, 
			 const long nullcol, const long n, const long m, \
			 const long basislen, const mpz_t mp_bdC, \
			 const FiniteField *basis, Double **CRNS, \
			 mpz_t *mp_b, mpz_t *mp_N, mpz_t mp_D,\
			 mpz_t *mp_NZ, mpz_t mp_DZ);

long certVerify (const long basislen, const long n, const long m, \
		 const FiniteField *basis, Double **ARNS, mpz_t mp_DZ, \
		 mpz_t *mp_NZ);

void minSolnoncompressLong (const long certflag, const long redflag, \
			    const long n, const long k, mpz_t *mp_Bb, \
			    const long* A,  mpz_t *mp_N, mpz_t mp_D, \
			    mpz_t *mp_NZ, mpz_t mp_DZ);

void minSolnoncompressRNS (const long certflag, const long redflag, \
			   const long n, const long k, const long basislen, \
			   const FiniteField *basis, mpz_t *mp_Bb, \
			   Double **ARNS, mpz_t *mp_N, mpz_t mp_D, \
			   mpz_t *mp_NZ, mpz_t mp_DZ);

void specialHermite (const long certflag, const long n, const long k, \
		     const long t, mpz_t *M, mpz_t *Cert);

mpz_t * mpz_init_array (const long n);

void migcdex (const mpz_t N, const mpz_t a, mpz_t *b, const long n, \
	      unsigned *c);

void kernelBasis (const long n, const long k, const long t, mpz_t *mp_M, mpz_t *mp_N);

#endif /* !CERTSOLV_H */


