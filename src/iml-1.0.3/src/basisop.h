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
 *
 * BasisOp.h includes routines performing basic operations
 * 
 * Functions:
 *
 *   - Dmod: perform mod p operation inplace over a double matrix/submatrix
 *
 *   - DCopy: copy double matrix/submatrix to another matrix/submatrix
 *
 *   - randomDb: generate a random Double dense matrix
 *
 *   - RandPrime: generate a random FintieField prime
 *
 *   - maxMagnLong: compute maximum magnitude of a signed long matrix/submatrix
 *
 *   - maxMagnMP: compute maximum magnitude of a mpz_t integer matrix/submatrix
 *
 *   - scalCpMP: copy a submatrix with scale from one matrix to the other one
 *
 */

#ifndef BASISOP_H
#define BASISOP_H 1

#include "gmp.h"
#include "common.h"
#include "iml.h"

void Dmod (const Double p, Double *A, const long n, const long m, \
	   const long lda);

void DCopy (const long n, const long m, const Double* A, const long lda, \
	    Double* B, const long ldb);

FiniteField RandPrime (const FiniteField lb, const FiniteField hb);

long maxMagnLong (const long *A, const long n, const long m, const long lda);

void maxMagnMP (mpz_t *mp_A, const long n, const long m, const long lda, \
		mpz_t mp_max);

Double *randomDb (const long n, const long m, const long bd);

void scalCpMP (const long n, const long m, const long lda, const long ldm, \
	       const mpz_t mp_a, mpz_t *mp_A, mpz_t *mp_M);


#endif /* !BASISOP_H */
