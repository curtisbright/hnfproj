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
 * RNSOp.h includes rountines related to the operations in Residue Number 
 * System (RNS).
 *
 * Functions:
 *   - basisExt: given a representation of a matrix/vector in some RNS, 
 *     extend to compute the representation in another positive integer
 *
 *   - basisExtPos: given a representation of a non-negative matrix/vector in 
 *     some RNS, extend to compute the representation in another positive 
 *     integer
 *
 *   - basisProd: compute the product of elements of a RNS basis
 *
 *   - ChineseRemainder: given a representation of an integer in some RNS, use 
 *     Chinese Remainder Algorithm to reconstruct the integer
 * 
 *   - ChineseRemainderPos: given a representation of a non-negative integer 
 *     in some RNS, use Chinese Remainder Algorithm to reconstruct the integer
 *
 *   - combBasis: compute the special combination of RNS basis
 *
 *   - cumProd: compute the representation of the combination of elements of 
 *     one RNS basis in another RNS basis
 *
 *   - findRNS: find a RNS basis and its special combination
 * 
 *   - maxInter: compute the maximum interval of positive and negative results
 *      of a matrix-matrix or matrix-vector product
 * 
 *   - repBound: compute the mix radix coefficients of a special integer in a
 *     RNS basis
 *
 * Note: the modular operations in these functions are implemented by function
 *   fmod, which requires the bit-length of operands(e.g., RNS basis) could 
 *   fit into mantissa of floating point numbers, i.e., less than 53 bits.
 *
 */

#ifndef RNSOP_H
#define RNSOP_H 1

#include "cblas.h"
#include "gmp.h"
#include "basisop.h"
#include "common.h"
#include "iml.h"

void basisExt(const long len, const long n, const FiniteField p, \
	      const FiniteField *basis, const FiniteField *cmbasis, \
	      const double cumprod, const FiniteField *bdcoeff, Double **R, \
	      Double *RE);

void basisExtPos(const long len, const long n, const FiniteField p, \
		 const FiniteField *basis, const FiniteField *cmbasis, \
		 Double **R, Double *RE);

void basisProd(const long len, const FiniteField *q, mpz_t mp_prod);

void ChineseRemainder(const long len, const mpz_t mp_prod, \
		      const FiniteField *basis, const FiniteField *cmbasis, \
		      const FiniteField *bdcoeff, Double *Ac, mpz_t mp_Ac); 

void ChineseRemainderPos(const long len, const FiniteField *basis, \
			 const FiniteField *cmbasis, Double *Ac, mpz_t mp_Ac);

FiniteField *combBasis(const long basislen, const FiniteField *basis);

double *cumProd(const long basislen, const FiniteField *basis, \
		const long extbasislen, const FiniteField *extbasis);

FiniteField **findRNS(const FiniteField RNS_bound, const mpz_t mp_maxInter, \
		      long *length);

void maxInter(const mpz_t mp_prod, const mpz_t mp_alpha, const long n, \
	      mpz_t mp_b);

void maxExtInter (const mpz_t mp_alpha, const long n, mpz_t mp_b);

FiniteField *repBound(const long len, const FiniteField *basis, \
		      const FiniteField *cmbasis);

FiniteField RNSbound(const long n);

#endif /* !RNSOP_H */

