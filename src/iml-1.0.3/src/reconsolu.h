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
 * reconsolu.c includes routines performing rational reconstructions
 *
 * Functions:
 *   - sumliftCoeff: recursively sum up the p-adic lifting coefficients
 *
 *   - find2exp: compute floor(log[2](len))
 *
 *   - iratrecon: perform rational number reconstruction
 *
 *   - soluRecon: try reconstructing the rational solution using p-adic 
 *     lifting coefficients
 *
 *   - findNumer: certify correctness of the input denominator and, upon 
 *     success, compute the numerator
 *
 */

#ifndef RECONSOLU_H
#define RECONSOLU_H 1

#include "gmp.h"
#include "common.h"
#include "iml.h"
#include "RNSop.h"

void sumliftCoeff(const mpz_t mp_basisprod, const long k, mpz_t* C, \
		  mpz_t mp_sum);

void sumCoeff_rec(long start, long len, mpz_t *C, mpz_t *mp_pow, long splflag,\
		  long savflag, long *idx, mpz_t *mp_left, mpz_t mp_right);

long find2exp(const long len);

long iratrecon(const mpz_t mp_u, const mpz_t mp_m, const mpz_t mp_nb, \
	       const mpz_t mp_db, mpz_t mp_N, mpz_t mp_D);

long soluRecon(const enum SOLU_POS solupos, const long k, const long basislen,\
	       const long n, const long m, const mpz_t mp_basisprod, \
	       const FiniteField *basis, const FiniteField *cmbasis, \
	       Double ***C, mpz_t mp_nb, mpz_t mp_db, mpz_t *mp_N, mpz_t mp_D);

long findNumer(const mpz_t mp_u, const mpz_t mp_m, const mpz_t mp_D, \
	       const mpz_t mp_nb, mpz_t mp_N);


#endif /* !RECONSOLU_H */
