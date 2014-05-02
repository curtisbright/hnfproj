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
 * nonsysolve.h includes routines solving the nonsingular system of linear 
 * equations.
 *
 * Functions:
 *
 *   - findliftbasisSmall: compute the p-adic lifting basis
 *
 *   - findliftbasisLarge: compute the p-adic lifting basis 
 *
 *   - adBasis: adjust the lifting basis if some element in the lifting basis 
 *     is bad
 *
 *   - liftbd: compute the initial number of lifting steps, initial solution 
 *     bounds, maximum number of lifting steps and maximum solution bounds
 *
 */


#ifndef NONSYSOLVE_H 
#define NONSYSOLVE_H 1

#include "gmp.h"
#include "basisop.h"
#include "common.h"
#include "padiclift.h"
#include "reconsolu.h"
#include "RNSop.h"
#include "iml.h"

void adBasis(const long idx, const long basislen, FiniteField *liftbasis);

FiniteField *findLiftbasisSmall(const long n, const mpz_t mp_alpha, \
				long *basislen);

FiniteField *findLiftbasisLarge(const long n, const mpz_t mp_alpha, \
				long *basislen);

void liftbd(const mpz_t mp_basisprod, const long n, const mpz_t mp_alpha, const mpz_t mp_beta, long *maxk, mpz_t mp_maxnb, mpz_t maxdb, long *k, mpz_t mp_nb, mpz_t mp_db);


#endif /* !NONSYSOLVE_H */





