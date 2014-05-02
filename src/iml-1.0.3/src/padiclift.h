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
 * padiclift.h includes routines performing p-adic lifting operations 
 *
 * Functions:
 *   - liftInit: perform initialization operations before lifting where the
 *     left hand side input matrix is a signed long matrix
 *
 *   - liftInitLlhs: perform initialization operations before lifting where 
 *     the left hand side input matrix is a mpz_t matrix
 *
 *   - liftInitRNS: perform initialization operations before lifting where the
 *     left hand side input matrix is represented in a RNS
 *
 *   - lift: compute p-adic lifting coefficients of system of linear equations 
 *
 */


#ifndef PADICLIFT_H 
#define PADICLIFT_H 1

#include "cblas.h"
#include "gmp.h"
#include "basisop.h"
#include "common.h"
#include "mtrans.h"
#include "RNSop.h"
#include "iml.h"

long liftInit(const long liftbasislen, const FiniteField *liftbasis, \
	      const long n, const long *A, mpz_t mp_basisprod, \
	      mpz_t mp_extbasisprod, long *extbasislen, FiniteField **cmbasis,\
	      FiniteField **extbdcoeff, Double **liftbasisInv, Double **AInv, \
	      FiniteField ***extbasis, Double ***ARNS);

long liftInitLlhs(const long liftbasislen, const FiniteField *liftbasis, \
		  const long n, mpz_t *mp_A, mpz_t mp_basisprod, \
		  mpz_t mp_extbasisprod, long *extbasislen, \
		  FiniteField **cmbasis, FiniteField **extbdcoeff, \
		  Double **liftbasisInv, Double **AInv, \
		  FiniteField ***extbasis, Double ***ARNS, double ** _Ainv, unsigned long * _X);

long liftInitRNS(const long liftbasislen, const FiniteField *liftbasis, \
		 const long basislen, const FiniteField *basis, const long n, \
		 Double **ARNS, mpz_t mp_liftbasisprod, mpz_t mp_extbasisprod,\
		 long *extbasislen, FiniteField **cmliftbasis, \
		 FiniteField **extbdcoeff, Double **liftbasisInv, \
		 Double **AInv, FiniteField ***extbasis, Double ***AExtRNS);


Double*** iml_lift(const enum SOLU_POS solupos, const long k, const long n, \
	       const long m, const long liftbasislen, const long extbasislen, \
	       const mpz_t mp_basisprod, const mpz_t mp_extbasisprod, \
	       const FiniteField *liftbasis, const FiniteField *cmbasis, \
	       const FiniteField *extbdcoeff, const Double *liftbasisInv,\
	       mpz_t *mp_r, FiniteField **extbasis, Double **AInv, \
	       Double **ARNS);


Double * invBasis(const long basislen, const FiniteField *basis, \
		  const mpz_t mp_basisprod);

#endif /* !PADICLIFT_H */
