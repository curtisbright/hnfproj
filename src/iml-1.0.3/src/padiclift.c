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


#include "padiclift.h"

/*
 * Calling Sequence:
 *   -1/i <-- liftInit(liftbasislen, liftbasis, n, A, mp_basisprod, 
 *            mp_extbasisprod, extbasislen, cmbasis, extbdcoeff, 
 *            liftbasisInv, AInv, extbasis, ARNS)
 *
 * Summary:
 *   Perform initialization operations before lifting, where the left hand 
 *   side input matrix is a signed long matrix
 *
 * Description:
 *   Initialization step before lifting computes necessory data and stores
 *   them to avoid recomputations if continuous lifting is needed. Seperating
 *   initialization step and lifting step makes it possible to lift adaptively.
 *   Pointers are passed into the calling sequence to store the outputs of 
 *   initialization.
 *
 * Input:
 *   liftbasislen: long, dimension of lifting basis
 *      liftbasis: 1-dim FiniteField array length liftbasislen, lifting basis
 *              n: long, dimension of input matrix A
 *              A: 1-dim long array length n*n, representation of n x n input
 *                 matrix
 *
 * Return:
 *   The first index i such that A^(-1) mod liftbasis[i] doesn't exist, where 
 *   i starts from 0. Otherwise, return -1 if A^(-1) mod liftbasis[i] exists 
 *   for any i
 *
 * Output: 
 *     mp_basisiprod: mpz_t, product of lifting basis
 *   mp_extbasisprod: mpz_t, product of extended RNS basis 
 *       extbasislen: pointer to long int, storing the dimension of extended 
 *                    RNS basis. Extended basis is used to compute the 
 *                    matrix-vector product AC_i during lifting
 *           cmbasis: pointer to a 1-dim FiniteField array, storing the 
 *                    special combination of lifting basis computed by 
 *                    function combBasis
 *        extbdcoeff: pointer to a 1-dim FiniteField array, storing the mix 
 *                    radix coefficients of a special integer, computed by 
 *                    function repBound, in extended RNS basis 
 *      liftbasisInv: pointer to a 1-dim Double array, storing 
 *                    (1/mp_basisprod) mod extbasis[i]
 *              AInv: pointer to a 1-dim Double array, storing the modp matrix
 *                    A^(-1) mod liftbasis[i] 
 *          extbasis: pointer to a 2-dim FiniteField array, where
 *                  - (*extbasis)[0] = extended RNS basis
 *                  - (*extbasis)[1] = the special combination of extended RNS
 *                    basis computed by function combBasis
 *              ARNS: pointer to a 2-dim Double array, where the last 
 *                    dimension (*ARNS)[i] stores A mod ith element in 
 *                    extended RNS basis
 */

long
liftInit (const long liftbasislen, const FiniteField *liftbasis, \
	  const long n, const long *A, mpz_t mp_basisprod, \
	  mpz_t mp_extbasisprod, long *extbasislen, FiniteField **cmbasis, \
	  FiniteField **extbdcoeff, Double **liftbasisInv, Double **AInv, \
	  FiniteField ***extbasis, Double ***ARNS)
{
  long i, j, alpha, minv, p, temp, len=0;
  mpz_t mp_maxInter, mp_alpha;
  FiniteField *q, *qinv;

  for (i = 0; i < liftbasislen; i++)
    {
      p = (long)liftbasis[i]; 
      for (j = 0; j < n*n; j++) 
	AInv[i][j] = (double)((temp = (A[j] % p)) >= 0 \
			      ? temp : p+temp);
      minv = mInverse(liftbasis[i], AInv[i], n);

      /* if fail to find inverse of A mod liftbasis[i] */
      if (minv == 0) { return i; }
    }
  *cmbasis = combBasis(liftbasislen, liftbasis);
  basisProd(liftbasislen, liftbasis, mp_basisprod);

  /* compute maximum intermediate result mp_maxInter */
  alpha = maxMagnLong(A, n, n, n);
  mpz_init_set_ui(mp_alpha, alpha);
  mpz_init(mp_maxInter);
  maxExtInter(mp_alpha, n, mp_maxInter);
  mpz_clear(mp_alpha);

  *extbasis = findRNS(liftbasis[liftbasislen-1]-1, mp_maxInter, &len);
  mpz_clear(mp_maxInter);
  *extbasislen = len;
  q = *(*extbasis);
  qinv = *((*extbasis)+1);
  *liftbasisInv = invBasis(len, q, mp_basisprod);
  basisProd(len, q, mp_extbasisprod);
  *extbdcoeff = repBound(len, q, qinv);
  *ARNS = XMALLOC(Double *, len);
  for (i = 0; i < len; i++)
    {
      p = (long)q[i];
      (*ARNS)[i] = XMALLOC(Double, n*n);
      for (j = 0; j < n*n; j++) 
	(*ARNS)[i][j] = (double)((temp = (A[j] % p)) >= 0 ? \
				 temp : p+temp);
    }
  return -1;
}




/*
 * Calling Sequence:
 *   -1/i <-- liftInitLlhs(liftbasislen, liftbasis, n, mp_A, mp_basisprod, 
 *                         mp_extbasisprod, extbasislen, cmbasis, extbdcoeff, 
 *                         liftbasisInv, AInv, extbasis, ARNS)
 *
 * Summary:
 *   Perform initialization operations before lifting, where the left hand 
 *   side input matrix is a mpz_t matrix
 *
 * Description:
 *   Initialization step before lifting computes necessory data and stores
 *   them to avoid recomputations if continuous lifting is needed. Seperating
 *   initialization step and lifting step makes it possible to lift adaptively.
 *   Pointers are passed into the calling sequence to store the outputs of 
 *   initialization.
 *
 * Input:
 *   liftbasislen: long, dimension of lifting basis
 *      liftbasis: 1-dim FiniteField array length liftbasislen, lifting basis
 *              n: long, dimension of input matrix A
 *           mp_A: 1-dim mpz_t array length n*n, representation of n x n input
 *                 matrix
 *
 * Return:
 *   The first index i such that A^(-1) mod liftbasis[i] doesn't exist, where 
 *   i starts from 0. Otherwise, return -1 if A^(-1) mod liftbasis[i] exists 
 *   for any i
 *
 * Output: 
 *     mp_basisiprod: mpz_t, product of lifting basis
 *   mp_extbasisprod: mpz_t, product of extended RNS basis 
 *       extbasislen: pointer to long int, storing the dimension of extended 
 *                    RNS basis. Extended basis is used to compute the 
 *                    matrix-vector product AC_i during lifting
 *           cmbasis: pointer to a 1-dim FiniteField array, storing the 
 *                    special combination of lifting basis computed by 
 *                    function combBasis
 *        extbdcoeff: pointer to a 1-dim FiniteField array, storing the mix 
 *                    radix coefficients of a special integer, computed by 
 *                    function repBound, in extended RNS basis 
 *      liftbasisInv: pointer to a 1-dim Double array, storing 
 *                    (1/mp_basisprod) mod extbasis[i]
 *              AInv: pointer to a 1-dim Double array, storing the modp matrix
 *                    A^(-1) mod liftbasis[i] 
 *          extbasis: pointer to a 2-dim FiniteField array, where
 *                  - (*extbasis)[0] = extended RNS basis
 *                  - (*extbasis)[1] = the special combination of extended RNS
 *                    basis computed by function combBasis
 *              ARNS: pointer to a 2-dim Double array, where the last 
 *                    dimension (*ARNS)[i] stores A mod ith element in 
 *                    extended RNS basis
 */

static long
myInverse (const FiniteField p, Double ** A, double ** _Ainv, unsigned long * _X)
{
  long i;
  if(_X && _Ainv) {
    for(i = 0; _X[i]; ++i) {
      if (p == _X[i]) {
        XFREE(*A);
        *A = _Ainv[i];
        return 1;
      }
    }
  }
  return 0;
}

long
liftInitLlhs (const long liftbasislen, const FiniteField *liftbasis, \
	      const long n, mpz_t *mp_A, mpz_t mp_basisprod, \
	      mpz_t mp_extbasisprod, long *extbasislen, \
	      FiniteField **cmbasis, FiniteField **extbdcoeff, \
	      Double **liftbasisInv, Double **AInv, FiniteField ***extbasis, \
	      Double ***ARNS, double ** _Ainv, unsigned long * _X)
{
  long i, j, minv, len=0;
  mpz_t mp_maxInter, mp_alpha;
  FiniteField *q, *qinv;

  for (i = 0; i < liftbasislen; i++)
    {
      minv = myInverse(liftbasis[i], &(AInv[i]), _Ainv, _X);
      if (minv == 1) { continue; }
      for (j = 0; j < n*n; j++) {
        AInv[i][j] = (Double)mpz_fdiv_ui(mp_A[j], liftbasis[i]);
      }
      minv = mInverse(liftbasis[i], AInv[i], n);

      /* if fail to find inverse of A mod liftbasis[i] */
      if (minv == 0) { return i; }
    }
  *cmbasis = combBasis(liftbasislen, liftbasis);
  basisProd(liftbasislen, liftbasis, mp_basisprod);

  /* compute maximum intermediate result mp_maxInter */
  mpz_init(mp_alpha);
  maxMagnMP(mp_A, n, n, n, mp_alpha);
  mpz_init(mp_maxInter);
  maxExtInter(mp_alpha, n, mp_maxInter);
  mpz_clear(mp_alpha);

  *extbasis = findRNS(liftbasis[liftbasislen-1]-1, mp_maxInter, &len);
  mpz_clear(mp_maxInter);
  *extbasislen = len;
  q = *(*extbasis);
  qinv = *((*extbasis)+1);
  *liftbasisInv = invBasis(len, q, mp_basisprod);
  basisProd(len, q, mp_extbasisprod);
  *extbdcoeff = repBound(len, q, qinv);
  *ARNS = XMALLOC(Double *, len);
  for (i = 0; i < len; i++)
    {
      (*ARNS)[i] = XMALLOC(Double, n*n);
      for (j = 0; j < n*n; j++) 
	(*ARNS)[i][j] = (Double)mpz_fdiv_ui(mp_A[j], q[i]);
    }

  return -1;
}


/*
 * Calling Sequence:
 * -1/i <-- liftInitRNS(liftbasislen, liftbasis, basislen, basis, n, ARNS, 
 *          mp_liftbasisprod, mp_extbasisprod, extbasislen, cmliftbasis, 
 *	    liftbasisInv, extbdcoeff, AInv, extbasis, AExtRNS)
 *
 * Summary:
 *   Perform initialization operations before lifting where the left hand 
 *   side input matrix is represented in a RNS
 *   
 * Description:
 *   Initialization step before lifting computes necessory data and stores
 *   them to avoid recomputations if continuous lifting is needed. Seperating
 *   initialization step and lifting step makes it possible to lift adaptively.
 *   Pointers are passed into the calling sequence to store the outputs of 
 *   initialization.
 *
 * Input:
 *   liftbasislen: long, dimension of lifting basis
 *      liftbasis: 1-dim FiniteField array length liftbasislen, lifting basis
 *       basislen: long, dimension of RNS basis used to represent A
 *          basis: 1-dim FiniteField array length basislen, RNS basis used to
 *                 represent A
 *              n: long, dimension of A
 *           ARNS: 2-dim Double array, dimension basislen x n^2,
 *                 representation of A in RNS, ARNS[i] = A mod basis[i]
 *
 * Return:
 *   The first index i such that A^(-1) mod liftbasis[i] doesn't exist, where 
 *   i starts from 0. Otherwise, return -1 if A^(-1) mod liftbasis[i] exists 
 *   for any i
 *
 * Output: 
 *   mp_liftbasisiprod: mpz_t, product of lifting basis
 *     mp_extbasisprod: mpz_t, product of extended RNS basis 
 *         extbasislen: pointer to long int, storing the dimension of extended 
 *                      RNS basis. Extended basis is used to compute the 
 *                      matrix-vector product AC_i during lifting
 *             cmbasis: pointer to a 1-dim FiniteField array, storing the 
 *                      special combination of lifting basis computed by 
 *                      function combBasis
 *          extbdcoeff: pointer to a 1-dim FiniteField array, storing the mix 
 *                      radix coefficients of a special integer, computed by 
 *                      function repBound, in extended RNS basis 
 *      liftbasisInv: pointer to a 1-dim Double array, storing 
 *                    (1/mp_basisprod) mod extbasis[i]
 *                AInv: pointer to a 1-dim Double array, storing the modp 
 *                      matrix A^(-1) mod liftbasis[i] 
 *            extbasis: pointer to a 2-dim FiniteField array, where
 *                    - (*extbasis)[0] = extended RNS basis
 *                    - (*extbasis)[1] = the special combination of extended 
 *                      RNS basis computed by function combBasis
 *                ARNS: pointer to a 2-dim Double array, where the last 
 *                      dimension (*ARNS)[i] stores A mod ith element in 
 *                      extended RNS basis
 */

long
liftInitRNS (const long liftbasislen, const FiniteField *liftbasis, \
	     const long basislen, const FiniteField *basis, const long n, \
	     Double **ARNS, mpz_t mp_liftbasisprod, mpz_t mp_extbasisprod, \
	     long *extbasislen, FiniteField **cmliftbasis, \
	     FiniteField **extbdcoeff, Double **liftbasisInv, Double **AInv, \
	     FiniteField ***extbasis, Double ***AExtRNS)
{
  long i, minv, len=0;
  double *cumprodRNS;
  mpz_t mp_maxInter, mp_alpha;
  FiniteField *q, *qinv, *cmbasis, *bdcoeff;

  cmbasis = combBasis(basislen, basis);
  bdcoeff = repBound(basislen, basis, cmbasis);
  cumprodRNS = cumProd(basislen, basis, liftbasislen, liftbasis);
  for (i = 0; i < liftbasislen; i++)
    {
      /* compute A mod liftbasis[i] from ARNS */
      basisExt(basislen, n*n, liftbasis[i], basis, cmbasis, cumprodRNS[i], \
	       bdcoeff, ARNS, AInv[i]);
      minv = mInverse(liftbasis[i], AInv[i], n);

      /* if fail to find inverse of A mod basis[i] */
      if (minv == 0) 
	{ XFREE(bdcoeff); XFREE(cmbasis); XFREE(cumprodRNS); return i; }
    }
  XFREE(cumprodRNS);
  *cmliftbasis = combBasis(liftbasislen, liftbasis);
  basisProd(liftbasislen, liftbasis, mp_liftbasisprod);

  /* compute maximum intermediate result mp_maxInter */
  mpz_init(mp_alpha);
  basisProd(basislen, basis, mp_alpha);
  mpz_init(mp_maxInter);
  maxExtInter(mp_alpha, n, mp_maxInter);
  mpz_clear(mp_alpha);

  *extbasis = findRNS(liftbasis[liftbasislen-1]-1, mp_maxInter, &len);
  mpz_clear(mp_maxInter);
  *extbasislen = len;
  q = *(*extbasis);
  qinv = *((*extbasis)+1);
  *liftbasisInv = invBasis(len, q, mp_liftbasisprod);
  basisProd(len, q, mp_extbasisprod);
  *extbdcoeff = repBound(len, q, qinv);
  *AExtRNS = XMALLOC(Double *, len);
  cumprodRNS = cumProd(basislen, basis, len, q);
  for (i = 0; i < len; i++)
    {
      (*AExtRNS)[i] = XMALLOC(Double, n*n);

      /* compute A mod extbasis[i] from ARNS */
      basisExt(basislen, n*n, q[i], basis, cmbasis, cumprodRNS[i], \
	       bdcoeff, ARNS, (*AExtRNS)[i]);
    }

  { XFREE(bdcoeff); XFREE(cmbasis); XFREE(cumprodRNS); }

  return -1;
}



/*
 * Calling Sequence:
 *   C <-- lift(solupos, k, n, m, liftbasislen, extbasislen, mp_basisprod, 
 *              mp_extbasisprod, liftbasis, cmbasis, extbdcoeff, 
 *              liftbasiInv, mp_r, extbasis, AInv, ARNS)
 *
 * Summary:
 *   Compute p-adic lifting coefficients of system of linear equations 
 *
 * Description:
 *   Given a system of linear equations AX = mp_r or Transpose(A)X = mp_r, 
 *   where A is a n x n nonsingular matrix and mp_r is a n x m matrix, the 
 *   function computes and stores lifting coefficients upto k lifting steps. 
 *
 *   The data computed from initialization function are reused each time this
 *   function is called. The right hand side matrix mp_r is updated such that
 *   we could continue to call this function to perform lifting using updated
 *   mp_r if we do not lift high enough.
 *
 * Input:
 *           solupos: enumerate, flag to indicate whether to transpose A or not
 *                  - solupos = LeftSolu: system be Transpose(A)X = mp_r
 *                  - solupos = RightSolu: system be AX = mp_r
 *                 k: long, number of lifting steps
 *                 n: long, dimension of A
 *                 m: long, column dimension of right hand side matrix mp_r
 *      liftbasislen: long, dimension of lifting basis
 *       extbasislen: long, dimension of extended RNS basis
 *     mp_basisiprod: mpz_t, product of lifting basis
 *   mp_extbasisprod: mpz_t, product of extended lifting basis
 *         liftbasis: 1-dim FiniteField array length liftbasislen, storing
 *                    lifting basis
 *           cmbasis: 1-dim FiniteField array length liftbasislen, storing the 
 *                    special combination of lifting basis computed in
 *                    initialization step
 *        extbdcoeff: 1-dim FiniteField array length liftbasislen, storing the
 *                    mix radix coefficients of a special integer in extended 
 *                    RNS basis computed in initialization step
 *      liftbasisInv: a 1-dim Double array, storing 
 *                    (1/mp_basisprod) mod extbasis[i]
 *              mp_r: 1-dim mpz_t array length n*m, representation of n x m 
 *                    right hand side lifting matrix
 *          extbasis: 2-dim FiniteField array, dimension 2 x extbasislen, 
 *                    computed in initialization step
 *              AInv: 1-dim Double array length liftbasislen, storing the modp 
 *                    matrix A^(-1) mod liftbasis[i] 
 *              ARNS: 2-dim Double array, dimension extbasislen x n^2, where 
 *                    the second dimension (*ARNS)[i] stores A mod ith element
 *                    in extended RNS basis
 *
 * Output:
 *   C: 3-dim Double array, dimension k x liftbasislen x n*m 
 *    - If solupos = RightSolu, then C[i][j] represents the n x m coefficient 
 *      matrix computed by A^(-1)mp_r mod liftbasis[j] at the ith lifting step.
 *    - If solupos = LeftSolu, then C[i][j] represents the n x m coefficient
 *      matrix computed by Transpose(A)^(-1)mp_r mod liftbasis[j] at the ith
 *      lifting step.
 *
 * Precondition: 
 *   Any element p in array liftbasis must satisfy n*(p-1)^2 <= 2^53-1. 
 */

Double ***
iml_lift (const enum SOLU_POS solupos, const long k, const long n, \
      const long m, const long liftbasislen, const long extbasislen, \
      const mpz_t mp_basisprod, const mpz_t mp_extbasisprod, \
      const FiniteField *liftbasis, const FiniteField *cmbasis, \
      const FiniteField *extbdcoeff, const Double *liftbasisInv, \
      mpz_t *mp_r, FiniteField **extbasis, Double **AInv, Double **ARNS)
{
  long i, j, l;
  mpz_t mp_r1;
  FiniteField *q, *qinv;
  Double *dtemp, *dtemp1, **Ac, ***C;

  /* initialize lifting coefficient matrix C[k][liftbasislen][n] */
  C = XMALLOC(Double **, k);
  for (i = 0; i < k; i++) 
    {
      C[i] = XMALLOC(Double *, liftbasislen);
      for (j = 0; j < liftbasislen; j++) 
	C[i][j] = XMALLOC(Double, m*n);
    }
  q = *extbasis;
  qinv = *(extbasis+1);
  mpz_init(mp_r1);
  Ac = XCALLOC(Double *, extbasislen);
  for (i = 0; i < extbasislen; i++) 
    Ac[i] = XCALLOC(Double, m*n);
  dtemp = XMALLOC(Double, m*n);
  dtemp1 = XMALLOC(Double, extbasislen);

  /* start lifting */
  for (i = 0; i < k; i++)
    {
      /* compute coefficients of p-adic lifting C[i][l] */
      for (l = 0; l < liftbasislen; l++)
	{
	  /* mod(mp_r, liftbasis[j]) */
	  for (j = 0; j < m*n; j++) 
	    dtemp[j] = (Double)mpz_fdiv_ui(mp_r[j], liftbasis[l]);

	  /* compute the coefficients of p-adic lifting */
	  if (solupos == LeftSolu)
	    {
	      if (m == 1)
		cblas_dgemv(CblasRowMajor, CblasTrans, n, n, 1.0, AInv[l], \
			    n, dtemp, 1, 0.0, C[i][l], 1);
	      else
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, \
			    n, 1.0, AInv[l], n, dtemp, m, 0.0, C[i][l], m);
	    }
	  else if (solupos == RightSolu)
	    {
	      if (m == 1)
		cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, AInv[l],\
			    n, dtemp, 1, 0.0, C[i][l], 1);
	      else
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, m, \
			    n, 1.0, AInv[l], n, dtemp, m, 0.0, C[i][l], m);
	    }
	  Dmod((double)liftbasis[l], C[i][l], n, m, m);
	}

      /* compute Ac mod extbasis[j] */
      for (j = 0; j < extbasislen; j++)
	{
	  basisExtPos(liftbasislen, m*n, q[j], liftbasis, cmbasis, C[i], \
		      dtemp);
	  if (solupos == LeftSolu)
	    {
	      if (m == 1)
		cblas_dgemv(CblasRowMajor, CblasTrans, n, n, 1.0, ARNS[j],\
			    n, dtemp, 1, 0.0, Ac[j], 1);
	      else
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, n, \
			    1.0, ARNS[j], n, dtemp, m, 0.0, Ac[j], m);
	    }
	  else if (solupos == RightSolu)
	    {
	      if (m == 1)
		cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, ARNS[j],\
			    n, dtemp, 1, 0.0, Ac[j], 1);
	      else
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, m,\
			    n, 1.0, ARNS[j], n, dtemp, m, 0.0, Ac[j], m);
	    }
	  Dmod((double)q[j], Ac[j], n, m, m);
	}

      /* compute r_quo_p+(r mod p-Ac)/p */
      for (j = 0; j < m*n; j++)
	{
	  /* mp_r[j] := Quo(mp_r[j], p), mp_r1 := Mod(mp_r[j], p) */
	  mpz_fdiv_qr(mp_r[j], mp_r1, mp_r[j], mp_basisprod);

	  /* compute ((r mod p) mod q[l] - Ac mod q[l])(1/p mod q[l]) */
	  for (l = 0; l < extbasislen; l++)
	    {
	      dtemp1[l] = (Double)mpz_fdiv_ui(mp_r1, q[l]);
	      dtemp1[l] = fmod(dtemp1[l]+(q[l]-1)*Ac[l][j], q[l]);
	      dtemp1[l] = fmod(dtemp1[l]*liftbasisInv[l], q[l]);
	    }

	  /* compute (r mod p-Ac)(1/p) by CRT */
	  ChineseRemainder(extbasislen, mp_extbasisprod, q, qinv, extbdcoeff, \
			   dtemp1, mp_r1);
	  mpz_add(mp_r[j], mp_r[j], mp_r1);
	}
    }

  mpz_clear(mp_r1); 
  { XFREE(dtemp); XFREE(dtemp1); }
  for (i = 0; i < extbasislen; i++) { XFREE(Ac[i]); } { XFREE(Ac); }

  return C;
}
  

/* return (1/mp_basisprod) mod basis[i] */

Double *
invBasis(const long basislen, const FiniteField *basis, \
	 const mpz_t mp_basisprod)
{
  long i;
  mpz_t mp_temp, mp_basis;
  Double *inv;

  { mpz_init(mp_temp); mpz_init(mp_basis); }
  inv = XMALLOC(Double, basislen);
  for (i = 0; i < basislen; i++)
    {
      mpz_set_ui(mp_basis, basis[i]);
      mpz_invert(mp_temp, mp_basisprod, mp_basis);
      inv[i] = mpz_get_d(mp_temp);
    }
  { mpz_clear(mp_temp); mpz_clear(mp_basis); }
  return inv;
}
