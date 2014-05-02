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


#include "RNSop.h"

/*
 *
 * Calling Sequence:
 *   basisExt(len, n, p, basis, cmbasis, cumprod, bdcoeff, R, RE)
 *
 * Summary:
 *   Given a representation of a matrix/vector in some RNS, extend to compute 
 *   the representation in another positive integer
 *
 * Description:
 *   Let R be the representation of a matrix/vector M in a residue basis
 *   'basis', i.e., R[i] = mod(M, basis[i]) (i = 0..len-1). The function 
 *   computes the representation of M in another positive integer p, 
 *   RE = mod(M, p) using Garner's algorithm. 'mod' represents positive modular
 *   operation.
 *
 *   Let q be product of all elements in the RNS basis. Every entry m in M 
 *   satisfies -(q-1)/2 <= m <= (q-1)/2. That is, M has both positive entries
 *   and negative entries.
 *
 *   To avoid repeat computations, the function takes some precomputed 
 *   informations as input, which are listed below.
 *
 * Input: 
 *       len: long, dimension of RNS basis
 *         n: long, length of array RE
 *         p: FiniteField, modulus
 *     basis: 1-dim FiniteField array length len, RNS basis
 *   cmbasis: 1-dim FiniteField array length len, computed by function 
 *            combBasis, inverses of special combination of RNS basis
 *   cumprod: double, computed by function cumProd, 
 *            (-basis[0]*basis[1]*...*basis[len-1]) mod p 
 *   bdcoeff: 1-dim FiniteField array length len, computed by function repBound
 *         R: 1-dim Double array length n*len, representation of a len x n 
 *            matrix, R[i]=mod(M, basis[i]) (i=0..len-1)
 * 
 * Output:
 *   RE: 1-dim Double array length n, RE = mod(M, p), the space of RE 
 *       should be allocated before calling the function.
 *
 * Precondition:
 *   t <= 2^53-1, where t is the maximal intermidiate value arised in this 
 *   function,
 *   t = max(2*(basis[len-1]-1)^2, (p-1)^2+basis[len-1]-1)
 *           
 */

void
basisExt (const long len, const long n, const FiniteField p, \
	  const FiniteField *basis, const FiniteField *cmbasis, \
	  const double cumprod, const FiniteField *bdcoeff, Double **R, \
	  Double *RE)
{
  long i, j;
  double temp;
  Double **U;
  const FiniteField *q, *qinv;

  q = basis;
  qinv = cmbasis;

  /* if p = q[i] then just copy the corresponding column to RE */
  for (i = 0; i < len; i++) 
    { if (p == q[i]) { cblas_dcopy(n, R[i], 1, RE, 1); return; } }
  U = XMALLOC(Double *, len);
  for (i = 0; i < len; i++) { U[i] = XMALLOC(Double, n); }

  /* compute the coefficients of mix radix in positive representation by 
     inplacing modular matrix U */
  cblas_dcopy(n, R[0], 1, U[0], 1);
  for (i = 1; i < len; i++)
    {
      cblas_dcopy(n, U[i-1], 1, U[i], 1);
      for (j = i-2; j >= 0; j--)
	{
	  cblas_dscal(n, (double)(q[j] % q[i]), U[i], 1);
	  cblas_daxpy(n, 1.0, U[j], 1, U[i], 1);
	  Dmod((double)q[i], U[i], 1, n, n);
	}
      temp = (double)qinv[i]*(double)(q[i]-1);
      temp = fmod(temp, (double)q[i]);
      cblas_dscal(n, temp, U[i], 1);
      cblas_daxpy(n, (double)qinv[i], R[i], 1, U[i], 1);
      Dmod((double)q[i], U[i], 1, n, n);
    }

  /* compute mod(r, p) in positive representation and store into RE */
  cblas_dcopy(n, U[len-1], 1, RE, 1);
  Dmod((double)p, RE, 1, n, n);
  for (i = len-2; i >= 0; i--)
    {
      cblas_dscal(n, (double)(q[i] % p), RE, 1);
      cblas_daxpy(n, 1.0, U[i], 1,  RE, 1);
      Dmod((double)p, RE, 1, n, n);
    }

  /* convert to symmetric representation */
  for (i = 0; i < n; i++)
    for (j = len-1; j >= 0; j--)
      {
	if (U[j][i] > bdcoeff[j])
	  {
	    RE[i] = fmod(RE[i]+cumprod, (double)p);
	    break;
	  }
	else if (U[j][i] < bdcoeff[j]) { break; }
      }
  for (i = 0; i < len; i++) { XFREE(U[i]); } { XFREE(U); }

  return;
}


/*
 *
 * Calling Sequence:
 *   basisExtPos(len, n, p, basis, cmbasis, R, RE)
 *
 * Summary:
 *   Given a representation of a non-negative matrix/vector in some RNS, 
 *   extend to compute the representation in another positive integer
 *
 * Description:
 *   Let R be the representation of a matrix/vector M in a residue basis
 *   'basis', i.e., R[i] = mod(M, basis[i]) (i = 0..len-1). The function 
 *   computes the representation of M in another positive integer p, 
 *   RE = mod(M, p) using Garner's algorithm. 'mod' represents positive modular
 *   operation.
 *
 *   Let q be product of all elements in the RNS basis. Every entry m in M 
 *   satisfies 0 <= m <= q-1. That is, M only contains non-negative entries.
 *
 *   To avoid repeat computations, the function takes some precomputed 
 *   informations as input, which are listed below.
 *
 * Input:
 *       len: long, dimension of RNS basis
 *         n: long, length of array RE
 *         p: FiniteField, modulus
 *     basis: 1-dim FiniteField array length len, RNS basis
 *   cmbasis: 1-dim FiniteField array length len, computed by function 
 *            combBasis, inverses of special combination of RNS basis
 *         R: 1-dim Double array length n*len, representation of a len x n 
 *            matrix, R[i]=mod(M, basis[i]) (i=0..len-1)
 * 
 * Output:
 *   RE: 1-dim Double array length n, RE = mod(M, p), the space of RE 
 *       should be allocated before calling the function.
 *
 * Precondition: 
 *   t <= 2^53-1, where t is the maximal intermidiate value arised in this 
 *   function, t = max(2*(basis[len-1]-1)^2, (p-1)^2+basis[len-1]-1)
 *           
 */

void
basisExtPos (const long len, const long n, const FiniteField p, \
	     const FiniteField *basis, const FiniteField *cmbasis, \
	     Double **R, Double *RE)
{
  long i, j;
  double temp;
  Double **U;
  const FiniteField *q, *qinv;

  q = basis;
  qinv = cmbasis;

  /* if p = q[i] then just copy the corresponding column to RE */
  for (i = 0; i < len; i++) 
    { if (p == q[i]) { cblas_dcopy(n, R[i], 1, RE, 1); return; } }
  U = XMALLOC(Double *, len);
  for (i = 0; i < len; i++) { U[i] = XMALLOC(Double, n); }

  /* compute the coefficients of mix radix in positive representation by 
     inplacing modular matrix U */
  cblas_dcopy(n, R[0], 1, U[0], 1);
  for (i = 1; i < len; i++)
    {
      cblas_dcopy(n, U[i-1], 1, U[i], 1);
      for (j = i-2; j >= 0; j--)
	{
	  cblas_dscal(n, (double)(q[j] % q[i]), U[i], 1);
	  cblas_daxpy(n, 1.0, U[j], 1, U[i], 1);
	  Dmod((double)q[i], U[i], 1, n, n);
	}
      temp = (double)qinv[i]*(double)(q[i]-1);
      temp = fmod(temp, (double)q[i]);
      cblas_dscal(n, temp, U[i], 1);
      cblas_daxpy(n, (double)qinv[i], R[i], 1, U[i], 1);
      Dmod((double)q[i], U[i], 1, n, n);
    }

  /* compute mod(r, p) in positive representation and store into RE */
  cblas_dcopy(n, U[len-1], 1, RE, 1);
  Dmod((double)p, RE, 1, n, n);
  for (i = len-2; i >= 0; i--)
    {
      cblas_dscal(n, (double)(q[i] % p), RE, 1);
      cblas_daxpy(n, 1.0, U[i], 1, RE, 1);
      Dmod((double)p, RE, 1, n, n);
    }
  for (i = 0; i < len; i++) { XFREE(U[i]); } { XFREE(U); }

  return;
}



/*
 *
 * Calling Sequence:
 *   basisProd(len, basis, mp_prod)
 *
 * Summary:
 *   Compute the product of elements of a RNS basis
 *
 * Description:
 *   Let a RNS basis be 'basis'. The function computes the product of its 
 *   elements basis[0]*basis[1]*...*basis[len-1].
 *
 * Input:
 *     len: long, dimension of RNS basis
 *   basis: 1-dim FiniteField array length len, RNS basis
 * 
 * Output:
 *   mp_prod: mpz_t, product of elements in 'basis'
 *
 */

void
basisProd (const long len, const FiniteField *basis, mpz_t mp_prod)
{
  long i;

  mpz_set_ui(mp_prod, basis[0]);
  for (i = 1; i < len; i++) { mpz_mul_ui(mp_prod, mp_prod, basis[i]); }
}



/*
 *
 * Calling Sequence:
 *   ChineseRemainder(len, mp_prod, basis, cmbasis, bdcoeff, Ac, mp_Ac)
 *
 * Summary:
 *   Given a representation of an integer in some RNS, use Chinese Remainder 
 *   Algorithm to reconstruct the integer
 * 
 * Description:
 *   Let A be an integer, and Ac contains the representation of A in a RNS
 *   basis 'basis', i.e. Ac[i] = mod(A, basis[i]), (i = 0..len). Here 'mod' 
 *   is in positive representation. The function reconstructs the integer A 
 *   given the RNS basis 'basis' and Ac. 
 *
 *   To avoid repeat computations, the function takes some precomputed 
 *   informations as input, which are listed below.
 *
 * Input: 
 *       len: long, dimension of RNS basis
 *   mp_prod: mpz_t, computed by function basisProd, product of RNS basis
 *     basis: 1-dim FiniteField array length len, RNS basis
 *   cmbasis: 1-dim FiniteField array length len, computed by function 
 *            combBasis, inverses of special combination of RNS basis
 *   bdcoeff: 1-dim FiniteField array length len, computed by function repBound
 *        Ac: 1-dim Double array length n, representation of A in RNS 
 *
 * Output:
 *   mp_Ac: mpz_t, reconstructed integer A
 *
 * Precondition:
 *   Let q be product of all elements in the RNS basis. Then A must satisfy 
 *   -(q-1)/2 <= A <= (q-1)/2. 
 *
 */

void
ChineseRemainder (const long len, const mpz_t mp_prod, \
		  const FiniteField *basis, const FiniteField *cmbasis, \
		  const FiniteField *bdcoeff, Double *Ac, mpz_t mp_Ac) 
{
  long i, j;
  double temp, tempq, tempqinv;
  Double *U;

  U = XMALLOC(Double, len);

  /* compute the coefficients of mix radix in positive representation by 
     inplacing modular matrix U */
  U[0] = Ac[0];
  for (i = 1; i < len; i++)
    {
      U[i] = U[i-1];
      tempq = (double)basis[i];
      tempqinv = (double)cmbasis[i];
      for (j = i-2; j >= 0; j--)
	{
	  U[i] = U[j] + U[i]*fmod((double)basis[j], tempq);
	  U[i] = fmod(U[i], tempq);
	}
      temp = fmod(tempqinv*(double)(basis[i]-1), tempq);
      U[i] = fmod(tempqinv*Ac[i]+temp*U[i], tempq);
    }
  /* compute Ac in positive representation */
  mpz_set_d(mp_Ac, U[len-1]); 
  for (j = len-2; j >= 0; j--)
    {
      mpz_mul_ui(mp_Ac, mp_Ac, basis[j]);
      mpz_add_ui(mp_Ac, mp_Ac, (FiniteField)U[j]);
    }
  /* transfer from positive representation to symmetric representation */
  for (j = len-1; j >= 0; j--)
    {
      if (U[j] > bdcoeff[j])
	{
	  mpz_sub(mp_Ac, mp_Ac, mp_prod);
	  break;
	}
      else if (U[j] < bdcoeff[j]) { break; }
    }
  XFREE(U);

  return;
}



/*
 *
 * Calling Sequence:
 *   ChineseRemainderPos(len, basis, cmbasis, Ac, mp_Ac)
 *
 * Summary:
 *   Given a representation of a non-negative integer in some RNS, use Chinese 
 *   Remainder Algorithm to reconstruct the integer
 *
 * Description:
 *   Let A be a non-negative integer, and Ac contains the representation of A
 *   in a RNS basis 'basis', i.e. Ac[i] = mod(A, basis[i]), (i = 0..len). 
 *   Here 'mod' is in positive representation. The function reconstructs the 
 *   integer A given the RNS basis 'basis' and Ac. 
 *
 *   To avoid repeat computations, the function takes some precomputed 
 *   informations as input, which are listed below.
 *
 * Input: 
 *       len: long, dimension of RNS basis
 *     basis: 1-dim FiniteField array length len, RNS basis
 *   cmbasis: 1-dim FiniteField array length len, computed by function 
 *            combBasis, inverses of special combination of RNS basis
 *        Ac: 1-dim Double array length n, representation of A in RNS 
 *
 * Output:
 *   mp_Ac: mpz_t, reconstructed integer A
 *
 * Precondition:
 *   Let q be product of all elements in the RNS basis. Then A must satisfy 
 *   0 <= A <= q-1. 
 *
 */

void
ChineseRemainderPos (const long len, const FiniteField *basis, \
		     const FiniteField *cmbasis, Double *Ac, mpz_t mp_Ac) 
{
  long i, j;
  double temp, tempq, tempqinv;
  Double *U;

  U = XMALLOC(Double, len);

  /* compute the coefficients of mix radix in positive representation by 
     inplacing modular matrix U */
  U[0] = Ac[0];
  for (i = 1; i < len; i++)
    {
      U[i] = U[i-1];
      tempq = (double)basis[i];
      tempqinv = (double)cmbasis[i];
      for (j = i-2; j >= 0; j--)
	{
	  U[i] = U[j] + U[i]*fmod((double)basis[j], tempq);
	  U[i] = fmod(U[i], tempq);
	}
      temp = fmod(tempqinv*(double)(basis[i]-1), tempq);
      U[i] = fmod(tempqinv*Ac[i]+temp*U[i], tempq);

    }
  /* compute Ac in positive representation */
  mpz_set_d(mp_Ac, U[len-1]); 
  for (j = len-2; j >= 0; j--)
    {
      mpz_mul_ui(mp_Ac, mp_Ac, basis[j]);
      mpz_add_ui(mp_Ac, mp_Ac, (FiniteField)U[j]);
    }
  XFREE(U);

  return;
}



/*
 *
 * Calling Sequence:
 *   cmbasis <-- combBasis(basislen, basis)
 * 
 * Summary:
 *   Compute the special combination of a RNS basis
 * 
 * Description:
 *   Let 'basis' be RNS basis. The function computes an array cmbasis 
 *   satisfying 
 *   cmbasis[0] = 0, cmbasis[i] = mod(1/(basis[0]*...*basis[i-1]), basis[i])
 *                   (i = 1..basislen-1)
 * 
 * Input:
 *   basislen: long, dimension of RNS basis
 *      basis: 1-dim FiniteField array length basislen, RNS basis
 * 
 * Return:
 *   cmbasis: 1-dim FiniteField array length basislen, shown as above
 *
 */

FiniteField *
combBasis (const long basislen, const FiniteField *basis)
{
  long i, j;
  double dtemp;
  mpz_t mp_prod, mp_q;
  FiniteField *cmbasis;

  cmbasis = XMALLOC(FiniteField, basislen);
  cmbasis[0] = 0;
  mpz_init(mp_prod);
  mpz_init(mp_q);
  for (i = 1; i < basislen; i++)
    {
      dtemp = fmod((double)basis[0], (double)basis[i]);
      for (j = 1; j <= i-1; j++)
	dtemp = fmod(dtemp*(double)basis[j], (double)basis[i]);
      mpz_set_ui(mp_q, basis[i]);
      mpz_set_d(mp_prod, dtemp);
      mpz_invert(mp_prod, mp_prod, mp_q);
      cmbasis[i] = mpz_get_ui(mp_prod);
    }
  mpz_clear(mp_prod);
  mpz_clear(mp_q);

  return cmbasis;
}



/*
 *
 * Calling Sequence:
 *   cumprod <-- cumProd(basislen, basis, extbasislen, extbasis)
 *
 * Summary:
 *   Compute the representation of the combination of elements of one RNS basis
 *   in another RNS basis
 * 
 * Description:
 *   Let 'basis' be one RNS basis with dimension basislen, and 'extbasis' be 
 *   another RNS basis with dimension extbasislen. The function computes an 
 *   array cumprod length extbasislen satisfying
 *   cumprod[i] = modp(-basis[0]*...*basis[basislen-1], extbasis[i]),
 *   i = 0..extbasislen-1
 *
 * Input: 
 *      basislen: long, dimension of RNS basis 'basis'
 *         basis: 1-dim FiniteField array length basislen, one RNS basis
 *   extbasislen: long, dimension of RNS basis 'extbasis'
 *      extbasis: 1-dim FiniteField array length basislen, another RNS basis
 * 
 * Return: 
 *   cumprod: 1-dim double array length extbasislen, shown above
 *
 */

double *
cumProd (const long basislen, const FiniteField *basis, \
	 const long extbasislen, const FiniteField *extbasis)
{
  long i, j;
  double dtemp, dextbasis;
  double *cumprod;

  cumprod = XMALLOC(double, extbasislen);
  for (i = 0; i < extbasislen; i++)
    {
      dextbasis = (double)extbasis[i];
      cumprod[i] = fmod((double)basis[0], dextbasis);
      for (j = 1; j < basislen; j++)
	{
	  dtemp = fmod((double)basis[j], dextbasis);
	  cumprod[i] = fmod(cumprod[i]*dtemp, dextbasis);
	}
      cumprod[i] = dextbasis-cumprod[i];
    }

  return cumprod;
}


/*
 *
 * Calling Sequence:
 *   basiscmb <-- findRNS(RNS_bound, mp_maxInter, len)
 *
 * Summary:
 *   Find a RNS basis and its special combination
 *   
 * Description:
 *   Given RNS_bound, the upper bound of the RNS basis, and mp_maxInter, the
 *   function finds a best RNS basis and a combination of that basis.
 *   
 *   The RNS basis 'basis' has the property:
 *   - its elements are all primes
 *   - basis[0] is the largest prime among all the primes at most RNS_bound
 *   - basis[i+1] is the next prime smaller than basis[i] (i = 0..len-2)
 *   - basis[0]*basis[1]*...*basis[len-1] >= mp_maxInter
 *
 *   After finding 'basis', the functions also computes the combination of
 *   'basis' as the operations in function combBasis.
 *
 * Input:  
 *     RNS_bound: FiniteField, the upper bound of the RNS basis
 *   mp_maxInter: mpz_t, the lower bound for the product of elements of basis
 *
 * Return:
 *   basiscmb: 2-dim FiniteField array, dimension 2 x len, where
 *           - basiscmb[0] represents the RNS basis
 *           - basiscmb[1] represents the special combination of basis 
 *
 * Output:
 *   len: pointer to a long int, storing the dimension of the computed
 *        RNS basis
 *
 */

FiniteField **
findRNS (const FiniteField RNS_bound, const mpz_t mp_maxInter, long *length)
{
  long i, j, len=0;
  double prod;
  mpz_t mp_l, mp_prod, mp_q;
  FiniteField **qqinv;

  mpz_init_set_ui(mp_prod, 1);
  mpz_init_set_ui(mp_l, RNS_bound);
  qqinv = XMALLOC(FiniteField *, 2);
  qqinv[0] = NULL;
  while (mpz_cmp(mp_maxInter, mp_prod) > 0)
    {
      ++len;
      qqinv[0] = XREALLOC(FiniteField, qqinv[0], len);
      while (mpz_probab_prime_p(mp_l, 10) == 0) { mpz_sub_ui(mp_l, mp_l, 1); }
      qqinv[0][len-1] = mpz_get_ui(mp_l);
      mpz_sub_ui(mp_l, mp_l, 1);
      mpz_mul_ui(mp_prod, mp_prod, qqinv[0][len-1]);
    }
  mpz_clear(mp_prod);
  mpz_clear(mp_l);
  qqinv[1] = XMALLOC(FiniteField, len);
  qqinv[1][0] = 0;
  mpz_init(mp_prod);
  mpz_init(mp_q);
  for (i = 1; i < len; i++)
    {  
      prod = (double)(qqinv[0][0] % qqinv[0][i]);
      for (j = 1; j <= i-1; j++) 
	prod = fmod(prod*(double)qqinv[0][j], (double)qqinv[0][i]);
      mpz_set_ui(mp_q, qqinv[0][i]);
      mpz_set_d(mp_prod, prod);
      mpz_invert(mp_prod, mp_prod, mp_q);
      qqinv[1][i] = mpz_get_ui(mp_prod);
    }
  mpz_clear(mp_prod);
  mpz_clear(mp_q);  
  *length = len;

  return qqinv;
}



/*
 *
 * Calling Sequence:
 *   maxInter(mp_prod, mp_alpha, n, mp_b)
 *
 * Summary:
 *   Compute the maximum interval of positive and negative results of a
 *   matrix-matrix or matrix-vector product
 * 
 * Description:
 *   Let mp_alpha be the maximum magnitude of a m x n matrix A, mp_prod-1 be
 *   the maximum magnitude of a n x k matrix C. The function computes the 
 *   maximum interval of positive and negative entries of A.C. That is, the 
 *   function computes mp_b satisfying
 *   (mp_b-1)/2 = n*mp_alpha*(mp_prod-1)
 * 
 * Input:
 *    mp_prod: mpz_t, mp_prod-1 be the maximum magnitude of matrix C
 *   mp_alpha: mpz_t, maximum magnitude of matrix A
 *          n: long, column dimension of A
 *
 * Output:
 *   mp_b: mpz_t, shown above
 *
 */

void
maxInter (const mpz_t mp_prod, const mpz_t mp_alpha, const long n, mpz_t mp_b)
{
  mpz_t mp_temp;

  mpz_init(mp_temp);
  mpz_sub_ui(mp_temp, mp_prod, 1);
  mpz_set(mp_b, mp_alpha);
  mpz_mul_ui(mp_b, mp_b, n);
  mpz_mul(mp_b, mp_b, mp_temp);
  mpz_mul_ui(mp_b, mp_b, 2);
  mpz_add_ui(mp_b, mp_b, 1);
  mpz_clear(mp_temp);
}



/*
 *
 * Calling Sequence:
 *   maxExtInter(mp_alpha, n, mp_b)
 *
 * Summary:
 *   Compute the maximum interval of positive and negative results for 
 *   lifting
 * 
 * Description:
 *   Let mp_alpha be the maximum magnitude of a m x n matrix A, 
 *   The function computes the mp_b satisfying
 *   (mp_b-1)/2 = n*mp_alpha+1
 * 
 * Input:
 *   mp_alpha: mpz_t, maximum magnitude of matrix A
 *          n: long, column dimension of A
 *
 * Output:
 *   mp_b: mpz_t, shown above
 *
 */

void
maxExtInter (const mpz_t mp_alpha, const long n, mpz_t mp_b)
{
  mpz_set_ui(mp_b, 1);
  mpz_addmul_ui(mp_b, mp_alpha, n);
  mpz_mul_ui(mp_b, mp_b, 2);
  mpz_add_ui(mp_b, mp_b, 1);
}



/*
 *
 * Calling Sequence:
 *   bdcoeff <-- repBound(len, basis, cmbasis)
 *
 * Summary:
 *   Compute the mix radix coefficients of a special integer in a RNS basis
 *
 * Description:
 *   Given a RNS basis, suppose the product of elements in the basis be q, 
 *   then this RNS basis is able to represent integers lying in 
 *   [-(q-1)/2, (q-1)/2] and [0, q-1] respectively with symmetric 
 *   representation and positive representation. To transfer the result from
 *   positive representation to symmetric representation, the function 
 *   computes the mix radix coefficients of the boundary value (q-1)/2 in the
 *   positive representation.
 *
 *   Let RNS basis be P. The function computes coefficient array U, such that
 * (q-1)/2 = U[0] + U[1]*P[0] + U[2]*P[0]*P[1] +...+ U[len-1]*P[0]*...*P[len-2]
 *
 * Input: 
 *       len: long, dimension of RNS basis
 *     basis: 1-dim FiniteField array length len, RNS basis
 *   cmbasis: 1-dim FiniteField array length len, computed by function 
 *            combBasis, inverses of special combination of RNS basis
 *
 * Output:
 *   bdcoeff: 1-dim FiniteField array length len, the coefficient array U above
 *
 */

FiniteField *
repBound (const long len, const FiniteField *basis, const FiniteField *cmbasis)
{
  long i, j;
  double dtemp;
  mpz_t mp_bd, mp_prod;
  FiniteField *bdcoeff;
  const FiniteField *q, *qinv;

  q = basis;
  qinv = cmbasis;

  /* set the bound of transformation from positive to negative */
  mpz_init(mp_prod);
  basisProd(len, q, mp_prod);
  mpz_init(mp_bd);
  mpz_sub_ui(mp_bd, mp_prod, 1);
  mpz_divexact_ui(mp_bd, mp_bd, 2);

  /* compute the coeffcients of bound of mix radix and store in bdcoeff */
  bdcoeff = XMALLOC(FiniteField, len);
  bdcoeff[0] = mpz_fdiv_ui(mp_bd, q[0]);
  for (i = 1; i < len; i++)
    {
      dtemp = (double)bdcoeff[i-1];
      for (j = i-2; j >= 0; j--)
	{ 
	  dtemp = fmod(dtemp*q[j], (double)q[i]);
	  dtemp = fmod(dtemp+bdcoeff[j], (double)q[i]);
	}
      bdcoeff[i] = mpz_fdiv_ui(mp_bd, q[i]);
      dtemp = fmod((double)bdcoeff[i]-dtemp, (double)q[i]);
      if (dtemp < 0) { dtemp = q[i]+dtemp; }
      bdcoeff[i] = (FiniteField)fmod(dtemp*qinv[i], (double)q[i]);
    }
  mpz_clear(mp_bd);
  mpz_clear(mp_prod);

  return bdcoeff;
}



/*
 * Calling Sequence:
 *   bd <-- RNSbound(n)
 *
 * Summary:
 *   Compute the upper bound of a RNS basis
 *
 * Description:
 *   Given a m x n mod p matrix A, and a n x k mod p matrix B, the maximum 
 *   magnitude of A.B is n*(p-1)^2. To use BLAS, it is needed that 
 *   n*(p-1)^2 <= 2^53-1 to make the result of product correct.
 *
 *   The function computes an integer bd, such that 
 *      n*(bd-1)^2 <= 2^53-1 and n*((bd+1)-1)^2 > 2^53-1
 * 
 * Input:
 *   n: long, column dimension of matrix A
 *
 * Output:
 *   bd: FiniteField, shown above
 *
 */

FiniteField
RNSbound (const long n)
{
  FiniteField bd;
  mpz_t mp_n, mp_d, mp_q;

  mpz_init(mp_n);
  mpz_init_set_ui(mp_d, n);
  mpz_init(mp_q);
  mpz_ui_pow_ui(mp_n, 2, 53);
  mpz_sub_ui(mp_n, mp_n, 1);
  mpz_fdiv_q(mp_q, mp_n, mp_d);
  mpz_sqrt(mp_q, mp_q);
  bd = mpz_get_ui(mp_q)+1;
  mpz_clear(mp_n);
  mpz_clear(mp_d);
  mpz_clear(mp_q);

  return bd;
}







