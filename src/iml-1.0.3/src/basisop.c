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



#include "basisop.h"

/*
 * Calling Sequence:
 *   Dmod(p, A, n, m, lda)
 * 
 * Summary:
 *   Perform mod p operation inplace over a double matrix/submatrix
 *
 * Description:
 *   Given a n x m matrix or submatrix A and an integer p, this function uses
 *   fmod function to do operation A mod p inplace. Although the type of input
 *   modulus p in this function is Double, it is a cast of an integer.
 *
 * Input:
 *     p: Double, module
 *     A: 1-dim Double array length n*m, representation of a n x m double 
 *        matrix/submatrix
 *     n: long, row dimension of A
 *     m: long, column dimension of A
 *   lda: long, number of entries between two continuous rows of A, help track
 *        the pointer if A is a submatrix
 *
 */

void
Dmod (const Double p, Double *A, const long n, const long m, const long lda)
{
  long i, j;
  Double temp;

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      *(A+i*lda+j) = (temp = fmod(*(A+i*lda+j), p)) >= 0 ? temp : p+temp;
}

  
/*
 * Calling Sequence:
 *   DCopy(n, m, A, lda, B, ldb)
 * 
 * Summary:
 *   Copy Double matrix/submatrix A to another matrix/submatrix B
 *
 * Input:
 *     n: long, row dimension of matrix/submatrix A and B
 *     m: long, column dimension of matrix/submatrix A and B
 *     A: 1-dim Double array length n*m, representation of a n x m Double 
 *        matrix/submatrix
 *   lda: long, stride of two continuous rows of A
 *     B: 1-dim Double array length n*m, representation of a n x m Double 
 *        matrix/submatrix
 *   ldb: long, stride of two continuous rows of B
 *
 */

void 
DCopy (const long n, const long m, const Double* A, const long lda, \
       Double* B, const long ldb)
{
  long i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      *(B+i*ldb+j) = *(A+i*lda+j);
}



/*
 * Calling Sequence:
 *   randomDb(n, m, bd)
 * 
 * Summary:
 *   Generate a random Double dense matrix
 * 
 * Description:
 *   Given the bound bd, this function generates a random Double dense matrix
 *   such that any entry e satisfies 0 <= e <= 2^bd-1
 *
 * Input:
 *    n: long, row dimension of the generated matrix
 *    m: long, column dimension of the generated matrix
 *   bd: long, magnitude bound of entries
 *
 * Return:
 *   M: 1-dim Double array length n*m, representation of a n x m random matrix
 *
 */

Double *
randomDb (const long n, const long m, const long bd)
{
  long i;
  mpz_t mp_rand;
  gmp_randstate_t state;
  unsigned long seed;
  FILE* devrandom;
  Double* M;
  static unsigned long inc = 0;
  size_t rc;

#if HAVE_TIME_H
  time_t t1;
#endif

  M = XMALLOC(Double, n*m);
  mpz_init(mp_rand);
  gmp_randinit_default(state);
  mpz_set_ui(mp_rand, 5);
  
  seed = 3828173;

  /* generate random seed using /dev/random */
  if ((devrandom = fopen("/dev/urandom", "r")) != NULL)
    {
      rc = fread(&seed, sizeof(seed), 1, devrandom);
      (void)rc;
      fclose(devrandom);
    }
  else
    {
#if HAVE_TIME_H
       time(&t1);
       seed = (unsigned long)t1;
#endif
    }

  seed += inc;
  inc += 1;
  gmp_randseed_ui(state, seed);
  for (i = 0; i < n*m; i++)
    {
       mpz_urandomb(mp_rand, state, bd);
       M[i] = mpz_get_d(mp_rand); 
    }
  mpz_clear(mp_rand);
  gmp_randclear(state);
  return M;
}



/*
 * Calling Sequence:
 *   RandomPrime(lb, hb)
 *
 * Summary:
 *   Generate a random prime p satisfying 2^lb <= p <= 2^hb-1
 *
 * Input:
 *   lb: FiniteField, lower bound of the prime
 *   hb: FiniteField, upper bound of the prime
 *
 * Return:
 *   p: FiniteField, randomly generated prime
 *
 * Note: 
 *   lb and hb can not exceed the bit-length of unsigned long.
 *
 */

FiniteField
RandPrime (const FiniteField lb, const FiniteField hb)
{
  mpz_t mp_rand, mp_temp, mp_lb, mp_hb;
  gmp_randstate_t state;
  FiniteField p;
  unsigned long seed;
  FILE* devrandom;
  static unsigned long inc = 0;
  size_t rc;

#if HAVE_TIME_H
  time_t t1;
#endif

  { mpz_init(mp_rand); mpz_init(mp_temp); } 
  { mpz_init(mp_lb);  mpz_init(mp_hb); }
  mpz_ui_pow_ui(mp_lb, 2, lb);
  mpz_ui_pow_ui(mp_hb, 2, hb);
  mpz_sub(mp_temp, mp_hb, mp_lb);
  gmp_randinit_default(state);

  seed = 3828173;

  /* generate random seed using /dev/urandom */
  if ((devrandom = fopen("/dev/urandom", "r")) != NULL)
    {
      rc = fread(&seed, sizeof(seed), 1, devrandom);
      (void)rc;
      fclose(devrandom);
    }
  else
    {
#if HAVE_TIME_H
      time(&t1);
      seed = (unsigned long)t1;
#endif
    }

  seed += inc;
  inc += 1;
  gmp_randseed_ui(state, seed);
  mpz_urandomm(mp_rand, state, mp_temp); /* 0 <= mp_rand <= 2^hb-2^lb-1 */
  mpz_add(mp_rand, mp_rand, mp_lb);      /* 2^lb <= mp_rand <= 2^hb-1 */
  while (mpz_probab_prime_p(mp_rand, 10) == 0) 
    mpz_sub_ui(mp_rand, mp_rand, 1); 
  p = mpz_get_ui(mp_rand);

  { mpz_clear(mp_rand); mpz_clear(mp_temp); }
  { mpz_clear(mp_lb); mpz_clear(mp_hb); }
  gmp_randclear(state);
  return p;
}



/*
 * Calling Sequence:
 *   max <-- maxMagnLong(A, n, m, lda)
 *
 * Summary:
 *   Compute maximum magnitude of a n x m signed long matrix/submatrix A
 *
 * Input:
 *     A: 1-dim long array length n*m, representation of a n x m 
 *        matrix/submatrix
 *     n: long, row dimension of A
 *     m: long, column dimension of A
 *   lda: long, number of entries of two continuous rows of A (normally m)
 *
 * Return: 
 *   max: long, maximum magnitude of A
 *
 */

long 
maxMagnLong (const long *A, const long n, const long m, const long lda)
{
  long i, j, temp, max=0;

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      if ((temp = labs(A[i*lda+j])) > max) { max = temp; }
  return max;
}



/*
 *
 * Calling Sequence:
 *   maxMagnMP(mp_A, n, m, lda, mp_max)
 *
 * Summary:
 *   Compute maximum magnitude of a n x m mpz_t matrix/submatrix mp_A
 *
 * Input:
 *   mp_A: 1-dim mpz_t array length n*m, representation of a n x m 
 *         matrix/submatrix
 *      n: long, row dimension of mp_A
 *      m: long, column dimension of mp_A
 *    lda: long, number of entries of two continuous rows of mp_A (normally m)
 *
 * Output: 
 *   mp_max: mpz_t, maximum magnitude of mp_A
 *
 */

void
maxMagnMP (mpz_t *mp_A, const long n, const long m, const long lda, \
	   mpz_t mp_max)
{
  long i, j;

  mpz_set_ui(mp_max, 0);
  for (i = 0; i < n; i++) 
    for (j = 0; j < m; j++)
      if (mpz_cmpabs(mp_A[i*lda+j], mp_max) > 0) 
	mpz_abs(mp_max, mp_A[i*lda+j]);
  return;
}


/*
 * Calling Sequence:
 *   scalCpMP(n, m, lda, ldm, mp_a, mp_A, mp_M)
 *
 * Summary:
 *   Copy a submatrix with scale from one matrix to the other one
 *
 * Description:
 *   Copy with scalar n x m submatrix from A to M. The submatrix starts
 *   from position pointed by mp_A with stride lda. Then the submatrix is 
 *   copied to position pointed by mp_M with stride ldm, meanwhile, the 
 *   entries of the submatrix are scaled by factor mp_a.
 * 
 * Input: 
 *      n: long, row dimension of the submatrix
 *      m: long, column dimension of the submatrix
 *    lda: long, stride of two continuous rows of the submatrix in A
 *    ldm: long, stride of two continuous rows of the submatrix in M
 *   mp_a: mpz_t, scalar factor
 *   mp_A: mpz_t pointer, start point of the submatrix in A
 *   mp_M: mpz_t pointer, start point of the submatrix in M
 *
 */

void 
scalCpMP (const long n, const long m, const long lda, const long ldm, \
	  const mpz_t mp_a, mpz_t *mp_A, mpz_t *mp_M)
{
  long i, j;

  if (mpz_cmp_ui(mp_a, 1) != 0)
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
	{
	  mpz_set(mp_M[i*ldm+j], mp_A[i*lda+j]);
	  mpz_mul(mp_M[i*ldm+j], mp_M[i*ldm+j], mp_a);
	}
  else
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
	mpz_set(mp_M[i*ldm+j], mp_A[i*lda+j]);
}
      
