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


/* test file for small entry (signed long) input matrix */

#include "gmp.h"

#include "../src/common.h"
#include "../src/iml.h"

long tstcertsolveLong(const long certflag, const long nullcol, const long n, \
		      const long m, const double den, const long alpha, \
		      const mpz_t mp_beta);
long *
randomLongMat (const long n, const long m, const double den, const long bd);
mpz_t *
randomMPMat (const long n, const long m, const double den, const mpz_t mp_bd);

int main(void)
{
  double den;
  long certflag, n, m, count, alpha, ret, nullcol;
  mpz_t mp_beta;
  mpz_init(mp_beta);
  alpha = 6; nullcol = 11;
  mpz_ui_pow_ui(mp_beta, 2, 6);

  /* test square dense matrix */
  n = 20; m = 20; den = 1; certflag = 1;
  ret = tstcertsolveLong(certflag, nullcol, n, m, den, alpha, mp_beta);
  if (ret == 1) { return 1; }

  /* test rectangular dense matrix */
  n = 10; m = 30;
  den = 1; 
  certflag = 1; 

  ret = tstcertsolveLong(certflag, nullcol, n, m, den, alpha, mp_beta);
  if (ret == 1) { return 1; }

  certflag = 0; 

  ret = tstcertsolveLong(certflag, nullcol, n, m, den, alpha, mp_beta);
  if (ret == 1) { return 1; }

  /* test rectangular sparse matrix */
  den = 0.1; 
  certflag = 1; 

  ret = tstcertsolveLong(certflag, nullcol, n, m, den, alpha, mp_beta);
  if (ret == 1) { return 1; }

  certflag = 0; 

  ret = tstcertsolveLong(certflag, nullcol, n, m, den, alpha, mp_beta);
  if (ret == 1) { return 1; }

  /* pass test */
  mpz_clear(mp_beta);
  return 0;
}


long tstcertsolveLong(const long certflag, const long nullcol, const long n, \
		      const long m, const double den, const long alpha, \
		      const mpz_t mp_beta)
{
  long i, j, l, r, ret, n1, m1;
  mpz_t mp_D, mp_DZ, mp_temp, mp_temp1;
  long *A;
  mpz_t *mp_b, *mp_N, *mp_NZ=NULL;

  A = randomLongMat(n, m, den, alpha);
  mp_b = randomMPMat(n, 1, den, mp_beta);
  mp_N = XMALLOC(mpz_t, m);
  for (i = 0; i < m; i++) { mpz_init(mp_N[i]); }
  mpz_init(mp_D); mpz_init(mp_DZ);
  if (certflag == 1) 
    { 
      mp_NZ = XMALLOC(mpz_t, n); 
      for (i = 0; i < n; i++) { mpz_init(mp_NZ[i]); }
    }

  ret = certSolveLong(certflag, n, m, A, mp_b, mp_N, mp_D, mp_NZ, mp_DZ);

  // test solution
  mpz_init(mp_temp); mpz_init(mp_temp1);
  if (ret != 3)
    {
      // solution exists, check solution
      for (l = 0; l < n; l++)
	{
	  mpz_mul_si(mp_temp, mp_N[0], A[l*m]);
	  for (j = 1; j < m; j++)
	    {
	      mpz_set_si(mp_temp1, A[l*m+j]);
	      mpz_addmul(mp_temp, mp_N[j], mp_temp1);
	    }
	  mpz_mul(mp_temp1, mp_D, mp_b[l]);
	  if (mpz_cmp(mp_temp, mp_temp1) != 0) { return 1; }
	}
      if ((certflag == 1) && (ret == 1))
	{
	  // check certificate z.A
	  for (l = 0; l < m; l++)
	    {
	      mpz_mul_si(mp_temp, mp_NZ[0], A[l]);
	      for (j = 1; j < n; j++)
		{ 
		    mpz_set_si(mp_temp1, A[j*m+l]);
		    mpz_addmul(mp_temp, mp_NZ[j], mp_temp1);
		}
	      if (!mpz_divisible_p(mp_temp, mp_DZ)) { return 1; }
	    }
	  // check denom(z.b)
	  mpz_mul(mp_temp, mp_NZ[0], mp_b[0]);
	  for (j = 1; j < n; j++)
	    mpz_addmul(mp_temp, mp_NZ[j], mp_b[j]);
	  mpz_gcd(mp_temp1, mp_temp, mp_DZ);
	  mpz_divexact(mp_temp1, mp_DZ, mp_temp1);
	  if (mpz_cmp(mp_temp1, mp_D) != 0) { return 1; }
	}
    }
  else if (certflag == 1) 
    {
      // solution doesn't exist
      for (l = 0; l < m; l++)
	{
	  mpz_mul_si(mp_temp, mp_NZ[0], A[l]);
	    for (j = 1; j < n; j++)
	      { 
		mpz_set_si(mp_temp1, A[j*m+l]);
		mpz_addmul(mp_temp, mp_NZ[j], mp_temp1);
	      }
	    if (mpz_cmp_ui(mp_temp, 0) != 0) { return 1; }
	  }
      mpz_mul(mp_temp, mp_NZ[0], mp_b[0]);
      for (j = 1; j < n; j++)
	mpz_addmul(mp_temp, mp_NZ[j], mp_b[j]);
      if (!mpz_cmp_ui(mp_temp, 0)) { return 1; }
    }
  mpz_clear(mp_temp); mpz_clear(mp_temp1);


  ret = certSolveRedLong(certflag, nullcol, n, m, A, mp_b, mp_N, mp_D, \
			 mp_NZ, mp_DZ);

  // test solution
  mpz_init(mp_temp); mpz_init(mp_temp1);
  if (ret != 3)
    {
      // solution exists, check solution
      for (l = 0; l < n; l++)
	{
	  mpz_mul_si(mp_temp, mp_N[0], A[l*m]);
	  for (j = 1; j < m; j++)
	    {
	      mpz_set_si(mp_temp1, A[l*m+j]);
	      mpz_addmul(mp_temp, mp_N[j], mp_temp1);
	    }
	  mpz_mul(mp_temp1, mp_D, mp_b[l]);
	  if (mpz_cmp(mp_temp, mp_temp1) != 0) { return 1; }
	}
      if ((certflag == 1) && (ret == 1))
	{
	  // check certificate z.A
	  for (l = 0; l < m; l++)
	    {
	      mpz_mul_si(mp_temp, mp_NZ[0], A[l]);
	      for (j = 1; j < n; j++)
		{ 
		    mpz_set_si(mp_temp1, A[j*m+l]);
		    mpz_addmul(mp_temp, mp_NZ[j], mp_temp1);
		}
	      if (!mpz_divisible_p(mp_temp, mp_DZ)) { return 1; }
	    }
	  // check denom(z.b)
	  mpz_mul(mp_temp, mp_NZ[0], mp_b[0]);
	  for (j = 1; j < n; j++)
	    mpz_addmul(mp_temp, mp_NZ[j], mp_b[j]);
	  mpz_gcd(mp_temp1, mp_temp, mp_DZ);
	  mpz_divexact(mp_temp1, mp_DZ, mp_temp1);
	  if (mpz_cmp(mp_temp1, mp_D) != 0) { return 1; }
	}
    }
  else if (certflag == 1) 
    {
      // solution doesn't exist
      for (l = 0; l < m; l++)
	{
	  mpz_mul_si(mp_temp, mp_NZ[0], A[l]);
	    for (j = 1; j < n; j++)
	      { 
		mpz_set_si(mp_temp1, A[j*m+l]);
		mpz_addmul(mp_temp, mp_NZ[j], mp_temp1);
	      }
	    if (mpz_cmp_ui(mp_temp, 0) != 0) { return 1; }
	  }
      mpz_mul(mp_temp, mp_NZ[0], mp_b[0]);
      for (j = 1; j < n; j++)
	mpz_addmul(mp_temp, mp_NZ[j], mp_b[j]);
      if (!mpz_cmp_ui(mp_temp, 0)) { return 1; }
    }
  mpz_clear(mp_temp); mpz_clear(mp_temp1);


  { mpz_clear(mp_D); mpz_clear(mp_DZ); }
  for (i = 0; i < m; i++) { mpz_clear(mp_N[i]); } { free(mp_N); }
  if (certflag == 1)
    for (i = 0; i < n; i++) { mpz_clear(mp_NZ[i]); } { free(mp_NZ); }
  free(A);
  for (i = 0; i < n; i++) { mpz_clear(mp_b[i]); } { free(mp_b); }
  return 0;
}



long *
randomLongMat (const long n, const long m, const double den, const long bd)
{
  long i, j, h,  flag;
  mpz_t mp_rand, mp_sign;
  gmp_randstate_t state;
  unsigned long seed;
  FILE* devrandom;
  long* M;
  static unsigned long inc = 0;

#if HAVE_TIME_H
  time_t t1;
#endif

#ifdef HAVE_FLOOR
  if (den <= 0.5) { h = (long)floor(1.0/den); flag = 1; }
  else { h = (long)floor(1.0/(1-den)); flag = 0; }
#else
  { h = 2; flag = 1; }
#endif

  M = XMALLOC(long, n*m);
  mpz_init(mp_sign);
  mpz_init(mp_rand);
  gmp_randinit_default(state);
  seed = 387439;


  /* generate random seed using /dev/random */
  if ((devrandom = fopen("/dev/urandom", "r")) != NULL)
    {
      fread(&seed, sizeof(seed), 1, devrandom);
      seed += inc;
      inc += 1;
      gmp_randseed_ui(state, seed);

      if (flag == 1)
	{
	  for (i = 0; i < n*m; i++)
	    {
	      mpz_urandomb(mp_rand, state, bd);
	      if ( (!mpz_divisible_ui_p(mp_rand, h)) || (den == 0.0)) 
		M[i] = 0;
	      else
		{
		  mpz_urandomb(mp_sign, state, 1);
		  if (mpz_sgn(mp_sign) == 0) { mpz_neg(mp_rand, mp_rand); }
		  M[i] = mpz_get_si(mp_rand); 
		}
	    }
	}
      else
	{
	  for (i = 0; i < n*m; i++)
	    {
	      mpz_urandomb(mp_rand, state, bd);
	      if ((mpz_divisible_ui_p(mp_rand, h) != 0) && (den != 1.0)) 
		M[i] = 0;
	      else
		{
		  mpz_urandomb(mp_sign, state, 1);
		  if (mpz_sgn(mp_sign) == 0) { mpz_neg(mp_rand, mp_rand); }
		  M[i] = mpz_get_si(mp_rand); 
		}
	    }
	}
      fclose(devrandom);
    }
  else
    {
#if HAVE_TIME_H
      time(&t1);
      seed = (unsigned long)t1;
#endif
      seed += inc;
      inc += 1;
      gmp_randseed_ui(state, seed);

      if (flag == 1)
	{
	  for (i = 0; i < n*m; i++)
	    {
	      mpz_urandomb(mp_rand, state, bd);
	      if ((!mpz_divisible_ui_p(mp_rand, h)) || (den == 0.0)) 
		M[i] = 0;
	      else
		{
		  mpz_urandomb(mp_sign, state, 1);
		  if (mpz_sgn(mp_sign) == 0) { mpz_neg(mp_rand, mp_rand); }
		  M[i] = mpz_get_si(mp_rand); 
		}
	    }
	}
      else
	{
	  for (i = 0; i < n*m; i++)
	    {
	      mpz_urandomb(mp_rand, state, bd);
	      if ((mpz_divisible_ui_p(mp_rand, h) != 0) && (den != 1.0))
		M[i] = 0;
	      else
		{
		  mpz_urandomb(mp_sign, state, 1);
		  if (mpz_sgn(mp_sign) == 0) { mpz_neg(mp_rand, mp_rand); }
		  M[i] = mpz_get_si(mp_rand); 
		}
	    }
	}
    }
  mpz_clear(mp_rand);
  gmp_randclear(state);
  mpz_clear(mp_sign);
  return M;
}




mpz_t *
randomMPMat (const long n, const long m, const double den, const mpz_t mp_bd)
{
  long i, j, h, flag;
  mpz_t mp_rand, mp_sign;
  gmp_randstate_t state;
  unsigned long seed;
  mpz_t *mp_B;
  FILE* devrandom;
  static unsigned long inc = 0;

#if HAVE_TIME_H
  time_t t1;
#endif

#ifdef HAVE_FLOOR
  if (den <= 0.5) { h = (long)floor(1.0/den); flag = 1; }
  else { h = (long)floor(1.0/(1-den)); flag = 0; }
#else
  { h = 2; flag = 1; }
#endif
  mp_B = XMALLOC(mpz_t, n*m);
  mpz_init(mp_sign);
  mpz_init(mp_rand);
  gmp_randinit_default(state);
  seed = 387439;

  /* generate random seed using /dev/random */
  if ((devrandom = fopen("/dev/urandom", "r")) != NULL)
    {
      fread(&seed, sizeof(seed), 1, devrandom);
      seed += inc;
      inc += 1;
      gmp_randseed_ui(state, seed);

      if (flag == 1)
	{
	  for (i = 0; i < n*m; i++)
	    {
	      mpz_urandomm(mp_rand, state, mp_bd);

	      /* with probability 1/h the entry is non-zero */
	      if ((!mpz_divisible_ui_p(mp_rand, h)) || (den == 0.0)) 
		mpz_init_set_ui(mp_B[i], 0);
	      else
		{
		  mpz_urandomb(mp_sign, state, 1);
		  if (mpz_sgn(mp_sign) == 0) { mpz_neg(mp_rand, mp_rand); }
		  mpz_init_set(mp_B[i], mp_rand); 
		}
	    }
	}
      else
	{
	  for (i = 0; i < n*m; i++)
	    {
	      mpz_urandomm(mp_rand, state, mp_bd);

	      /* with probability 1/h the entry is non-zero */
	      if ((mpz_divisible_ui_p(mp_rand, h) != 0) && (den != 1.0)) 
		mpz_init_set_ui(mp_B[i], 0);
	      else
		{
		  mpz_urandomb(mp_sign, state, 1);
		  if (mpz_sgn(mp_sign) == 0) { mpz_neg(mp_rand, mp_rand); }
		  mpz_init_set(mp_B[i], mp_rand); 
		}
	    }
	}
      fclose(devrandom);
    }
  else
    {
#if HAVE_TIME_H
      time(&t1);
      seed = (unsigned long)t1;
#endif
      seed += inc;
      inc += 1;
      gmp_randseed_ui(state, seed);

      if (flag == 1)
	{
	  for (i = 0; i < n*m; i++)
	    {
	      mpz_urandomm(mp_rand, state, mp_bd);

	      /* with probability 1/h the entry is non-zero */
	      if ((!mpz_divisible_ui_p(mp_rand, h)) || (den == 0.0)) 
		mpz_init_set_ui(mp_B[i], 0);
	      else
		{
		  mpz_urandomb(mp_sign, state, 1);
		  if (mpz_sgn(mp_sign) == 0) { mpz_neg(mp_rand, mp_rand); }
		  mpz_init_set(mp_B[i], mp_rand); 
		}
	    }
	}
      else
	{
	  for (i = 0; i < n*m; i++)
	    {
	      mpz_urandomm(mp_rand, state, mp_bd);

	      /* with probability 1/h the entry is non-zero */
	      if ((mpz_divisible_ui_p(mp_rand, h) != 0) && (den != 1.0)) 
		mpz_init_set_ui(mp_B[i], 0);
	      else
		{
		  mpz_urandomb(mp_sign, state, 1);
		  if (mpz_sgn(mp_sign) == 0) { mpz_neg(mp_rand, mp_rand); }
		  mpz_init_set(mp_B[i], mp_rand); 
		}
	    }
	}
    }
  mpz_clear(mp_rand);
  mpz_clear(mp_sign);
  gmp_randclear(state);
  return mp_B;
}
