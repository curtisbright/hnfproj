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


/* this example shows how to call the nonsingular solving function to solve
 * a randomly generated nonsingular system
 */


#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "iml.h"

long * randomLongMat (const long n, const long m, const long bd);
mpz_t * randomMPMat (const long n, const long m, const mpz_t mp_bd);

int main(void)
{
  long i, j, n, bd;
  long *A;
  mpz_t mp_bd, mp_D;
  mpz_t *mp_B, *mp_N;

  /* generate a n x n random left hand side matrix A */
  n = 5;   /* set the input matrix dimension */
  bd = 7;  /* any entry e in A satisfying -2^bd < e < 2^bd */
  A = randomLongMat(n, n, bd);
  
  /* generate a n x 1 random right hand side matrix mp_B */
  mpz_init(mp_bd);
  mpz_ui_pow_ui(mp_bd, 2, bd);     /* any entry e in mp_B satisfying */
  mp_B = randomMPMat(n, 1, mp_bd); /* -2^bd < e < 2^bd               */

  /* print the input system */
  fprintf(stdout, "Input Matrices:\n");
  fprintf(stdout, "A:\n");
  for (i = 0; i < n; i++)
    {
      fprintf(stdout, "  ");
      for (j = 0; j < n; j++)
	fprintf(stdout, "%ld\t", A[i*n+j]);
      fprintf(stdout, "\n");
    }
  fprintf(stdout, "\n B:\n");
  for (i = 0; i < n; i++)
    gmp_fprintf(stdout, "  %Zd\n", mp_B[i]);

  /* allocate space for numerator vector mp_N and denominator mp_D */
  mpz_init(mp_D);
  mp_N = (mpz_t *) malloc(n*sizeof(mpz_t));
  for (i = 0; i < n; i++) { mpz_init(mp_N[i]); }

  /* solve system AX = mp_B */
  nonsingSolvMM(RightSolu, n, 1, A, mp_B, mp_N, mp_D);

  /* print the computation result */
  fprintf(stdout, "\nSolution:\n");
  fprintf(stdout, "  Denominator: \n");
  gmp_fprintf(stdout, "  %Zd\n", mp_D);
  fprintf(stdout, "\n  Numerator: \n");
  for (i = 0; i < n; i++)
    gmp_fprintf(stdout, "  %Zd\n", mp_N[i]); 

  /* free the memory */
  for (i = 0; i < n; i++) { mpz_clear(mp_B[i]); mpz_clear(mp_N[i]); }
  { free(mp_N); free(mp_B); free(A); }
  mpz_clear(mp_D);

  return 0;
}



/* generate a n x m random dense signed long matrix with entries lying in
 * (-2^bd, 2^bd) 
 */
long *
randomLongMat (const long n, const long m, const long bd)
{
  long i, j;
  mpz_t mp_rand, mp_sign;
  gmp_randstate_t state;
  unsigned long seed;
  FILE* devrandom;
  long* M;
  static unsigned long inc = 0;

  M = (long *) malloc(n*m*sizeof(long));
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

      for (i = 0; i < n*m; i++)
	{
	  mpz_urandomb(mp_rand, state, bd);
	  mpz_urandomb(mp_sign, state, 1);
	  if (mpz_sgn(mp_sign) == 0) { mpz_neg(mp_rand, mp_rand); }
	  M[i] = mpz_get_si(mp_rand); 
	}
      fclose(devrandom);
    }
  mpz_clear(mp_rand);
  gmp_randclear(state);
  mpz_clear(mp_sign);
  return M;
}


/* generate a n x m random dense mpz_t matrix with entries lying in 
 * (-mp_bd, mp_bd)
 */
mpz_t *
randomMPMat (const long n, const long m, const mpz_t mp_bd)
{
  long i, j;
  mpz_t mp_rand, mp_sign;
  gmp_randstate_t state;
  unsigned long seed;
  mpz_t *mp_B;
  FILE* devrandom;

  mp_B = (mpz_t *) malloc(n*m*sizeof(mpz_t));
  mpz_init(mp_sign);
  mpz_init(mp_rand);
  gmp_randinit_default(state);

  /* generate random seed using /dev/random */
  if ((devrandom = fopen("/dev/urandom", "r")) != NULL)
    {
      for (i = 0; i < n*m; i++)
	{
	  fread(&seed, sizeof(seed), 1, devrandom);
	  gmp_randseed_ui(state, seed);
	  mpz_urandomm(mp_rand, state, mp_bd);
	  mpz_urandomb(mp_sign, state, 1);
	  if (mpz_sgn(mp_sign) == 0) { mpz_neg(mp_rand, mp_rand); }
	  mpz_init_set(mp_B[i], mp_rand); 
	}
      fclose(devrandom);
    }
  mpz_clear(mp_rand);
  mpz_clear(mp_sign);
  gmp_randclear(state);
  return mp_B;
}
