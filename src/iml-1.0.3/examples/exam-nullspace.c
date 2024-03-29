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

#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "iml.h"

long *randomLongMat (const long n, const long m, const long bd);
mpz_t *randomMPMat (const long n, const long m, const mpz_t mp_bd);

int
main (void)
{
  long i, j, n, m, bd, s, *A;
  mpz_t *mp_B, *mp_N;

  /* generate a n x m random left hand side matrix A */
  n = 5;
  m = 10;
  bd = 34;			/* entris of A satisfying -2^bd < e < 2^bd */
  A = randomLongMat (n, m, bd);
  fprintf (stdout, "Input system:\n");
  for (i = 0; i < n; i++)
    {
      fprintf (stdout, "  ");
      for (j = 0; j < m; j++)
	fprintf (stdout, "%ld\t", A[i * m + j]);
      fprintf (stdout, "\n");
    }
  s = nullspaceLong (n, m, A, &mp_N);
  fprintf (stdout, "Dimension of nullspace: ");
  fprintf (stdout, " %ld\n", s);
  for (i = 0; i < m; i++)
    {
      for (j = 0; j < s; j++)
	gmp_fprintf (stdout, "  %Zd", mp_N[i * s + j]);
      fprintf (stdout, "\n");
    }
  free (A);
  for (i = 0; i < m * s; i++)
    mpz_clear (mp_N[i]);
  free (mp_N);
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
  FILE *devrandom;
  long *M;
  static unsigned long inc = 0;

  M = (long *) malloc (n * m * sizeof (long));
  mpz_init (mp_sign);
  mpz_init (mp_rand);
  gmp_randinit_default (state);
  seed = 387439;


  /* generate random seed using /dev/random */
  if ((devrandom = fopen ("/dev/urandom", "r")) != NULL)
    {
      fread (&seed, sizeof (seed), 1, devrandom);
      fclose (devrandom);
    }
  seed += inc;
  inc += 1;
  gmp_randseed_ui (state, seed);


  for (i = 0; i < n * m; i++)
    {
      mpz_urandomb (mp_rand, state, bd);
      mpz_urandomb (mp_sign, state, 1);
      if (mpz_sgn (mp_sign) == 0)
	mpz_neg (mp_rand, mp_rand);
      M[i] = mpz_get_si (mp_rand);
    }
  mpz_clear (mp_rand);
  gmp_randclear (state);
  mpz_clear (mp_sign);
  return M;
}
