/* ---------------------------------------------------------------------
 *
 * -- Integer Matrix Library (IML)
 *    (C) Copyright 2004, 2006 All Rights Reserved
 *
 * -- IML routines -- Version 1.0.1 -- November, 2006
 *
 * Author         : Zhuliang Chen
 * Contributor(s) : Arne Storjohann, Cory Fletcher
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


#include "latreduce.h"

/* Input: A -  an array of length n x m of initialized mpz_t
            
Note: A should have rank n
    
Output:  A is LLL reduced inplace

*/
void 
LLL (mpz_t *AA, int n, int m)
{
  long i, j, l, k, flag, kmax;
  mpz_t *BB, *DD, *dd, *d0, *d1, *d2, *c;
  static mpz_t t, t1, t2, t3, dn, q, r, swap;

  { mpz_init(t);  mpz_init(t1);  mpz_init(t2);  mpz_init(t3); }
  { mpz_init(dn);  mpz_init(q);  mpz_init(r);  mpz_init(swap); }
  BB = XCALLOC(mpz_t, n*n);
  DD = XCALLOC(mpz_t, n+1);
  dd = DD + 1;
  mpz_init_set_si(d(-1), 1);
  for (i = 0; i < n; i++) 
    {
      mpz_init(d(i));
      for (j = i + 1; j < n; j++)
	mpz_init(B(i, j));
    }
  for (l = 0; l < m; l++) 
    {
      mpz_mul(t1, A(0, l), A(0, l));
      mpz_add(d(0), d(0), t1);
    }
  k = 1;
  kmax = 0;
  flag = 1;
  while (k < n) 
    {
      if (k > kmax) 
	{
	  for (i = 0; i <= k; i++) 
	    {
	      mpz_set_ui(t, 0);
	      for (l = 0; l < m; l++) 
		{
		  mpz_mul(t1, A(i, l), A(k, l));
		  mpz_add(t, t, t1);
		}
	      for (l = 0; l < i; l++) 
		{
		  mpz_mul(t1, t, d(l));
		  mpz_mul(t2, B(l, i), B(l, k));
		  mpz_sub(t3, t1, t2);
		  mpz_divexact(t, t3, d(l - 1));
		}
	      if (i == k) { mpz_set(d(k), t); }
	      else { mpz_set(B(i, k), t); }
	    }
	  kmax++;
	}

      for (i = k - 1; i >= 0 && flag; i--) 
	{
	  mpz_fdiv_qr(q, r, B(i, k), d(i));
	  mpz_mul_2exp(t, r, 1);
	  if (mpz_cmp(t, d(i)) > 0) 
	    {
	      mpz_add_ui(q, q, 1);
	      mpz_sub(r, r, d(i));
	    }
	  if (!mpz_sgn(q)) { continue; }
	  for (j = 0; j < m; j++) 
	    {
	      mpz_mul(t, q, A(i, j));
	      mpz_sub(A(k, j), A(k, j), t);
	    }
	  for (j = 0; j < i; j++) 
	    {
	      mpz_mul(t, q, B(j, i));
	      mpz_sub(B(j, k), B(j, k), t);
	    }
	  mpz_set(B(i, k), r);
	}

      d0 = dd + k;
      d1 = dd + k - 1;
      d2 = dd + k - 2;
      c = BB + (k - 1) * n + k;
      mpz_mul(t, *d0, *d2);
      mpz_mul_2exp(t1, t, 1);
      mpz_mul(t2, *d1, *d1);
      if (mpz_cmp(t1, t2) < 0) 
	{
	  mpz_mul(t2, *c, *c);
	  mpz_add(t, t2, t);
	  mpz_divexact(dn, t, *d1);
	  for (j = 0; j < m; j++) 
	    {
	      mpz_set(swap, A(k, j));
	      mpz_set(A(k, j), A(k - 1, j));
	      mpz_set(A(k - 1, j), swap);
	    }
	  for (i = 0; i < k - 1; i++) 
	    {
	      mpz_set(swap, B(i, k - 1));
	      mpz_set(B(i, k - 1), B(i, k));
	      mpz_set(B(i, k), swap);
	    }
	  for (j = k + 1; j <= kmax; j++) 
	    {
	      mpz_set(t, B(k - 1, j));
	      mpz_mul(t1, B(k, j), *d2);
	      mpz_mul(t2, *c, B(k - 1, j));
	      mpz_add(t1, t2, t1);
	      mpz_divexact(B(k - 1, j), t1, *d1);
	      mpz_mul(t1, dn, t);
	      mpz_mul(t2, *c, B(k - 1, j));
	      mpz_sub(t1, t1, t2);
	      mpz_divexact(B(k, j), t1, *d2);
	    }
	  mpz_set(d(k - 1), dn);
	  if (k > 1) 
	    {
	      k--;
	      flag = 0;
	    }
	  else { flag = 1; }
	}
      else 
	{
	  k++;
	  flag = 1;
	}
    }

  mpz_clear(d(-1));
  for (i = 0; i < n; i++) 
    {
      mpz_clear(d(i));
      for (j = i + 1; j < n; j++) { mpz_clear(B(i, j)); }
    }
  { XFREE(BB); XFREE(DD); }
  { mpz_clear(t);  mpz_clear(t1);  mpz_clear(t2);  mpz_clear(t3); }
  { mpz_clear(dn);  mpz_clear(q);  mpz_clear(r);  mpz_clear(swap); }
}


mpz_t	mpz_t_tmp[NTMP];
long	mpz_t_ntmp = 0;


void ired(mpz_t *A, long n, long m, long n1)
{
  mpz_t	*T, *U, *B, d, *dd, dk2, dd2k, M, t1, t2;
  long	i, j, k;
#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  double tt;
#endif

#if HAVE_TIME_H
  clock_t ti, ti1;
#endif

  if (n < 2 || n1 < 1) return;
  mpz_initall_tmp();

  /* Phase 1: Fraction Free Gaussian Elimination */

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  ti = clock();
#endif

  /* set T = A x A^tr */
  T = XMALLOC(mpz_t, n * n);
  for (i = 0; i < n1; i++)
    for (j = 0; j < n1; j++)
      {
	mpz_init(T[i*n1+j]);
	for (k = 0; k < m; k++) 
	  mpz_addmul(T[i*n1+j], A[i*m+k], A[j*m+k]);
      }
  for (i = n1*n1; i < n*n; i++) mpz_init(T[i]);

  /* compute T = FF(A x A^tr) x (A x A^tr) */
  mpz_init_set_ui(d, 1);
  for (j = 0; j < n1; j++) {
    for (i = j+1; i < n1; i++) {
      for (k = j + 1; k < n1; k++) {
	mpz_mul(T[i*n1+k], T[j*n1+j], T[i*n1+k]);
	mpz_submul(T[i*n1+k], T[i*n1+j], T[j*n1+k]);
	mpz_divexact(T[i*n1+k], T[i*n1+k], d);
      }
      mpz_set_ui(T[i*n1+j], 0);
    }
    mpz_set(d, T[j*n1+j]);
  }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("                        time for first phase: %f\n", tt);
  ti = clock();
#endif

  /* Phase 2: 2-reduction */

  mpz_init(t1); mpz_init(t2); mpz_init(M);
  mpz_init(dk2); mpz_init(dd2k);
  dd = XMALLOC(mpz_t, n);
  for (i = 0; i < n; i++) mpz_init(dd[i]);
	
  U = XMALLOC(mpz_t, n * n);
  for (i = 0; i < n1; i++)
    for (j = 0; j < n1; j++)
      mpz_init_set_ui(U[i*n1+j], (i == j) ? 1 : 0);
  for (i = n1*n1; i < n*n; i++) mpz_init(U[i]);

  while (1) {

    UpdateMdd(M, dd, n1, T);

    /* find next k */
    k = 1;
    mpz_mul(dk2, T[(k-1)*n1+(k-1)], T[(k-1)*n1+(k-1)]);
    mpz_set(dd2k, T[k*n1+k]);

    for (i = 2; i < n1; i++) {
      mpz_mul(t1, T[i*n1+i], T[(i-2)*n1+(i-2)]); mpz_mul(t1, t1, dk2);
      mpz_mul(t2, T[(i-1)*n1+(i-1)], T[(i-1)*n1+(i-1)]); mpz_mul(t2, t2, dd2k);

      if (mpz_cmp(t1, t2) < 0) {
	k = i;
	mpz_mul(dk2, T[(k-1)*n1+(k-1)], T[(k-1)*n1+(k-1)]);
	mpz_mul(dd2k, T[k*n1+k], T[(k-2)*n1+(k-2)]);
      }
    }		


    /* check if finished */
    mpz_mul_2exp(t1, dd2k, 1);
    if (mpz_cmp(t1, dk2) >= 0 || n1 < 2) break;

    TwoReduce(U, T, n1, M, dd, k);
  }	
#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("                        time for phase2 : %f\n", tt);
  ti = clock();
#endif

  /* Phase 3: Size-reduction */

  for (k = 1; k < n1; k++)
    for (j = k-1; j >= 0; j--) {
      mpz_div_round(t1, T[j*n1+k], T[j*n1+j]);
      ModSubtractRow(U, T, n1, M, dd, k, j, t1);
    }

  /* set B = U x A */
  B = XMALLOC(mpz_t, n1 * m);
  for (i = 0; i < n1; i++)
    for (j = 0; j < m; j++) {
      mpz_init(B[i*m+j]);
      for (k = 0; k < n1; k++) 
	mpz_addmul(B[i*m+j], U[i*n1+k], A[k*m+j]);
      mpz_mods(B[i*m+j], B[i*m+j], M);
    }

  /* set T = B x B^tr */
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	mpz_set_ui(T[i*n+j], 0);
	for (k = 0; k < m; k++) 
	  mpz_addmul(T[i*n+j], 
		     (i >= n1) ? A[i*m+k] : B[i*m+k], 
		     (j >= n1) ? A[j*m+k] : B[j*m+k]);
      }

  /* compute T = FF(B x B^tr) x (B x B^tr) */
  mpz_set_ui(d, 1);
  for (j = 0; j < n; j++) {
    for (i = j+1; i < n; i++) {
      for (k = j + 1; k < n; k++) {
	mpz_mul(T[i*n+k], T[j*n+j], T[i*n+k]);
	mpz_submul(T[i*n+k], T[i*n+j], T[j*n+k]);
	mpz_divexact(T[i*n+k], T[i*n+k], d);
      }
      mpz_set_ui(T[i*n+j], 0);
    }
    mpz_set(d, T[j*n+j]);
  }

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      mpz_set_ui(U[i*n+j], (i == j) ? 1 : 0);

  UpdateMdd(M, dd, n, T);

  for (k = n1; k < n; k++)
    for (j = n1-1; j >= 0; j--) {
      mpz_div_round(t1, T[j*n+k], T[j*n+j]);
      ModSubtractRow(U, T, n, M, dd, k, j, t1);
    }

  /* Compute R2 = A2 + U2 x A1 (mods M) */

  for (i = n1; i < n; i++)
    for (j = 0; j < m; j++) {
      for (k = 0; k < n1; k++)
	mpz_addmul(A[i*m+j], U[i*n+k], B[k*m+j]);
      mpz_mods(A[i*m+j], A[i*m+j], M);
    }

   /* added stuff here */
  for (i = 0; i < n1; i++)
    for (j = 0; j < m; j++) 
      mpz_set(A[i*m+j],B[i*m+j]);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("                        time for  phase 3: %f\n", tt);
#endif

  /* Clear gmp variables */
	
  { mpz_clear(d); mpz_clear(M); mpz_clear(t1); mpz_clear(t2); }
  { mpz_clear(dk2); mpz_clear(dd2k); }
  for (i = 0; i < n*n; i++) mpz_clear(T[i]); XFREE(T);
  for (i = 0; i < n*n; i++) mpz_clear(U[i]); XFREE(U);
  for (i = 0; i < n1*m; i++) mpz_clear(B[i]); XFREE(B);
  for (i = 0; i < n; i++) mpz_clear(dd[i]); XFREE(dd);
  mpz_freeall_tmp();

  return;
}

void mpz_initall_tmp() { long i; for (i = 0; i < NTMP; i++) mpz_init(mpz_t_tmp[i]); }

void mpz_freeall_tmp() { long i; for (i = 0; i < NTMP; i++) mpz_clear(mpz_t_tmp[i]); }

void mpz_mods(mpz_t r, mpz_t n, mpz_t d)
{
  mpz_mod(r, n, d);
  mpz_tdiv_q_2exp(d, d, 1);
  if (mpz_cmp(r, d) > 0) {
    mpz_mul_2exp(d, d, 1);
    mpz_sub(r, r, d);
  } else 
    mpz_mul_2exp(d, d, 1);
}

void mpz_div_round(mpz_t r, mpz_t n, mpz_t d)
{
  long u = (mpz_sgn(n) < 0), v = (mpz_sgn(d) < 0); 
  if (u) mpz_neg(n, n); if (v) mpz_neg(d, d);
  mpz_set(r, d); 
  mpz_addmul_ui(r, n, 2);
  mpz_mul_2exp(d, d, 1);
  mpz_fdiv_q(r, r, d); if (u ^ v) mpz_neg(r, r);
  mpz_fdiv_q_2exp(d, d, 1);
  if (u) mpz_neg(n, n); if (v) mpz_neg(d, d);
}

void SubtractRow(mpz_t *A, mpz_t *T, long n, long k, long r, mpz_t q)
{
  long	i;

  for (i = 0; i < n; i++) {
    mpz_submul(A[k*n+i], A[r*n+i], q);
    mpz_submul(T[i*n+k], T[i*n+r], q);
  }
}

void SwitchRow(mpz_t *A, mpz_t *T, long n, long k)
{
  long	i;

  for (i = 0; i < n; i++) mpz_swap(A[k*n+i], A[(k-1)*n+i]);
  for (i = 0; i < n; i++) {
    if (k >= 2) mpz_mul(T[k*n+i], T[k*n+i], T[(k-2)*n+(k-2)]);
    mpz_addmul(T[k*n+i], T[(k-1)*n+k], T[(k-1)*n+i]);
    mpz_divexact(T[k*n+i], T[k*n+i], T[(k-1)*n+(k-1)]);
  }
  for (i = 0; i < n; i++) mpz_swap(T[k*n+i], T[(k-1)*n+i]);
  for (i = 0; i < n; i++) mpz_swap(T[i*n+k], T[i*n+(k-1)]);
  for (i = 0; i < n; i++) {
    mpz_mul(T[k*n+i], T[k*n+i], T[(k-1)*n+(k-1)]);
    mpz_submul(T[k*n+i], T[(k-1)*n+k], T[(k-1)*n+i]);
    if (k >= 2) mpz_divexact(T[k*n+i], T[k*n+i], T[(k-2)*n+(k-2)]);
  }		
}

void ModSubtractRow(mpz_t *A, mpz_t *T, long n, mpz_t M, mpz_t *dd, long k, 
                    long r, mpz_t q)
{
  long	i;

  SubtractRow(A, T, n, k, r, q);
  for (i = 0; i < k-1; i++) mpz_mods(T[i*n+k], T[i*n+k], dd[i]);
  for (i = 0; i < n; i++) mpz_mods(A[k*n+i], A[k*n+i], M);
}

void ModSwitchRow(mpz_t *A, mpz_t *T, long n, mpz_t M, mpz_t *dd, long k)
{
  long	i;

  SwitchRow(A, T, n, k);
  mpz_mul(dd[k], T[k*n+k], M); 
  mpz_mul(dd[k], dd[k], T[(k-1)*n+(k-1)]);
  mpz_mul(dd[k-1], T[(k-1)*n+(k-1)], M); 
  if (k >= 2) mpz_mul(dd[k-1], dd[k-1], T[(k-2)*n+(k-2)]);
  for (i = 0; i < k-2; i++) mpz_mods(T[i*n+(k-1)], T[i*n+(k-1)], dd[i]);
  for (i = 0; i < k-1; i++) mpz_mods(T[i*n+k], T[i*n+k], dd[i]);
  for (i = k; i < n; i++) mpz_mods(T[(k-1)*n+i], T[(k-1)*n+i], dd[k-1]);
  for (i = k+1; i < n; i++) mpz_mods(T[k*n+i], T[k*n+i], dd[k]);
}

void GetNextU(mpz_t *U[], mpz_t t00, mpz_t t11, mpz_t t12, mpz_t t22)
{
  mpz_t	*s00, *s11, *s12, *s22, *q, *a, *b;

  s00 = mpz_next_tmp(); mpz_set(*s00, t00);
  s11 = mpz_next_tmp(); mpz_set(*s11, t11);
  s12 = mpz_next_tmp(); mpz_set(*s12, t12);
  s22 = mpz_next_tmp(); mpz_set(*s22, t22);
  q = mpz_next_tmp(); a = mpz_next_tmp(); b = mpz_next_tmp();

  mpz_set_ui(*U[0], 1); mpz_set_ui(*U[1], 0);
  mpz_set_ui(*U[2], 0); mpz_set_ui(*U[3], 1);

  while (1) {
    mpz_mul(*a, *s22, *s00); mpz_mul_2exp(*a, *a, 1);
    mpz_mul(*b, *s11, *s11);
    if (mpz_cmp(*a, *b) >= 0) break;
    mpz_tdiv_q_2exp(*a, *a, 1);

    mpz_div_round(*q, *s12, *s11);
    mpz_submul(*U[2], *q, *U[0]); mpz_submul(*U[3], *q, *U[1]);
    mpz_swap(*U[0], *U[2]); mpz_swap(*U[1], *U[3]);
    mpz_submul(*s12, *q, *s11);
    mpz_addmul(*a, *s12, *s12); mpz_divexact(*s11, *a, *s11);
  }

  mpz_free_tmp(7);	
}

void TwoReduce(mpz_t *A, mpz_t *T, long n, mpz_t M, mpz_t *dd, long k)
{
  mpz_t	*U[4], *a, *t;
  long	i;

  t = mpz_next_tmp(); 
  for (i = 0; i < 4; i++) U[i] = mpz_next_tmp();
  a = mpz_next_tmp(); (k >= 2) ? mpz_set(*a, T[(k-2)*n+(k-2)]) : mpz_set_ui(*a, 1);

  GetNextU(U, *a, T[(k-1)*n+(k-1)], T[(k-1)*n+k], T[k*n+k]);

  /* update A */
  for (i = 0; i < n; i++) {
    mpz_set(*t, A[(k-1)*n+i]);
    mpz_mul(A[(k-1)*n+i], A[(k-1)*n+i], *U[0]);
    mpz_addmul(A[(k-1)*n+i], *U[1], A[k*n+i]); 
    mpz_mods(A[(k-1)*n+i], A[(k-1)*n+i], M);
    mpz_mul(A[k*n+i], A[k*n+i], *U[3]);
    mpz_addmul(A[k*n+i], *U[2], *t);	
    mpz_mods(A[k*n+i], A[k*n+i], M);
  }

  /* update T */
  for (i = k-1; i < n; i++) {
    mpz_mul(T[k*n+i], T[k*n+i], *a);
    mpz_addmul(T[k*n+i], T[(k-1)*n+k], T[(k-1)*n+i]);
    mpz_divexact(T[k*n+i], T[k*n+i], T[(k-1)*n+(k-1)]);
  }
  for (i = k-1; i < n; i++) {
    mpz_set(*t, T[(k-1)*n+i]);
    mpz_mul(T[(k-1)*n+i], T[(k-1)*n+i], *U[0]);
    mpz_addmul(T[(k-1)*n+i], *U[1], T[k*n+i]); 
    mpz_mul(T[k*n+i], T[k*n+i], *U[3]);
    mpz_addmul(T[k*n+i], *U[2], *t);	
  }
  for (i = 0; i <= k; i++) {
    mpz_set(*t, T[i*n+k-1]);
    mpz_mul(T[i*n+k-1], T[i*n+k-1], *U[0]);
    mpz_addmul(T[i*n+k-1], *U[1], T[i*n+k]); 
    mpz_mul(T[i*n+k], T[i*n+k], *U[3]);
    mpz_addmul(T[i*n+k], *U[2], *t);	
  }
  for (i = k-1; i < n; i++) {
    mpz_mul(T[k*n+i], T[k*n+i], T[(k-1)*n+(k-1)]);
    mpz_submul(T[k*n+i], T[(k-1)*n+k], T[(k-1)*n+i]);
    mpz_divexact(T[k*n+i], T[k*n+i], *a);
  }
  mpz_mul(dd[k], T[k*n+k], M); 
  mpz_mul(dd[k], dd[k], T[(k-1)*n+(k-1)]);
  mpz_mul(dd[k-1], T[(k-1)*n+(k-1)], M); 
  if (k >= 2) mpz_mul(dd[k-1], dd[k-1], T[(k-2)*n+(k-2)]);
  for (i = 0; i < k-2; i++) mpz_mods(T[i*n+(k-1)], T[i*n+(k-1)], dd[i]);
  for (i = 0; i < k-1; i++) mpz_mods(T[i*n+k], T[i*n+k], dd[i]);
  for (i = k; i < n; i++) mpz_mods(T[(k-1)*n+i], T[(k-1)*n+i], dd[k-1]);
  for (i = k+1; i < n; i++) mpz_mods(T[k*n+i], T[k*n+i], dd[k]);

  mpz_free_tmp(6);
}

void UpdateMdd(mpz_t M, mpz_t *dd, long n, mpz_t *T)
{
  long	i, j;
  mpz_t	*t1, *t2;
  t1 = mpz_next_tmp(); t2 = mpz_next_tmp();

  /* compute M */
  mpz_set(*t1, T[0]);
  for (i = 1; i < n; i++) {
    mpz_cdiv_q(*t2, T[i*n+i], T[(i-1)*n+(i-1)]);
    if (mpz_cmp(*t2, *t1) > 0) mpz_set(*t1, *t2);
  }
  mpz_mul_ui(*t1, *t1, n);
		
  for (j = 1; mpz_cmp_ui(*t1, 1) > 0; j++) mpz_fdiv_q_2exp(*t1, *t1, 1);
  j = j / 2 + 5;
		
  mpz_set_ui(M, 1); mpz_mul_2exp(M, M, j);

  /* precompute dd[i] = d[i]*d[i-1]*M */
  mpz_mul(dd[0], T[0], M);
  for (i = 1; i < n; i++) {
    mpz_mul(dd[i], T[(i-1)*n+(i-1)], T[i*n+i]);
    mpz_mul(dd[i], dd[i], M);
  }

  mpz_free_tmp(2);
}
