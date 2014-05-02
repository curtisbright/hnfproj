/* ---------------------------------------------------------------------
 *
 * -- Integer Matrix Library (IML)
 *    (C) Copyright 2004, 2006 All Rights Reserved
 *
 * -- IML routines -- Version 1.0.1 -- November, 2004
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


#include "nullspace.h"

/*
 * Calling Sequence:
 *   nullspaceLong(n, m, A, mp_N_pass)
 *
 * Summary:
 *   Compute the right nullspace of A.
 *
 * Input:  n: long, row dimension of A
 *         m: long, column dimension of A
 *         A: 1-dim signed long array length n*m, representing n x m matrix
 *            in row major order
 *
 * Output:
 *   - *mp_N_pass: points to a 1-dim mpz_t array of length m*s, where s is the 
 *                dimension of the right nullspace of A
 *   - the dimension s of the nullspace is returned
 *
 * Notes:
 *   - The matrix A is represented by one-dimension array in row major order.
 *   - Space for what mp_N_points to is allocated by this procedure: if the
 *     nullspace is empty, mp_N_pass is set to NULL.
 */
long
nullspaceLong(const long n, const long m, const long *A, mpz_t * *mp_N_pass)
{
    long i, j, k, r, s, *P, *rp, *Pt, *rpt, *C, flag, temp;
    double *DA;
    FiniteField p, d;
    mpz_t *mp_B, *mp_N, mp_D, mp_t1, mp_t2;

    P = XCALLOC(long, n + 1);
    rp = XCALLOC(long, n + 1);
    while (1) {
	p = RandPrime(15, 19);
	DA = XCALLOC(double, n * m);
	for (i = 0; i < n * m; i++)
	    DA[i] =
		(double) ((temp =
			   (A[i] % ((long) p))) >=
			  0 ? temp : ((long) p) + temp);
	for (i = 0; i < n + 1; i++) {
	    P[i] = i;
	    rp[i] = 0;
	}
	d = 1;
	RowEchelonTransform(p, DA, n, m, 1, 1, 0, 0, P, rp, &d);
	XFREE(DA);
	r = rp[0];
	s = m - r;
	if (s == 0) {
	    *mp_N_pass = NULL;
	} else if (r == 0) {
	    flag = 1;
	    for (i = 0; i < n * m; i++)
		if (A[i] != 0)
		    flag = 0;
	    if (!flag)
		continue;
	    mp_N = XCALLOC(mpz_t, m * m);
	    for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++)
		    mpz_init_set_ui(mp_N[i * m + j], 0);
		mpz_init_set_ui(mp_N[i * m + i], 1);
	    }
	    *mp_N_pass = mp_N;
	} else {		/* r>0 and s>0 */

	    Pt = revseq(r, n, P);
	    rpt = revseq(r, m, rp);

	    C = XCALLOC(long, r * r);
	    for (i = 0; i < r; i++)
		for (j = 0; j < r; j++)
		    C[i * r + j] = A[Pt[i] * m + rpt[j]];

	    mp_B = XCALLOC(mpz_t, r * s);
	    for (i = 0; i < r; i++)
		for (j = 0; j < s; j++)
		    mpz_init_set_si(mp_B[i * s + j],
				    A[Pt[i] * m + rpt[r + j]]);

	    mpz_init(mp_D);
	    mp_N = XCALLOC(mpz_t, m * s);
	    for (i = 0; i < m * s; i++)
		mpz_init(mp_N[i]);

	    nonsingSolvMM(RightSolu, r, s, C, mp_B, mp_N, mp_D);
	    mpz_neg(mp_D, mp_D);
	    for (i = 0; i < s; i++)
		mpz_set(mp_N[(r + i) * s + i], mp_D);

	    XFREE(C);
	    for (i = 0; i < r * s; i++)
		mpz_clear(mp_B[i]);
	    XFREE(mp_B);
	    mpz_clear(mp_D);

	    for (i = r; i >= 1; i--)
		for (j = 0; j < s; j++)
		    mpz_swap(mp_N[(i - 1) * s + j],
			     mp_N[(rp[i] - 1) * s + j]);

	    *mp_N_pass = mp_N;

	    flag = 1;
	    mpz_init(mp_t1);
	    mpz_init(mp_t2);
	    for (i = r; i < n && flag; i++) {
		for (j = 0; j < s && flag; j++) {
		    mpz_set_ui(mp_t2, 0);
		    for (k = 0; k < m; k++) {
			mpz_mul_si(mp_t1, mp_N[k * s + j],
				   A[Pt[i] * m + k]);
			mpz_add(mp_t2, mp_t2, mp_t1);
		    }
		    if (mpz_sgn(mp_t2) != 0)
			flag = 0;
		}
	    }
	    mpz_clear(mp_t1);
	    mpz_clear(mp_t2);

	    XFREE(Pt);
	    XFREE(rpt);

	    if (!flag) {
		for (i = 0; i < m * s; i++)
		    mpz_clear(mp_N[i]);
		XFREE(mp_N);
		continue;
	    }
	}
	break;
    }
    XFREE(P);
    XFREE(rp);

    return s;

}


/*
 * Calling Sequence:
 *   nullspaceMP(n, m, A, mp_N_pass)
 *
 * Summary: Compute the right nullspace of A. In this function A is a
 *          1-dimensional mpz_t array.
 *
 * Input:  n: long, row dimension of A
 *         m: long, column dimension of A
 *         A: 1-dim mpz_t array length n*m, representing n x m matrix
 *            in row major order
 *
 * Output:
 *   - *mp_N_pass: points to a 1-dim mpz_t array of length m*s, where s is the 
 *                dimension of the right nullspace of A
 *   - the dimension s of the nullspace is returned
 *
 * Notes:
 *   - The matrix A is represented by one-dimension array in row major order.
 *   - Space for what mp_N_points to is allocated by this procedure: if the
 *     nullspace is empty, mp_N_pass is set to NULL.
 */

long
nullspaceMP(const long n, const long m, const mpz_t *A, mpz_t * *mp_N_pass)
{
    long i, j, k, r, s, *P, *rp, *Pt, *rpt, flag;
    double *DA;
    FiniteField p, d = 1;
    mpz_t *mp_B, *mp_N, mp_D, mp_t1, mp_t2, *C, mp_r;

    mpz_init(mp_r);

    P = XCALLOC(long, n + 1);
    rp = XCALLOC(long, n + 1);
    while (1) {
	p = RandPrime(15, 19);
	DA = XCALLOC(double, n * m);
	for (i = 0; i < n * m; i++) {
  	    mpz_mod_ui (mp_r, A[i], p);
	    DA[i] = mpz_get_d(mp_r);
	}
	for (i = 0; i < n + 1; i++) {
	    P[i] = i;
	    rp[i] = 0;
	}
	d = 1;
	RowEchelonTransform(p, DA, n, m, 1, 1, 0, 0, P, rp, &d);
	XFREE(DA);
	r = rp[0];
	s = m - r;
	if (s == 0) {
	    *mp_N_pass = NULL;
	} else if (r == 0) {
	    flag = 1;
	    for (i = 0; i < n * m; i++)
		if ( mpz_cmp_si(A[i],0) )
		    flag = 0;
	    if (!flag)
		continue;
	    mp_N = XCALLOC(mpz_t, m * m);
	    for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++)
		    mpz_init_set_ui(mp_N[i * m + j], 0);
		mpz_init_set_ui(mp_N[i * m + i], 1);
	    }
	    *mp_N_pass = mp_N;
	} else {		/* r>0 and s>0 */

	    Pt = revseq(r, n, P);
	    rpt = revseq(r, m, rp);

	    C = XCALLOC(mpz_t, r * r);
	    for (i = 0; i < r; i++)
		for (j = 0; j < r; j++)
		    mpz_init_set(C[i * r + j], A[Pt[i] * m + rpt[j]]);

	    mp_B = XCALLOC(mpz_t, r * s);
	    for (i = 0; i < r; i++)
		for (j = 0; j < s; j++)
		    mpz_init_set(mp_B[i * s + j], A[Pt[i] * m + rpt[r + j]]);

	    mpz_init(mp_D);
	    mp_N = XCALLOC(mpz_t, m * s);
	    for (i = 0; i < m * s; i++)
		mpz_init(mp_N[i]);


	    nonsingSolvLlhsMM(RightSolu, r, s, C, mp_B, mp_N, mp_D);

	    mpz_neg(mp_D, mp_D);
	    for (i = 0; i < s; i++)
		mpz_set(mp_N[(r + i) * s + i], mp_D);

	    for (i = 0; i < r*r; i++)
	        mpz_clear(C[i]);
	    XFREE(C);

	    for (i = 0; i < r * s; i++)
		mpz_clear(mp_B[i]);
	    XFREE(mp_B);
	    mpz_clear(mp_D);

	    for (i = r; i >= 1; i--)
		for (j = 0; j < s; j++)
		    mpz_swap(mp_N[(i - 1) * s + j],
			     mp_N[(rp[i] - 1) * s + j]);

	    *mp_N_pass = mp_N;

	    flag = 1;
	    mpz_init(mp_t1);
	    mpz_init(mp_t2);
	    for (i = r; i < n && flag; i++) {
		for (j = 0; j < s && flag; j++) {
		    mpz_set_ui(mp_t2, 0);
		    for (k = 0; k < m; k++) {
			mpz_mul(mp_t1, mp_N[k * s + j],  A[Pt[i] * m + k]);
			mpz_add(mp_t2, mp_t2, mp_t1);
		    }
		    if (mpz_sgn(mp_t2) != 0)
			flag = 0;
		}
	    }
	    mpz_clear(mp_t1);
	    mpz_clear(mp_t2);

	    XFREE(Pt);
	    XFREE(rpt);

	    if (!flag) {
		for (i = 0; i < m * s; i++)
		    mpz_clear(mp_N[i]);
		XFREE(mp_N);
		continue;
	    }
	}
	break;
    }
    XFREE(P);
    XFREE(rp);

    mpz_clear(mp_r);
    return s;

}


/*
 * Calling Sequence:
 *   kernelLong(n, m, A, mp_N_pass, reduce)
 *
 * Summary:
 *   Compute the right kernel of A. 
 *   If reduce is true the kernel will be reduced using LLL.
 *
 * Input:  n: long, row dimension of A
 *         m: long, column dimension of A
 *         A: 1-dim signed long array length n*m, representing n x m matrix
 *            in row major order
 *
 * Output:
 *   - *mp_N_pass: points to a 1-dim mpz_t array of length m*s, where s is the 
 *                dimension of the right kernel of A
 *   - the dimension s of the kernel is returned
 *
 * Notes:
 *   - The matrix A is represented by one-dimension array in row major order.
 *   - Space for what mp_N_points to is allocated by this procedure: if the
 *     nullspace is empty, mp_N_pass is set to NULL.
 */
long
kernelLong(const long n, const long m, const long *A, mpz_t * *mp_N_pass, const long reduce)
{
    long i, j, k, r, s, *P, *rp, *Pt, *rpt, *C, flag, temp;
    double *DA;
    FiniteField p, d;
    mpz_t *mp_B, *mp_N, *mp_M, mp_D, mp_t1, mp_t2;

    P = XCALLOC(long, n + 1);
    rp = XCALLOC(long, n + 1);
    while (1) {
	p = RandPrime(15, 19);
	DA = XCALLOC(double, n * m);
	for (i = 0; i < n * m; i++)
	    DA[i] =
		(double) ((temp =
			   (A[i] % ((long) p))) >=
			  0 ? temp : ((long) p) + temp);
	for (i = 0; i < n + 1; i++) {
	    P[i] = i;
	    rp[i] = 0;
	}
	d = 1;
	RowEchelonTransform(p, DA, n, m, 1, 1, 0, 0, P, rp, &d);
	XFREE(DA);
	r = rp[0];
	s = m - r;
	if (s == 0) {
	    *mp_N_pass = NULL;
	} else if (r == 0) {
	    flag = 1;
	    for (i = 0; i < n * m; i++)
		if (A[i] != 0)
		    flag = 0;
	    if (!flag)
		continue;
	    mp_N = XCALLOC(mpz_t, m * m);
	    for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++)
		    mpz_init_set_ui(mp_N[i * m + j], 0);
		mpz_init_set_ui(mp_N[i * m + i], 1);
	    }
	    *mp_N_pass = mp_N;
	} else {		/* r>0 and s>0 */

	    Pt = revseq(r, n, P);
	    rpt = revseq(r, m, rp);

	    C = XCALLOC(long, r * r);
	    for (i = 0; i < r; i++)
		for (j = 0; j < r; j++)
		    C[i * r + j] = A[Pt[i] * m + rpt[j]];

	    mp_B = XCALLOC(mpz_t, r * s);
	    for (i = 0; i < r; i++)
		for (j = 0; j < s; j++)
		    mpz_init_set_si(mp_B[i * s + j],
				    A[Pt[i] * m + rpt[r + j]]);

	    mpz_init(mp_D);
	    mp_N = XCALLOC(mpz_t, m * s);
	    for (i = 0; i < m * s; i++)
		mpz_init(mp_N[i]);

	    nonsingSolvMM(RightSolu, r, s, C, mp_B, mp_N, mp_D);
	    mpz_neg(mp_D, mp_D);
	    for (i = 0; i < s; i++)
		mpz_set(mp_N[(r + i) * s + i], mp_D);

	    XFREE(C);
	    for (i = 0; i < r * s; i++)
		mpz_clear(mp_B[i]);
	    XFREE(mp_B);
	    mpz_clear(mp_D);

	    /* Convert nullspace to kernel */
	    mp_M = XCALLOC(mpz_t, m * s);
	    for (i = 0; i < m * s; i++)
		mpz_init_set(mp_M[i], mp_N[i]);

	    specialHermite(0, m-s, s, 0, mp_M, NULL);

	    kernelBasis(m-s, s, 0, mp_M, mp_N);

	    if (reduce==1)  {
               for (i=0; i<m; i++)
	           for (j=0; j<s; j++)
		       mpz_set(mp_M[j*m+i], mp_N[i*s+j]);

	       ired(mp_M, s, m, s);

               for (i=0; i<m; i++)
	           for (j=0; j<s; j++)
		       mpz_set(mp_N[i*s+j], mp_M[j*m+i]);
            }

	    for (i = 0; i < m * s; i++)
		mpz_clear(mp_M[i]);
	    XFREE(mp_M);

	    for (i = r; i >= 1; i--)
		for (j = 0; j < s; j++)
		    mpz_swap(mp_N[(i - 1) * s + j],
			     mp_N[(rp[i] - 1) * s + j]);

	    *mp_N_pass = mp_N;

	    flag = 1;
	    mpz_init(mp_t1);
	    mpz_init(mp_t2);
	    for (i = r; i < n && flag; i++) {
		for (j = 0; j < s && flag; j++) {
		    mpz_set_ui(mp_t2, 0);
		    for (k = 0; k < m; k++) {
			mpz_mul_si(mp_t1, mp_N[k * s + j],
				   A[Pt[i] * m + k]);
			mpz_add(mp_t2, mp_t2, mp_t1);
		    }
		    if (mpz_sgn(mp_t2) != 0)
			flag = 0;
		}
	    }
	    mpz_clear(mp_t1);
	    mpz_clear(mp_t2);

	    XFREE(Pt);
	    XFREE(rpt);

	    if (!flag) {
		for (i = 0; i < m * s; i++)
		    mpz_clear(mp_N[i]);
		XFREE(mp_N);
		continue;
	    }
	}
	break;
    }
    XFREE(P);
    XFREE(rp);

    return s;

}
