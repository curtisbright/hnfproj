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



#include "mtrans.h"

/*
 * Calling Sequence:
 *   RowEchelonTransform(p, A, n, m, frows, lrows, redflag, eterm, Q, rp, d)
 *
 * Summary:
 *   Compute a mod p row-echelon transform of a mod p input matrix
 *
 * Description:
 *   Given a n x m mod p matrix A, a row-echelon transform of A is a 4-tuple 
 *   (U,P,rp,d) with rp the rank profile of A (the unique and strictly 
 *   increasing list [j1,j2,...jr] of column indices of the row-echelon form 
 *   which contain the pivots), P a permutation matrix such that all r leading
 *   submatrices of (PA)[0..r-1,rp] are nonsingular, U a nonsingular matrix 
 *   such that UPA is in row-echelon form, and d the determinant of 
 *   (PA)[0..r-1,rp].
 *
 *   Generally, it is required that p be a prime, as inverses are needed, but
 *   in some cases it is possible to obtain an echelon transform when p is 
 *   composite. For the cases where the echelon transform cannot be obtained
 *   for p composite, the function returns an error indicating that p is 
 *   composite.
 *
 *   The matrix U is structured, and has last n-r columns equal to the last n-r
 *   columns of the identity matrix, n the row dimension of A.
 *
 *   The first r rows of UPA comprise a basis in echelon form for the row 
 *   space of A, while the last n-r rows of U comprise a basis for the left 
 *   nullspace of PA.
 *
 *   For efficiency, this function does not output an echelon transform 
 *   (U,P,rp,d) directly, but rather the expression sequence (Q,rp,d).
 *   Q, rp, d are the form of arrays and pointers in order to operate inplace,
 *   which require to preallocate spaces and initialize them. Initially, 
 *   Q[i] = i (i=0..n), rp[i] = 0 (i=0..n), and *d = 1. Upon completion, rp[0]
 *   stores the rank r, rp[1..r] stores the rank profile. i<=Q[i]<=n for 
 *   i=1..r. The input Matrix A is modified inplace and used to store U. 
 *   Let A' denote the state of A on completion. Then U is obtained from the
 *   identity matrix by replacing the first r columns with those of A', and P
 *   is obtained from the identity matrix by swapping row i with row Q[i], for
 *   i=1..r in succession.
 *
 *   Parameters flrows, lrows, redflag, eterm control the specific operations
 *   this function will perform. Let (U,P,rp,d) be as constructed above. If 
 *   frows=0, the first r rows of U will not be correct. If lrows=0, the last
 *   n-r rows of U will not be correct. The computation can be up to four 
 *   times faster if these flags are set to 0.
 *
 *   If redflag=1, the row-echelon form is reduced, that is (UPA)[0..r-1,rp] 
 *   will be the identity matrix. If redflag=0, the row-echelon form will not
 *   be reduced, that is (UPA)[1..r,rp] will be upper triangular and U is unit
 *   lower triangular. If frows=0 then redflag has no effect.
 *
 *   If eterm=1, then early termination is triggered if a column of the 
 *   input matrix is discovered that is linearly dependant on the previous
 *   columns. In case of early termination, the third return value d will be 0
 *   and the remaining components of the echelon transform will not be correct.
 *
 * Input:
 *         p: FiniteField, modulus
 *         A: 1-dim Double array length n*m, representation of a n x m input
 *            matrix
 *         n: long, row dimension of A
 *         m: long, column dimension of A
 *     frows: 1/0, 
 *          - if frows = 1, the first r rows of U will be correct
 *          - if frows = 0, the first r rows of U will not be correct
 *     lrows: 1/0,
 *          - if lrows = 1, the last n-r rows of U will be correct
 *          - if lrows = 0, the last n-r rows of U will not be correct
 *   redflag: 1/0,
 *          - if redflag = 1, compute row-echelon form
 *          - if redflag = 0, not compute reow-echelon form
 *     eterm: 1/0,
 *          - if eterm = 1, terminate early if not in full rank
 *          - if eterm = 0, not terminate early
 *         Q: 1-dim long array length n+1, compact representation of 
 *            permutation vector, initially Q[i] = i, 0 <= i <= n
 *        rp: 1-dim long array length n+1, representation of rank profile, 
 *            initially rp[i] = 0, 0 <= i <= n
 *         d: pointer to FiniteField, storing determinant of the matrix, 
 *            initially *d = 1
 *
 * Precondition:
 *   ceil(n/2)*(p-1)^2+(p-1) <= 2^53-1 = 9007199254740991 (n >= 2)
 *
 */

void 
RowEchelonTransform (const FiniteField p, Double *A, const long n, \
		     const long m, const long frows, const long lrows, \
		     const long redflag, const long eterm, long *Q, \
		     long *rp, FiniteField *d)
{
  mpz_t mp_a, mp_p; /* preallocated temporary storage */

  { mpz_init(mp_a); mpz_init_set_ui(mp_p, p); }

  RowEchelonTransform_rec(p, A, n, m, 1, m, 0, 0, frows, lrows, redflag, \
			  eterm, Q, rp, d, mp_a, mp_p);

  { mpz_clear(mp_a); mpz_clear(mp_p); }
}


long 
RowEchelonTransform_rec (const FiniteField p, Double *A, const long n, \
			 const long m, long m1, long m2, long k, \
			 const long ks, long frows, long lrows, long redflag, \
			 long eterm, long *P, long *rp, FiniteField *d, \
			 mpz_t mp_a, mpz_t mp_p)
{
  long i, j, r1, r2, r, ri, mm, inv;
  FiniteField a;
  Double *A1;
  Double b;

  if (m1 == m2)
    {
      for (i = k+1; i <= n; i++) { if (*(A+(i-1)*m+m1-1) != 0) { break; } }
      if ((i > n) && (eterm == 0)) { return 0; }
      else if ((i > n) && (eterm == 1))
	{
	  *d = 0;
	  return 0;
	}
      if (i > k+1)
	cblas_dswap(m-m1+1, A+k*m+m1-1, 1, A+(i-1)*m+m1-1, 1);
      if (k-ks > 0)
	cblas_dswap(k-ks, A+k*m, 1, A+(i-1)*m, 1);
      P[k+1] = i;
      mpz_set_d(mp_a, *(A+m*k+m1-1));
      inv = mpz_invert(mp_a, mp_a, mp_p);

      /* mp_a is not relatively prime with modulus */
      if (!inv) { iml_fatal("in RowEchelonTransform: modulus is composite"); }
      a = mpz_get_ui(mp_a);
      b = fmod(*(A+k*m+m1-1), p);
      if (b < 0) { b = p+b; }
      if ((frows == 1) && (redflag == 1))
	{
	  for (j = 1; j <= n; j++)
	    *(A+(j-1)*m+k-ks) = *(A+(j-1)*m+m1-1)*(p-a);
	  Dmod(p, A+k-ks, n, 1, m);
	  *(A+k*m+k-ks) = a;	  
	}
      else
	{
	  if (k+2 <= n) 
	    {
	      for (j = k+2; j <= n; j++)
		*(A+(j-1)*m+k-ks) = *(A+(j-1)*m+m1-1)*(p-a);
	      Dmod(p, A+(k+1)*m+k-ks, n-k-1, 1, m);
	    }
	  if (k > 0)
	    {
	      for (j = 1; j <= k; j++)
		*(A+(j-1)*m+k-ks) = 0;
	    }
	  *(A+k*m+k-ks) = 1;	  
	}
      ++rp[0];
      *d = fmod( (*d)*b, p);
      ri = rp[0];
      rp[ri] = m1;
      return 1;
    }
  /* Recursively solve the first subproblem */
  mm = m1+(long)((m2-m1)/2);
  r1 = RowEchelonTransform_rec(p, A, n, m, m1, mm, k, ks, frows, 1, redflag,\
			       eterm, P, rp, d, mp_a, mp_p);
  if ((eterm == 1) && (r1 < mm-m1+1))
    {
      *d = 0;
      return 0;
    }
  /* If r1=0 then don't need to construct second subproblem */
  if (r1 > 0)
    {
      /* Compute U1.A2 by submatrix multiply */
      if (k+r1 < n)
	{
	  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n-k-r1, \
		      m2-mm, r1, 1.0, A+(k+r1)*m+k-ks, m, A+k*m+mm, m, \
		      1.0, A+(k+r1)*m+mm, m);
	  Dmod(p, A+(k+r1)*m+mm, n-k-r1, m2-mm, m);
	}
      if ((frows == 1) && (redflag == 1))
	{
	  if (k > 0)
	    {
	      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, k, \
			  m2-mm, r1, 1.0, A+k-ks, m, A+k*m+mm, m, 1.0, \
			  A+mm, m);
	      Dmod(p, A+mm, k, m2-mm, m);
	    }
	  A1 = XMALLOC(Double, r1*(m2-mm));
	  DCopy(r1, m2-mm, A+k*m+mm, m, A1, m2-mm);
	  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, r1, m2-mm, \
		      r1, 1.0, A+k*m+k-ks, m, A1, m2-mm, 0.0, A+k*m+mm, m);
	  XFREE(A1);
	  Dmod(p, A+k*m+mm, r1, m2-mm, m);
	}
    }
  /* Recursively solve the second subproblem */
  r2 = RowEchelonTransform_rec(p, A, n, m, mm+1, m2, k+r1, ks, frows, lrows, \
			       redflag, eterm, P, rp, d, mp_a, mp_p);
  r = r1+r2;
  if ((eterm == 1) && (r < m2-m1+1))
    {
      *d = 0;
      return 0;
    }
  /* Only need to combine when both subproblems nontrivial */
  if ((r2 > 0) && (r1 > 0))
    {
      if ((k+r+1 <= n) && (lrows == 1))
	{
	  /* Bottom block of U */
	  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n-k-r, \
		      r1, r-r1, 1.0, A+(k+r)*m+k-ks+r1, m, A+(k+r1)*m+k-ks, \
		      m, 1.0, A+(k+r)*m+k-ks, m);
	  Dmod(p, A+(k+r)*m+k-ks, n-k-r, r1, m);
	}
      if (frows == 1)
	{
	  if (redflag == 1)
	    i = 1;
	  else
	    i = k+1;

	  /* Rows i..k of top block of U plus first r1 rows of middle block */
	  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, k+r1-i+1, \
		      r1, r-r1, 1.0, A+(i-1)*m+k-ks+r1, m, A+(k+r1)*m+k-ks, \
		      m, 1.0, A+(i-1)*m+k-ks, m);
	  Dmod(p, A+(i-1)*m+k-ks, k+r1-i+1, r1, m);

	  /* Last r2 rows of middle block */
	  A1 = XMALLOC(Double, r1*(r-r1));
	  DCopy(r-r1, r1, A+(k+r1)*m+k-ks, m, A1, r1);
	  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, r-r1, r1, \
		      r-r1, 1.0, A+(k+r1)*m+k-ks+r1, m, A1, r1, \
		      0.0, A+(k+r1)*m+k-ks, m);
	  XFREE(A1);
	  Dmod(p, A+(k+r1)*m+k-ks, r-r1, r1, m);
	}
    }
  return r;
}



/* 
 * Calling Sequence:
 *   Adj <-- mAdjoint(p, A, n)
 *
 * Summary:
 *   Compute the adjoint of a mod p square matrix
 *  
 * Description:
 *   Given a n x n mod p matrix A, the function computes adjoint of A. Input
 *   A is not modified upon completion.
 *
 * Input:
 *   p: FiniteField, prime modulus
 *      if p is a composite number, the routine will still work if no error 
 *      message is returned
 *   A: 1-dim Double array length n*n, representation of a n x n mod p matrix.
 *      The entries of A are casted from integers
 *   n: long, dimension of A
 *
 * Return:
 *   1-dim Double matrix length n*n, repesentation of a n x n mod p matrix,
 *   adjoint of A
 *
 * Precondition:
 *   n*(p-1)^2 <= 2^53-1 = 9007199254740991
 *
 */

Double *
mAdjoint (const FiniteField p, Double *A, const long n)
{
  long i, j, k, r, count=0;
  long *P, *rp;
  FiniteField det, d[1];
  Double p1, *B, *C;
 
  p1 = (Double)p;
  P = XMALLOC(long, n+1);
  for (i = 0; i < n+1; i++) { P[i] = i; }
  rp = XCALLOC(long, n+1);
  d[0] = 1;
  B = XMALLOC(Double, n*n);
  DCopy(n, n, A, n, B, n);
  RowEchelonTransform(p, B, n, n, 1, 1, 1, 0, P, rp, d);
  det = d[0];
  r = rp[0];
  if (r < n-1)
    {
      for (i = 0; i < n*n; i++) { B[i] = 0; }
      { XFREE(P); XFREE(rp); }
      return B;
    }
  if (r == n)
    {
      for (i = r; i > 0; i--)
	{
	  if (P[i] != i)
	    {
	      ++count;
	      cblas_dswap(n, B+i-1, n, B+P[i]-1, n);
	    }
	}
      if (count % 2 == 0)
	for (i = 0; i < n*n; i++) { B[i] = fmod(det*B[i], p1); }
      else
	for (i = 0; i < n*n; i++) { B[i] = fmod((p-det)*B[i], p1); }	
      { XFREE(P); XFREE(rp); }
      return B;
    }
  else
    {
      if (n == 1) 
	{ 
	  B[0] = 1; 
	  { XFREE(P); XFREE(rp); }
	  return B;
	}
      for (i = 0; i < n; i++) { B[i*n+n-1] = 0; }
      B[(n-1)*n+(n-1)] = 1;
      for (i = r; i > 0; i--)
	{
	  if (P[i] != i)
	    {
	      ++count;
	      cblas_dswap(n, B+i-1, n, B+P[i]-1, n);
	    }
	}
      for (j = 1; j < r+1; j++) { if (j != rp[j]) { break; } }
      C = XMALLOC(Double, n);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, B, n, \
		  A+j-1, n, 0.0, C, 1);
      Dmod(p, C, 1, n, 1);
      for (i = 0; i < n-1; i++) { C[i] = fmod((p-1)*C[i], p1); }
      for (i = 0; i < n-1; i++) { for (k = 0; k < n; k++) { B[i*n+k] = 0; } }
      cblas_dger(CblasRowMajor, n-1, n, 1.0, C, 1, B+(n-1)*n, 1, B, n);
      Dmod(p, B, n-1, n, n);
      for (i = n-2; i > j-2; i--)
	{
	  ++count;
	  cblas_dswap(n, B+i*n, 1, B+(i+1)*n, 1);
	}
      if (count % 2 == 0)
	for (i = 0; i < n*n; i++) { B[i] = fmod(det*B[i], p1); }
      else
	for (i = 0; i < n*n; i++) { B[i] = fmod((p-det)*B[i], p1); }	
      { XFREE(P); XFREE(rp); XFREE(C); }
      return B;
    }
}
  



/* 
 * Calling Sequence:
 *   r/-1 <-- mBasis(p, A, n, m, basis, nullsp, B, N)
 *
 * Summary:
 *   Compute a basis for the rowspace and/or a basis for the left nullspace 
 *   of a mod p matrix
 *  
 * Description:
 *   Given a n x m mod p matrix A, the function computes a basis for the
 *   rowspace B and/or a basis for the left nullspace N of A. Row vectors in 
 *   the r x m matrix B consist of basis of A, where r is the rank of A in
 *   Z/pZ. If r is zero, then B will be NULL. Row vectors in the n-r x n
 *   matrix N consist of the left nullspace of A. N will be NULL if A is full
 *   rank.
 *
 *   The pointers are passed into argument lists to store the computed basis 
 *   and nullspace. Upon completion, the rank r will be returned. The 
 *   parameters basis and nullsp control whether to compute basis and/or
 *   nullspace. If set basis and nullsp in the way that both basis and 
 *   nullspace will not be computed, an error message will be printed and 
 *   instead of rank r, -1 will be returned.
 *  
 * Input:
 *        p: FiniteField, prime modulus
 *           if p is a composite number, the routine will still work if no 
 *           error message is returned
 *        A: 1-dim Double array length n*m, representation of a n x m mod p 
 *           matrix. The entries of A are casted from integers
 *        n: long, row dimension of A
 *        m: long, column dimension of A
 *    basis: 1/0, flag to indicate whether to compute basis for rowspace or 
 *           not
 *         - basis = 1, compute the basis
 *         - basis = 0, not compute the basis
 *   nullsp: 1/0, flag to indicate whether to compute basis for left nullspace
 *           or not
 *         - nullsp = 1, compute the nullspace
 *         - nullsp = 0, not compute the nullspace
 *
 * Output:
 *   B: pointer to (Double *), if basis = 1, *B will be a 1-dim r*m Double
 *      array, representing the r x m basis matrix. If basis = 1 and r = 0, 
 *      *B = NULL
 *  
 *   N: pointer to (Double *), if nullsp = 1, *N will be a 1-dim (n-r)*n Double
 *      array, representing the n-r x n nullspace matrix. If nullsp = 1 and
 *      r = n, *N = NULL.
 *
 * Return:
 *   - if basis and/or nullsp are set to be 1, then return the rank r of A 
 *   - if both basis and nullsp are set to be 0, then return -1
 *
 * Precondition:
 *   n*(p-1)^2 <= 2^53-1 = 9007199254740991
 *
 * Note:
 *   - In case basis = 0, nullsp = 1, A will be destroyed inplace. Otherwise,
 *     A will not be changed.
 *   - Space of B and/or N will be allocated in the function
 *
 */

long
mBasis (const FiniteField p, Double *A, const long n, const long m, \
	const long basis, const long nullsp, Double **B, Double **N)
{
  long i, r;
  long *P, *rp;
  FiniteField d[1];
  Double *A1, *U;

  P = XMALLOC(long, n+1);
  for (i = 0; i < n+1; i++) { P[i] = i; }
  rp = XCALLOC(long, n+1);
  d[0] = 1;
  if ((basis == 1) && (nullsp == 1))
    {
      A1 = XMALLOC(Double, n*m);
      DCopy(n, m, A, m, A1, m);
      RowEchelonTransform(p, A1, n, m, 1, 1, 1, 0, P, rp, d);
      r = rp[0];
      U = XCALLOC(Double, n*n);
      for (i = 0; i < n; i++) { U[i*n+i] = 1; }
      if (r != 0) { DCopy(n, r, A1, m, U, n); }
      for (i = r; i > 0; i--)
	{ if (P[i] != i) { cblas_dswap(n, U+i-1, n, U+P[i]-1, n); } }
      if (r == 0) { *B = NULL; }
      else
	{
	  *B = XMALLOC(Double, r*m);
	  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, r, m, n, \
		      1.0, U, n, A, m, 0.0, *B, m);
	  Dmod(p, *B, r, m, m);
	}
      if (r ==n) { *N = NULL; }
      else
	{
	  *N = XMALLOC(Double, (n-r)*n);
	  DCopy(n-r, n, U+r*n, n, *N, n);;
	}
      { XFREE(A1); XFREE(U); XFREE(P); XFREE(rp); }
      return r;
    }
  else if ((basis == 1) && (nullsp == 0))
    {
      A1 = XMALLOC(Double, n*m);
      DCopy(n, m, A, m, A1, m);
      RowEchelonTransform(p, A1, n, m, 1, 0, 1, 0, P, rp, d);
      r = rp[0];
      if (r == 0) 
	{
	  *B = NULL; 
	  { XFREE(A1); XFREE(P); XFREE(rp); }
	  return r;
	}
      U = XCALLOC(Double, r*n);
      DCopy(r, r, A1, m, U, n);
      for (i = r; i > 0; i--)
	{ if (P[i] != i) { cblas_dswap(r, U+i-1, n, U+P[i]-1, n); } }
      *B = XMALLOC(Double, r*m);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, r, m, n, \
		  1.0, U, n, A, m, 0.0, *B, m);
      Dmod(p, *B, r, m, m);
      { XFREE(A1); XFREE(U); XFREE(P); XFREE(rp); }
      return r;
    }
  else if ((basis == 0) && (nullsp == 1))
    {
      RowEchelonTransform(p, A, n, m, 0, 1, 1, 0, P, rp, d);
      r = rp[0];
      if (r == n)
	{
	  *N = NULL;
	  { XFREE(P); XFREE(rp); }
	  return r; 
	}
      *N = XCALLOC(Double, (n-r)*n);
      if (r != 0) { DCopy(n-r, r, A+r*m, m, *N, n); }
      for (i = 0; i < n-r; i++) { (*N)[i*n+r+i] = 1; }
      for (i = r; i > 0; i--)
	{ if (P[i] != i) { cblas_dswap(n-r, *N+i-1, n, *N+P[i]-1, n); } }
      { XFREE(P); XFREE(rp); }
      return r;
    }
  else
    {
      iml_error("In mBasis, both basis and nullsp are zero.");
      return -1;
    }
}



/*
 * Calling Sequence:
 *   det <-- mDeterminant(p, A, n)
 * 
 * Summary:
 *   Compute the determinant of a square mod p matrix
 *
 * Input:
 *   p: FiniteField, prime modulus
 *      if p is a composite number, the routine will still work if no error 
 *      message is returned
 *   A: 1-dim Double array length n*n, representation of a n x n mod p matrix.
 *      The entries of A are casted from integers
 *   n: long, dimension of A
 *
 * Output:
 *   det(A) mod p, the determinant of square matrix A
 *
 * Precondition:
 *   ceil(n/2)*(p-1)^2+(p-1) <= 2^53-1 = 9007199254740991 (n >= 2)
 *
 * Note:
 *   A is destroyed inplace
 *
 */
long 
mDeterminant (const FiniteField p, Double *A, const long n)
{
  long i, count=0;
  long *P, *rp;
  FiniteField det, d[1];

  P = XMALLOC(long, n+1);
  for (i = 0; i < n+1; i++) { P[i] = i; }
  rp = XCALLOC(long, n+1);
  d[0] = 1;
  RowEchelonTransform(p, A, n, n, 0, 0, 0, 1, P, rp, d);
  det = d[0];
  if (det != 0)
    {
      for (i = 1; i < n+1; i++)	{ if (P[i] != i) { ++count; } }
      if (count % 2 == 0) 
	{ 
	  { XFREE(P); XFREE(rp); }
	  return det;
	}
      else
	{
	  { XFREE(P); XFREE(rp); }
	  return p-det;
	}
    }
  { XFREE(P); XFREE(rp); }     
  return det;
}


/*
 * Calling Sequence:
 *   1/0 <-- mInverse(p, A, n)
 * 
 * Summary:
 *   Certified compute the inverse of a mod p matrix inplace
 *
 * Description:
 *   Given a n x n mod p matrix A, the function computes A^(-1) mod p 
 *   inplace in case A is a nonsingular matrix in Z/Zp. If the inverse does
 *   not exist, the function returns 0.
 *
 *   A will be destroyed at the end in both cases. If the inverse exists, A is
 *   inplaced by its inverse. Otherwise, the inplaced A is not the inverse.
 *
 * Input:
 *   p: FiniteField, prime modulus
 *      if p is a composite number, the routine will still work if no error 
 *      message is returned
 *   A: 1-dim Double array length n*n, representation of a n x n mod p matrix.
 *      The entries of A are casted from integers
 *   n: long, dimension of A
 *
 * Return: 
 *   - 1, if A^(-1) mod p exists
 *   - 0, if A^(-1) mod p does not exist
 *
 * Precondition:
 *   ceil(n/2)*(p-1)^2+(p-1) <= 2^53-1 = 9007199254740991 (n >= 2)
 *
 * Note:
 *   A is destroyed inplace
 *
 */

long
mInverse (const FiniteField p, Double *A, const long n)
{
  long i;
  long *P, *rp;
  FiniteField d[1];

  P = XMALLOC(long, n+1);
  for (i = 0; i < n+1; i++) { P[i] = i; }
  rp = XCALLOC(long, n+1);
  d[0] = 1;
  RowEchelonTransform(p, A, n, n, 1, 1, 1, 0, P, rp, d);
  if (rp[0] == n) 
    {
      for (i = n; i > 0; i--)
	if (P[i] != i)
	  cblas_dswap(n, A+i-1, n, A+P[i]-1, n);
      { XFREE(P); XFREE(rp); }
      return 1;
    }
  { XFREE(P);  XFREE(rp); }
  return 0;
}




/*
 * Calling Sequence:
 *   r <-- mRank(p, A, n, m)
 *
 * Summary:
 *   Compute the rank of a mod p matrix
 *
 * Input:
 *   p: FiniteField, prime modulus
 *      if p is a composite number, the routine will still work if no 
 *      error message is returned
 *   A: 1-dim Double array length n*m, representation of a n x m mod p 
 *      matrix. The entries of A are casted from integers
 *   n: long, row dimension of A
 *   m: long, column dimension of A
 *   
 * Return:
 *   r: long, rank of matrix A
 *
 * Precondition:
 *   ceil(n/2)*(p-1)^2+(p-1) <= 2^53-1 = 9007199254740991 (n >= 2)
 *
 * Note:
 *   A is destroyed inplace
 *
 */

long 
mRank (const FiniteField p, Double *A, const long n, const long m)
{
  long i, r;
  long *P, *rp;
  FiniteField d[1];

  P = XMALLOC(long, n+1);
  for (i = 0; i < n+1; i++) { P[i] = i; }
  rp = XCALLOC(long, n+1);
  d[0] = 1;
  RowEchelonTransform(p, A, n, m, 0, 0, 0, 0, P, rp, d);
  r = rp[0];
  { XFREE(P); XFREE(rp); }
  return r;
}




/*
 * Calling Sequence:
 *   rp <-- mRankProfile(p, A, n, m)
 *
 * Summary:
 *   Compute the rank profile of a mod p matrix
 *
 * Input:
 *   p: FiniteField, prime modulus
 *      if p is a composite number, the routine will still work if no 
 *      error message is returned
 *   A: 1-dim Double array length n*m, representation of a n x m mod p 
 *      matrix. The entries of A are casted from integers
 *   n: long, row dimension of A
 *   m: long, column dimension of A
 *   
 * Return:
 *   rp: 1-dim long array length n+1, where
 *     - rp[0] is the rank of matrix A
 *     - rp[1..r] is the rank profile of matrix A
 *
 * Precondition:
 *   ceil(n/2)*(p-1)^2+(p-1) <= 2^53-1 = 9007199254740991 (n >= 2)
 *
 * Note:
 *   A is destroyed inplace
 *
 */

long *
mRankProfile (const FiniteField p, Double *A, const long n, const long m)
{
  long i;
  long *P, *rp;
  FiniteField d[1];

  P = XMALLOC(long, n+1);
  for (i = 0; i < n+1; i++) { P[i] = i; }
  rp = XCALLOC(long, n+1);
  d[0] = 1;
  RowEchelonTransform(p, A, n, m, 0, 0, 0, 0, P, rp, d);
  XFREE(P);
  return rp;
}


