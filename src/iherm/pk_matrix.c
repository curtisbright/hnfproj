#include "pk_matrix.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "thresholds.h"

void pkMatrix_print(FILE * stream, pk_t const * T)
{
  long i;
  mpzMatrix_t A;
  A.nrows = T->n;
  A.ncols = T->k;
  A.data = T->data;
  if(A.data) {
    mpzMatrix_print(stream, &A);
  }
  fprintf(stream, "  -----\n");
  fprintf(stream, "  ");
  for (i = 0; i< T->k; ++i) {
    fprintf(stream, "%ld\t", T->idxs[i]);
  }
  fprintf(stream, "\n");
}

/*static void pkMatrix_max(mpz_t m, pk_t const * T)
{
  mpzMatrix_t A;
  A.nrows = T->n;
  A.ncols = T->k;
  A.data = T->data;
  if(A.data) {
    mpzMatrix_max(m, &A);
  }
}*/

mpz_srcptr pkgetc(pk_t const * Z, int row, int col)
{
  long i;
  assert(0 <= row);
  assert(0 <= col);
  assert(row < Z->n);
  assert(col < Z->k);
  /*assert(row <= Z->idxs[col]);*/
  i = Z->k * row + col;
  return Z->data[i];
}
mpz_ptr pkget(pk_t * Z, int row, int col)
{
  return (mpz_ptr)pkgetc(Z, row, col);
}

static long * getNonTrivialIdxs(mpzMatrix_t const * A, long * k_ret)
{
  long i, c;
  long * idxs;
  long k = 0;
  for (i = 0; i < A->nrows; ++i) {
    if (mpz_cmp_ui(mmgetc(A, i, i), 1) != 0) {
      k += 1;
    }
  }
  idxs = malloc(k * sizeof(long));

  c = 0;
  for (i = 0; i < A->nrows; ++i) {
    if (mpz_cmp_ui(mmgetc(A, i, i), 1) != 0) {
      idxs[c] = i;
      c += 1;
    }
  }

  *k_ret = k;
  return idxs;
}

static long * makeRevIdxs(long n, long k, long * idxs)
{
  long i;
  long * rev_idxs = malloc(n * sizeof(long));
  memset(rev_idxs, -1, n*sizeof(long));
  for (i = 0; i < k; ++i) {
    rev_idxs[idxs[i]] = i;
  }
  return rev_idxs;
}

static long * makeMinIdxs(long n, long k, long * idxs)
{
  long i;
  long j = 0;
  long * min_idxs = malloc(n * sizeof(long));
  for (i = j = 0; i < n && j < k; ++i) {
    min_idxs[i] = j;
    if (i >= idxs[j]) { ++j; }
  }
  for ( ; i < n; ++i) {
    min_idxs[i] = k;
  }
  return min_idxs;
}

static void pkMatrix_initContents(pk_t * Z, long * idxs, long k, long n)
{
  long i;

  Z->n = n;
  Z->k = k;
  Z->idxs = idxs;
  Z->rev_idxs = makeRevIdxs(n, k, idxs);
  Z->min_idxs = makeMinIdxs(n, k, idxs);
  Z->data = malloc(Z->n*Z->k*sizeof(mpz_t));
  for (i = 0; i < Z->n*Z->k; ++i) {
    mpz_init(Z->data[i]);
  }
}

static void pkMatrix_finiContents(pk_t * Z)
{
  long i;

  free(Z->idxs);
  free(Z->rev_idxs);
  free(Z->min_idxs);
  for (i = 0; i < Z->n*Z->k; ++i) {
    mpz_clear(Z->data[i]);
  }
  free(Z->data);
  memset(Z, 0, sizeof(pk_t));
}

pk_t * pkMatrix_init(long n)
{
  return pkMatrix_initFromIdxs(NULL, 0, n);
}

pk_t * pkMatrix_initFromIdxs(long * idxs, long k, long n)
{
  pk_t * Z = malloc(sizeof(pk_t));
  pkMatrix_initContents(Z, idxs, k, n);
  return Z;
}

void pkMatrix_reinit(pk_t * Z, long n)
{
  pkMatrix_realloc(Z, NULL, 0, n);
}

#if 0
void pkMatrix_realloc(char const * str, pk_t * Z, long * idxs, long k, long n)
{
  pkMatrix_finiContents(Z);
  pkMatrix_initContents(Z, idxs, k, n);
}
#else
void pkMatrix_realloc(pk_t * Z, long * idxs, long k, long n)
{
  long i;

  assert(Z->n == n);
  if(k <= Z->k) {
    for (i = Z->n * k; i < Z->n*Z->k; ++i) { mpz_clear(Z->data[i]); } 
    Z->data = realloc(Z->data, Z->n*k * sizeof(mpz_t));
    for (i = 0; i < Z->n*k; ++i) { mpz_set_ui(Z->data[i], 0); }
    Z->k = k;
    free(Z->idxs);
    free(Z->rev_idxs);
    free(Z->min_idxs);
    Z->idxs = idxs;
    Z->rev_idxs = makeRevIdxs(n, k, idxs);
    Z->min_idxs = makeMinIdxs(n, k, idxs);
  } else {
    pkMatrix_finiContents(Z);
    pkMatrix_initContents(Z, idxs, k, n);
  }
}
#endif

void pkMatrix_fini(pk_t * Z)
{
  pkMatrix_finiContents(Z);
  free(Z);
}

int pkMatrix_equal(pk_t const * A, pk_t const * B)
{
  long i;

  if(A->n != B->n) { return 0; }
  if(A->k != B->k) { return 0; }

  for (i = 0; i < A->k; ++i) {
    if(A->idxs[i] != B->idxs[i]) { return 0; }
  }

  for (i = 0; i < A->n*A->k; ++i) {
    if (0 != mpz_cmp(A->data[i], B->data[i])) { return 0; }
  }

  return 1;
}

pk_t * pkMatrix_fromFull(mpzMatrix_t const * A)
{
  long k, row, acol, zcol;
  long * idxs = getNonTrivialIdxs(A, &k);
  pk_t * Z = pkMatrix_initFromIdxs(idxs, k, A->nrows);

  for (zcol = 0; zcol < Z->k; ++zcol) {
    acol = Z->idxs[zcol];
    for (row = 0; row <= Z->idxs[zcol]; ++row) {
      mpz_init_set(pkget(Z, row, zcol), mmgetc(A, row, acol));
    }
  }
  return Z;
}

mpzMatrix_t * pkMatrix_toFull(pk_t const * Z)
{
  long row, zcol, acol;
  mpzMatrix_t * A = mpzMatrix_init(Z->n, Z->n);
  mpzMatrix_identity(A);
  for (zcol = 0; zcol < Z->k; ++zcol) {
    acol = Z->idxs[zcol];
    for (row = 0; row < Z->n; ++row) {
      mpz_set(mmget(A, row, acol), pkgetc(Z, row, zcol));
    }
  }
  return A;
}

static long * mergeIdxs(long * S, long nS, long * T, long nT, long * k_ret)
{
  long i, j, k, val;

  long * idxs = malloc((nS + nT)*sizeof(long));
  i = 0;
  j = 0;
  k = 0;
  val = 0;
  while(i < nS || j < nT) {
    if (i < nS && j < nT) {
      if (S[i] <= T[j]) {
        val = S[i]; ++i;
      } else {
        val = T[j]; ++j;
      }
    } else if (i < nS) {
      val = S[i]; ++i;
    } else if (j < nT) {
      val = T[j]; ++j;
    } else {
      assert(0);
    }

    if( k == 0 || val != idxs[k-1] ) {
      idxs[k] = val;
      ++k;
    }
  }
  *k_ret = k;
  return idxs;
}

static int rev_idx(pk_t const * Z, long q)
{
  long ret = Z->rev_idxs[q];
  assert(ret != -1);
  assert(Z->idxs[ret] == q);
  return ret;
}

static int in_idx(pk_t const * A, long q)
{
  return A->rev_idxs[q] != -1;
}

void pkMatrix_hnf(pk_t * Z)
{
  /*        |--j-|
   *       col
   *        v
   *  [ H * * * *]
   *  [   H * * *]  <- row
   *  [     * * *]
   *  [     H * *] <- src_row
   *  [       * *]
   *  [       H *]
   *  [         H]
   */

  long row, col, j, col_start, row_start, src_row;
  mpz_t quo;

  if (Z->k == 0) { return; }

  mpz_init(quo);
  col_start = Z->k-1;
  row_start = Z->idxs[Z->k-1]-1;

  for (row = row_start; row >= 0; --row) {
    for (col = col_start; col < Z->k; ++col) {
      src_row = Z->idxs[col];
      mpz_fdiv_q(quo, pkget(Z, row, col), pkget(Z, src_row, col));

      for (j = col; j < Z->k; ++j) {
        mpz_submul( pkget(Z, row, j), quo, pkget(Z, src_row, j));
      }
    }

    if (col_start != 0 && Z->idxs[col_start-1] == row) {
      --col_start;
    }
  }

  mpz_clear(quo);
}

void pkMatrix_gemm_block(pk_t * Z, pk_t const * S, pk_t const * T)
{
  long j, k;
  long col, row, zcol;
  long * idxs = mergeIdxs(S->idxs, S->k, T->idxs, T->k, &k);


  mpzMatrix_t * C = mpzMatrix_init(S->n, T->k);
  mpzMatrix_t * B = mpzMatrix_init(S->k, T->k);
  mpzMatrix_t A;
  A.nrows = S->n;
  A.ncols = S->k;
  A.data = S->data;

  pkMatrix_realloc(Z, idxs, k, S->n);

  for (row = 0; row < S->k; ++row) {
    /*for (j = 0; j < T->k; ++j) {*/
    for (j = T->min_idxs[S->idxs[row]]; j < T->k; ++j) {
      mpz_set(mmget(B, row, j), pkgetc(T, S->idxs[row], j));
    }
  }

  /*mpzMatrix_gemm(C, &A, B);*/
  mpzMatrix_rnsGemm(C, &A, B);

  for (col = 0; col < C->ncols; ++col) {
    zcol = rev_idx(Z, T->idxs[col]);
    for(row = 0; row <= T->idxs[col]; ++row) {
      mpz_set( pkget(Z, row, zcol), mmgetc(C, row, col));
    }
  }

  for (zcol = 0; zcol < Z->k; ++zcol) {
    col = Z->idxs[zcol];

    if (in_idx(T, col)) {
      for (row = 0; row <= col; ++row) {
        if (in_idx(S, row)) { continue; }
          mpz_add(pkget(Z, row, zcol),
                  pkget(Z, row, zcol),
                  pkgetc(T, row, rev_idx(T, col)));
      }
    } else if (in_idx(S, col)) {
      for (row = 0; row <= col; ++row) {
          mpz_set( pkget(Z, row, zcol),
                   pkgetc(S, row, rev_idx(S, col)));
      }
    }
  }
}

void pkMatrix_gemm_iter(pk_t * Z, pk_t const * S, pk_t const * T)
{
  long row, col, zcol, scol, tcol, i, k;
  long * idxs;
  assert(Z != S && Z != T); /* no aliased arguments */

  if(S->k == 0 && T->k == 0) { 
    pkMatrix_realloc(Z, 0, 0, S->n);
    return;
  }

  idxs = mergeIdxs(S->idxs, S->k, T->idxs, T->k, &k);
  pkMatrix_realloc(Z, idxs, k, S->n);

  for (row = Z->n-1; row >= 0; --row) {
    for (zcol = 0; zcol < Z->k; ++zcol) {
      col = Z->idxs[zcol];
      if (col < row) { continue; }
      if (in_idx(T, col)) {
        tcol = rev_idx(T, col);
        for (i = 0; i < S->k; ++i) {
          mpz_addmul( pkget(Z, row, zcol),
                      pkgetc(S, row, i),
                      pkgetc(T, S->idxs[i], tcol));
        }
        if (!in_idx(S,row)) {
          mpz_add( pkget(Z, row, zcol), pkget(Z, row, zcol), pkgetc(T, row, tcol) );
        }
      } else if (in_idx(S, col)) {
        scol = rev_idx(S, col);
        mpz_set(pkget(Z, row, zcol), pkgetc(S, row, scol));
      }
    }
  }
}

void pkMatrix_gemm(pk_t * Z, pk_t const * S, pk_t const * T)
{
  if(S->k >= PK_GEMM_BLOCK_THRESHOLD && T->k >= PK_GEMM_BLOCK_THRESHOLD) {
    pkMatrix_gemm_block(Z, S, T);
  } else {
    pkMatrix_gemm_iter(Z, S, T);
  }
}

void pkMatrix_applyVector_simple(mpzMatrix_t * X, pk_t const * Z, long start, long count)
{
  int row, col, j, xrow;

  for (j = 0; j < count; ++j) {
    for (col = 0; col < Z->k; ++col) {
      xrow = Z->idxs[col];
      for (row = 0; row < xrow; ++row) {
        mpz_addmul(mmget(X, row, start+j), mmget(X, xrow, start+j), pkgetc(Z, row, col));
      }
      mpz_mul(mmget(X, xrow, start+j), mmget(X, xrow, start+j), pkgetc(Z, xrow, col));
    }
  }
}

void pkMatrix_applyVector_block(mpzMatrix_t * X, pk_t const * Z, long start, long count)
{
  long row, col, j;

  mpzMatrix_t * xx = mpzMatrix_init(Z->k, count);
  mpzMatrix_t * Y = mpzMatrix_init(Z->n, count);
  mpzMatrix_t A;

  A.nrows = Z->n;
  A.ncols = Z->k;
  A.data = Z->data;

  for (row = 0; row < Z->k; ++row) {
    for (j = 0; j < count; ++j) {
      mpz_set(mmget(xx, row, j), mmgetc(X, Z->idxs[row], j+start));
    }
  }
   mpzMatrix_rnsGemm(Y, &A, xx);

  for (row = 0; row < Y->nrows; ++row) {
    if (in_idx(Z, row)) { continue; }
    for (j = 0; j < count; ++j) {
      mpz_add(mmget(Y, row, j), mmget(Y, row, j), mmgetc(X, row, j+start));
    }
  }

  for (row = 0; row < Y->nrows; ++row) {
    for (col = 0; col < count; ++col) {
      mpz_set(mmget(X, row, col+start), mmgetc(Y, row, col));
    }
  }

  mpzMatrix_fini(xx);
  mpzMatrix_fini(Y);
}

void pkMatrix_applyHinv(mpzMatrix_t * B, mpzMatrix_t const * A, pk_t const * H)
{
  int hcol, row, i, bcol;
  mpz_t t;

  assert(B->nrows == A->ncols && B->nrows == A->ncols);
  mpzMatrix_set(B, A);
  mpz_init(t);
  for (hcol = 0; hcol < H->k; ++hcol) {
    for (row = 0; row < H->n; ++row) {
      mpz_set_ui(t, 0);
      for (i = 0; i < H->idxs[hcol]; ++i) {
        mpz_addmul(t, mmget(B, row, i), pkgetc(H, i, hcol));
      }
      bcol = H->idxs[hcol];
      mpz_sub(mmget(B, row, bcol), mmget(B, row, bcol), t);
      mpz_tdiv_q( mmget(B, row, bcol),
                  mmget(B, row, bcol),
                  pkgetc(H, bcol, hcol) );

    }
  }
  mpz_clear(t);
}
