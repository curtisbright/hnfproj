#include "iherm_utils.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "gmp.h"
#include "iml.h"

#include "arith_utils.h"
#include "residue.h"

static void addRow(mpzMatrix_t * A, long dst, long src, long k)
{
  int i;
  for (i = 0; i < A->ncols; ++i) {
    mpz_addmul_ui(mmget(A, dst, i), mmgetc(A, src, i), k);
  }
}
static void addCol(mpzMatrix_t * A, long dst, long src, long k)
{
  int i;
  for (i = 0; i < A->nrows; ++i) {
    mpz_addmul_ui(mmget(A, i, dst), mmgetc(A, i, src), k);
  }
}

static void swapRow(mpzMatrix_t * A, long i, long j)
{
  long k;
  for (k = 0; k < A->ncols; ++k) {
    mpz_swap(mmget(A, i, k), mmget(A, j, k));
  }
}
static void swapCol(mpzMatrix_t * A, long i, long j)
{
  long k;
  for (k = 0; k < A->nrows; ++k) {
    mpz_swap(mmget(A, k, i), mmget(A, k, j));
  }
}

static void shuffleRows(mpzMatrix_t * A, int trans)
{
  long i,j;
  long n = trans ? A->ncols : A->nrows;
  void (*swap)(mpzMatrix_t *, long, long) = trans ? swapCol : swapRow;
  for (i = n-1; i >= 0; --i) {
    j = rand() % (i+1); /* 0 <= j <= i */
    swap(A, i, j);
  }
}

static mpzMatrix_t * randHerm_rect(long nrows, long ncols, long rank, long bitlen)
{
  long i, j;
  mpzMatrix_t * A = mpzMatrix_init(nrows, ncols);
  mpzMatrix_t * R = mpzMatrix_init(rank, rank);
  mpzMatrix_rand(R, bitlen);
  assert(nrows && ncols && rank && bitlen);
  for (i = 0; i < rank; ++i) {
    for (j = 0; j < rank; ++j) {
      mpz_set(mmget(A, i, j), mmgetc(R, i, j));
    }
  }
  for (i = rank; i < nrows; ++i) {
    addRow(A, i, rand()%rank, (rand()%10)+1);
  }
  for (i = rank; i < ncols; ++i) {
    addCol(A, i, rand()%rank, (rand()%10)+1);
  }

  shuffleRows(A, 0);
  shuffleRows(A, 1);
  mpzMatrix_fini(R);
  return A;
}

/*****************************************************************************/


static void uniShuffleRows(mpzMatrix_t * A, int trans)
{
  long n = A->nrows;
  long r1 = 0;
  long r2 = 0;
  long i,j;
  void (*func)(mpz_ptr, mpz_srcptr, mpz_srcptr);

  for (i = 0; i < 10*n; ++i) {
    do {
      r1 = rand() % n;
      r2 = rand() % n;
    } while (r1 == r2);
    func = (rand() % 2) ? mpz_add : mpz_sub;
    for (j = 0; j < n; ++j) {
      if (trans) {
        func(mmget(A, j, r1), mmget(A, j, r1), mmget(A, j, r2));
      } else {
        func(mmget(A, r1, j), mmget(A, r1, j), mmget(A, r2, j));
      }
    }
  }
}

static mpzMatrix_t * diagcat(mpzMatrix_t const * S, mpzMatrix_t const * T)
{
  long i, j;
  long Srows = S->nrows;
  long Scols = S->ncols;

  long Trows = T->nrows;
  long Tcols = T->ncols;

  long Arows = Srows + Trows;
  long Acols = Scols + Tcols;

  mpzMatrix_t * A = mpzMatrix_init(Arows, Acols);
  for (i = 0; i < Srows; ++i) {
    for (j = 0; j < Scols; ++j) {
      mpz_set(mmget(A, i, j), mmgetc(S, i, j));
    }
  }
  for (i = 0; i < Trows; ++i) {
    for (j = 0; j < Tcols; ++j) {
      mpz_set(mmget(A, i+Srows, j+Scols), mmgetc(T, i, j));
    }
  }
  return A;
}


static mpzMatrix_t * nontriv_steel(long n, long l)
{
  mpzMatrix_t * A = mpzMatrix_init(n, n);
  long i;

  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, rand());
  for (i = 0; i < n; ++i) {
    mpz_urandomb(A->data[i*n + i], state, l);
    mpz_add_ui(A->data[i*n + i], A->data[i*n+i], 1);
  }
  gmp_randclear(state);
  uniShuffleRows(A, 0);
  uniShuffleRows(A, 1);

  return A;
}

static mpzMatrix_t * randHerm_double(long n, long bitlen)
{
  mpzMatrix_t * A;
  mpzMatrix_t * A1 = mpzMatrix_init(n/2, n/2);
  mpzMatrix_t * A2 = mpzMatrix_init(n/2, n/2);
  mpzMatrix_rand(A1, bitlen);
  mpzMatrix_rand(A2, bitlen);

  A = diagcat(A1, A2);

  uniShuffleRows(A, 0);

  mpzMatrix_fini(A1);
  mpzMatrix_fini(A2);

  return A;
}

static mpzMatrix_t * nontriv_jaeger(long n)
{
  long i, j;
  mpz_t m, s, t, x;
  mpzMatrix_t * A;
  mpz_inits(m, s, t, x, 0);

  n = nextprime(n);
  mpz_set_ui(m, n);

  A = mpzMatrix_init(n, n);
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      mpz_set_ui(s, i);
      mpz_set_ui(t, j);
      mpz_powm(x, s, t, m);
      mpz_set(mmget(A, i, j), x);
    }
  }
  mpz_clears(s, t, x, m, 0);

  return A;
}

static int check_det(mpzMatrix_t const * A, mpzMatrix_t * H)
{
  long i, d1, d2;
  int rslt = 1;
  long p = (1L << 10);
  long n = A->nrows;
  residue_t * Ap;
  mpz_t d;
  mpz_init_set_ui(d, 1);
  for (i = 0; i < n; ++i) {
    mpz_mul(d, d, H->data[i*n+i]);
  }

  for (i = 0; i < 10; ++i) {
    p = nextprime(p);
    Ap = atlas_init(n, n, p);
    atlas_fromMpzMatrix(Ap, A);

    d1 = atlas_determinant(Ap);

    d2 = mpz_fdiv_ui(d, p);

    rslt &= ( (d1 == d2) || (d1 == (p-d2)) );

    atlas_fini(Ap);
  }
  mpz_clear(d);

  return rslt;
}

static int check_full(mpzMatrix_t const * A, mpzMatrix_t * H)
{
  long i;
  long d;
  int rc;
  residue_t * Ap;
  residue_t * Hp;
  residue_t * U;
  long p = (1L << 20);
  long n = A->nrows;

  int rslt = 1;
  for (i = 0; i < 10; ++i) {
    p = nextprime(p);
    Ap = atlas_init(n, n, p);
    Hp = atlas_init(n, n, p);
    U = atlas_init(n, n, p);
    atlas_fromMpzMatrix(Ap, A);
    atlas_fromMpzMatrix(Hp, H);
    rc = atlas_inverse(Ap);
    if (!rc) { continue; }
    atlas_gemm(U, Ap, Hp);

    d = atlas_determinant(U);
    rslt &= ( (d == 1) || (d == p-1) );

    atlas_fini(Ap);
    atlas_fini(Hp);
    atlas_fini(U);
  }

  return rslt;
}

int iherm_check(mpzMatrix_t const * A, mpzMatrix_t * H)
{
  int rc = check_det(A, H);
  rc &= check_full(A, H);

  return rc;
}

mpzMatrix_t * iherm_input(char const * type, long nrows, long ncols, long rank, long bitlen)
{
  mpzMatrix_t * A = NULL;

  if (0 == strcmp(type, "random")) {
    /* random with l-bit entries */
    A = mpzMatrix_init(nrows, ncols);
    mpzMatrix_rand(A, bitlen);
  } else if (0 == strcmp(type, "steel")) {
    /* nontrivial diagonal with smooth entries (cf. Steel) */
    A = nontriv_steel(nrows, bitlen);
  } else if (0 == strcmp(type, "jaeger")) {
    /* a_ij = (i-1)^(j-1) mod n */
    A = nontriv_jaeger(nrows);
  } else if (0 == strcmp(type, "diagcat")) {
    /* [1 1 1 1 X] [1 1 1 1 X] */
    A = randHerm_double(nrows, bitlen);
  } else if (0 == strcmp(type, "rect")) {
    A = randHerm_rect(nrows, ncols, rank, bitlen);
  }

  return A;
}
