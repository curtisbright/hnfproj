#include "unicert_utils.h"

#include <assert.h>
#include <stdlib.h>

#include "gmp.h"

#include "arith_utils.h"
#include "imlsolve.h"
#include "mpz_matrix.h"
#include "residue.h"

/* TODO: share this? */
static void shuffleRows(mpzMatrix_t * A, int trans)
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

mpzMatrix_t * unicert_input(long n, long l, int u)
{
  mpzMatrix_t * A = mpzMatrix_init(n, n);
  if (u) {
    mpzMatrix_identity(A);
    shuffleRows(A, 0);
    shuffleRows(A, 1);
  } else {
    mpzMatrix_rand(A, l);
  }
  return A;
}

int unicert_check(mpzMatrix_t const * A, int uni)
{
  long i, d;
  int rslt = 1;
  long p = (1L << 20);
  long n = A->nrows;
  residue_t * Ap;

  for (i = 0; i < 10; ++i) {
    p = nextprime(p);
    Ap = atlas_init(n, n, p);
    atlas_fromMpzMatrix(Ap, A);

    d = atlas_determinant(Ap);

    rslt &= ( (d == 1) || (d == (p-1)) );
  }

  return (rslt == uni);
}

int hor_check(mpzMatrix_t const * A, mpzMatrix_t const * R)
{
  long cR, i, n = A->nrows;
  mpzMatrix_t * x = mpzMatrix_init(n, n);
  mpz_t d;
  mpz_init(d);

  imlSolve(x, d, A, R);

  cR = 0;
  for(i = 0; i < n*n; ++i) {
    if (mpz_cmpabs(x->data[i], d) > 0) { cR += 1; }
  }

  mpz_clear(d);
  mpzMatrix_fini(x);

  return (cR == 0);
}
