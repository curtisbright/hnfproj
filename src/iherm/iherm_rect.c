
#include <assert.h>
#include <stdlib.h>

#include "arith_utils.h"
#include "iml.h"
#include "imlsolve.h"
#include "iherm.h"
#include "residue.h"
#include "mpz_matrix.h"
#include "pk_matrix.h"
#include "timer.h"

/*****************************************************************************/

#define MAX(a,b) ( (a) >= (b) ? (a) : (b) )
static long * rankProfile(mpzMatrix_t const * A, long * rank)
{
  long * rp, *rpt;
  long i;
  long p = prevprime(pickStartModulus(MAX(A->nrows, A->ncols)));
  residue_t * Ap = residue_init(A->nrows, A->ncols, p);

  residue_fromMpzMatrix(Ap, A);
  rpt = mRankProfile(p, Ap->_data, Ap->nrows, Ap->ncols);
  *rank = rpt[0];
  rp = malloc(sizeof(long) * (*rank));
  for (i = 0; i < *rank; ++i) { rp[i] = rpt[i+1]-1; }

  free(rpt);
  residue_fini(Ap);

  return rp;
}

static long * columnPermutation(long rank, long const * rp, long ncols)
{
  long * Q = malloc(sizeof(long)*ncols);
  long i,j, t;

  for (i = 0; i < ncols; ++i) { Q[i] = i; }

  for (i = 0; i < rank; ++i) {
    j = rp[i];
    if (i == j) { continue; }
    t = Q[j]; Q[j] = Q[i]; Q[i] = t;
  }

  return Q;
}

static void splitColumns(mpzMatrix_t const * A, long rank, long * Q, mpzMatrix_t ** pAr, mpzMatrix_t ** pAn)
{
  long i, col;
  mpzMatrix_t * Ar = mpzMatrix_init(A->nrows, rank);
  mpzMatrix_t * An = mpzMatrix_init(A->nrows, A->ncols-rank);

  for (col = 0; col < rank; ++col) {
    for (i = 0; i < Ar->nrows; ++i) {
      mpz_set( mmget(Ar, i, col), mmgetc(A,i,Q[col])); }
  }

  for (col = 0; col < A->ncols-rank; ++col) {
    for (i = 0; i < An->nrows; ++i) {
      mpz_set( mmget(An, i, col), mmgetc(A,i,Q[rank+col])); }
  }

  *pAr = Ar;
  *pAn = An;
}

#if 0
static mpzMatrix_t * extendToSquare_random(mpzMatrix_t const * Ar)
{
  /*      [    | * ]
   * B := [ Ar | * ]
   *      [    | * ]
   */
  long row, col;
  long r = Ar->ncols;
  long k = Ar->nrows - Ar->ncols;
  mpzMatrix_t * B = mpzMatrix_init(r+k, r+k);

  for (row = 0; row < r+k; ++row) {
    for (col = 0; col < r; ++col) {
      mpz_set( mmget(B, row, col), mmgetc(Ar, row, col) );
    }
  }

  /* fill last k columns randomly */
  for (row = 0; row < r+k; ++row) {
    for (col = r; col < r+k; ++col) {
      mpz_set_ui( mmget(B, row, col), rand() % 10);
    }
  }
  return B;
}
#endif

static long * nonRankCols(long const * rp, long rank, long n)
{
  long k = n - rank;
  long * nrp = malloc(k * sizeof(long));
  long i;
  long col = 0;
  long j = 0;
  for (i = 0; i < rank; ++i,++col) {
    while(col != rp[i]) {
      nrp[j] = col;
      ++j; ++col;
    }
  }
  for( ; j < k; ++j, ++col) {
    nrp[j] = col;
  }
  return nrp;
}

static mpzMatrix_t * extendToSquare(mpzMatrix_t const * Ar)
{
  long row, col, rank, i;
  long *rp, *nrp;
  long r = Ar->ncols;
  long k = Ar->nrows - Ar->ncols;
  mpzMatrix_t * B = mpzMatrix_init(r, r+k);
  mpzMatrix_t * C;

  /* B = column reversal of A^T */
  for (row = 0; row < r+k; ++row) {
    for (col = 0; col < r; ++col) {
      mpz_set( mmget(B, col, r+k-1-row), mmgetc(Ar, row, col) );
    }
  }

  TIMER("rp", rp = rankProfile(B, &rank);)
  mpzMatrix_fini(B);
  assert(rank == r);
  nrp = nonRankCols(rp, rank, r+k);

  C = mpzMatrix_init(r+k, r+k);
  for (row = 0; row < r+k; ++row) {
    for (col = 0; col < r; ++col) {
      mpz_set( mmget(C, row, col), mmgetc(Ar, row, col) );
    }
  }

  for(i = 0; i < k; ++i) {
    row = r+(k-nrp[i]-1);
    mpz_set_ui( mmget(C, row, r+i), 1);
  }

  free(rp);
  free(nrp);
  return C;
}

static mpzMatrix_t * nonRankHermiteCols(mpzMatrix_t const * U, mpzMatrix_t const * An)
{
  mpzMatrix_t * Hn = mpzMatrix_init(An->nrows, An->ncols);
  mpz_t denom; mpz_init(denom);
  if(An->ncols != 0) {
    imlSolve(Hn, denom, U, An);
    assert(0 == mpz_cmp_ui(denom, 1));
  }

  mpz_clear(denom);
  return Hn;
}

static mpzMatrix_t * assembleHermite(long rank, mpzMatrix_t * const Hb, mpzMatrix_t * const Hn, long * Q)
{
  long col, i;
  mpzMatrix_t * H = mpzMatrix_init(Hn->nrows, rank + Hn->ncols);
  for (col = 0; col < rank; ++col) {
    for (i = 0; i < Hb->nrows; ++i) {
      mpz_set(mmget(H, i, Q[col]), mmget(Hb, i, col));
    }
  }
  for (col = rank; col < H->ncols; ++col) {
    for (i = 0; i < Hn->nrows; ++i) {
      mpz_set(mmget(H, i, Q[col]), mmget(Hn, i, col-rank));
    }
  }
  return H;
}

static int checkShape(long rank, long const * rp, mpzMatrix_t const * H)
{
  /* [ * * * * * * * ]
   * [     * * * * * ]
   * [           * * ]
   * [               ]
   *   ^   ^     ^     rank profile cols
   */

  long nrows = H->nrows;
  long ncols = H->ncols;
  long i, j;

  for (i = rank; i < nrows; ++i) {
    for (j = 0; j < ncols; ++j) {
      if (0 != mpz_cmp_ui(mmgetc(H, i, j), 0)) {
        return 0;
      }
    }
  }

  for (i = 0; i < rank; ++i) {
    for (j = 0; j < ncols; ++j) {
      if (0 != mpz_cmp_ui(mmgetc(H, i, j), 0)) {
        break;
      }
    }
    if (rp[i] != j) { return 0;}
  }

  return 1;
}

static void compress(mpzMatrix_t * Ar, mpzMatrix_t * An)
{
  mpzMatrix_t * L, *La;
  const long K = 10;
  if(Ar->nrows - Ar->ncols <= K) { return; }

  L = mpzMatrix_init(Ar->ncols+K, Ar->nrows);
  mpzMatrix_rand(L, 1);

  /* Ar := L.Ar */
  La = mpzMatrix_init(Ar->ncols+K, Ar->ncols);
  mpzMatrix_rnsGemm(La, L, Ar);
  mpzMatrix_swap(La, Ar);
  mpzMatrix_fini(La);

  /* An := L.An */
  La = mpzMatrix_init(Ar->ncols+K, An->ncols);
  mpzMatrix_rnsGemm(La, L, An);
  mpzMatrix_swap(La, An);
  mpzMatrix_fini(La);

  mpzMatrix_fini(L);
}

static void preamble(mpzMatrix_t const * A, mpzMatrix_t ** pB, mpzMatrix_t ** pAn,
                     long * prank, long ** prp, long ** pQ)
{
  mpzMatrix_t * Ar, *An, *B;
  long *rp, *Q;
  long rank;

  TIMER("Rank profile",
  rp = rankProfile(A, &rank);)

  TIMER("permute",
  Q = columnPermutation(rank, rp, A->ncols);)

  TIMER("split",
  splitColumns(A, rank, Q, &Ar, &An);)

  TIMER("compress",
  compress(Ar, An);)

  TIMER("extend",
  B = extendToSquare(Ar);)

  mpzMatrix_fini(Ar);
  *pB = B;
  *pAn = An;
  *prank = rank;
  *prp = rp;
  *pQ = Q;
}

mpzMatrix_t * hermiteRect(mpzMatrix_t const * A)
{
  int rc;
  long * rp, * Q;
  long rank;
  mpzMatrix_t *An, *B, *U, *Hb, *Hn, *H;

  TIMER("Preamble",
  preamble(A, &B, &An, &rank, &rp, &Q);
  )

  TIMER("Nonsing hermite",
  Hb = hermiteWithTransform(B, &U);)

  TIMER("unisolve",
  Hn = uniSolve(U, An);)
  mpzMatrix_fini(Hn);

  TIMER("Nonrank cols",
  Hn = nonRankHermiteCols(U, An);)

  TIMER("assemble",
  H = assembleHermite(rank, Hb, Hn, Q);)

  TIMER("check shape",
  rc = checkShape(rank, rp, H);)
  printf("rc: %d\n", rc);

  mpzMatrix_fini(An);
  mpzMatrix_fini(B);
  mpzMatrix_fini(Hb);
  mpzMatrix_fini(Hn);
  mpzMatrix_fini(U);
  free(Q);
  free(rp);


  return H;
}



