#include "libhnfproj.h"
#include "iherm.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gmp.h"

#include "arith_utils.h"
#include "dbg_print.h"
#include "imlsolve.h"
#include "highorder.h"
#include "hnfproj.h"
#include "reconstruct.h"
#include "timer.h"

static void mpzMatrix_realloc(mpzMatrix_t * M, long nrows, long ncols) /*TODO: don't realloc if (m,n) == (m,n) */
{
  long i;
  long nelms = nrows * ncols;
  for (i = 0; i < mpzMatrix_numElems(M); ++i) {
    mpz_clear(M->data[i]);
  }
  free(M->data);
  M->data = (mpz_t *)malloc(nelms*sizeof(mpz_t));
  for (i = 0; i < nelms; ++i) {
    mpz_init(M->data[i]);
  }
  M->nrows = nrows;
  M->ncols = ncols;
}

/** (X/d) := (1/B).R.v; uses precomputed (1/B) if available. */
static void projection(mpzMatrix_t * X, mpz_t d, mpzMatrix_t const * B, rnsMatrix_t const * R, mpzMatrix_t const * v, rnsMatrix_t const * Binv)
{
  mpzMatrix_t * Rv = NULL;

  mpzMatrix_realloc(X, v->nrows, v->ncols);

  if (!R) {
    /* (X/d) := (1/B).v */
    imlSolveInv(X, d, B, v, Binv);
  } else if (v->ncols != v->nrows) {
    /* (X/d) := (1/B).(R.v)
     * Keep R in rns representation. */
    Rv = mpzMatrix_init(v->nrows, v->ncols);
    mpzMatrix_gemmRnsMpz(Rv, R, v);
    imlSolveInv(X, d, B, Rv, Binv);
    mpzMatrix_fini(Rv);
  } else {
    /* (X/d) := (1/B).(R)
     * as v is identity matrix. */
    Rv = mpzMatrix_init(v->nrows, v->ncols);
    mpzMatrix_reconstruct(Rv, R);
    imlSolveInv(X, d, B, Rv, Binv);
    mpzMatrix_fini(Rv);
  }
}

/** Choose random vector for projection right hand size. */
static void chooseRHS(mpzMatrix_t * v, long nrows, long proj_size)
{
  mpzMatrix_realloc(v, nrows, proj_size);
  if (proj_size == nrows) {
    mpzMatrix_identity(v);
  } else {
    mpzMatrix_rand(v, 10);
  }
}

/** High order residue of matrix B.
  * Also returns (1/B). */
static void hor(rnsMatrix_t ** R, rnsMatrix_t ** Binv, mpzMatrix_t const * B)
{
  if (*Binv) { rnsMatrix_fini(*Binv); }
  if (*R) { rnsMatrix_fini(*R); }
  *R = highOrderResidue_rns(B, Binv);
}

/** HNF algorithm state. */
struct iherm_info {
  mpzMatrix_t const * A;
  mpzMatrix_t * B;
  pk_t * H;
  rnsMatrix_t * R;
  rnsMatrix_t * Binv;
};
typedef struct iherm_info iherm_info_t;

/** Initialize HNF algorithm; return new algorithm state. */
static iherm_info_t * iherm_init(mpzMatrix_t const * A)
{
  iherm_info_t * s = malloc(sizeof(iherm_info_t));
  s->A = A;
  s->B = mpzMatrix_init(A->nrows, A->ncols);
  s->H = pkMatrix_init(A->nrows);
  s->R = NULL;
  s->Binv = NULL;

  mpzMatrix_set(s->B, A);

  return s;
}

/** Finialize HNF algorithm state. */
static void iherm_fini(iherm_info_t * s)
{
  pkMatrix_fini(s->H);
  if (s->B) mpzMatrix_fini(s->B);
  rnsMatrix_fini(s->Binv);
  rnsMatrix_fini(s->R);
  free(s);
}

/** Perform one iteration of HNF algorithm.
  * Algorithm state s modified in place.
  * \return 1 if process has completed; 0 otherwise. */
static int iherm_iteration(iherm_info_t * s, int proj_size, int need_B)
{
  int ret = 0;
  long nrows = s->A->nrows;
  mpzMatrix_t * X = mpzMatrix_init(0, 0);
  mpzMatrix_t * v = mpzMatrix_init(0, 0);
  pk_t * Htmp = pkMatrix_init(nrows);
  pk_t * T = pkMatrix_init(nrows);
  pk_t tmp;
  int need_hor = (proj_size < nrows);
  int need_apply = (need_hor || need_B);

  mpz_t denom;
  mpz_init(denom);

  chooseRHS(v, nrows, proj_size);

  TIMER("iml", 
  projection(X, denom, s->B, s->R, v, s->Binv);)

  TIMER("hermiteOfProj_balanced",
  hermiteOfProjection_balanced(T, X, denom);)

  /*TIMER("hermiteOfProj",
  hermiteOfProjection_simple(T, X, denom);)*/

  TIMER("prod of hermite",
  pkMatrix_gemm_iter(Htmp, T, s->H);) /* dont use rnsgemm */
  tmp = *Htmp; *Htmp = *(s->H); *(s->H) = tmp; /*TODO: swap */

  TIMER("hnf of prod",
  pkMatrix_hnf(s->H);)

  if (need_apply) {
    TIMER("applyHinv", pkMatrix_applyHinv(s->B, s->A, s->H);)
  }
  if (need_hor) {
    TIMER("hor", hor(&(s->R), &(s->Binv), s->B);)
    /* Remaining matrix is unimodular; stop.*/
    if (rnsMatrix_isZero(s->R)) { ret = 1; }
  } else {
    /* Skip hor check: entirety of inverse must have been captured. */
    ret = 1;
  }

  mpz_clear(denom);
  mpzMatrix_fini(X);
  mpzMatrix_fini(v);
  pkMatrix_fini(Htmp);
  pkMatrix_fini(T);

  return ret;
}


/* External entry points for HNF.
 */
#define MAX(a,b) ( (a) >= (b) ? (a) : (b) )
static mpzMatrix_t * _hermite(mpzMatrix_t const * A, mpzMatrix_t ** U)
{
  int i = 0;
  int done = 0;
  long proj_size[5] = {0};
  mpzMatrix_t * Hfull;
  iherm_info_t * s = iherm_init(A);

  proj_size[0] = 8;
  proj_size[1] = MAX(A->nrows/10, 8);
  proj_size[2] = A->nrows;

  for(i = 0; !done; ++i) {
    dprintf(1, "\nITERATION %d (proj_size = %ld; k = %ld; n = %ld)", i,  proj_size[i], s->H->k, A->ncols);
    TIMERF("Iteration %d total",
    done = iherm_iteration(s, proj_size[i], !!U);, i)
  }

  Hfull = pkMatrix_toFull(s->H);
  if(U) { *U = s->B; s->B = 0; }

  iherm_fini(s);

  return Hfull;
}

mpzMatrix_t * hermite(mpzMatrix_t const * A)
{
  return _hermite(A, NULL);
}

mpz_t * hermiteLong(long const * A, long n)
{
  mpzMatrix_t * A_mpz = mpzMatrix_initSet(n, n, A);
  mpzMatrix_t * H = hermite(A_mpz);

  mpzMatrix_fini(A_mpz);

  return mpzMatrix_data(H);
}

mpz_t * hermiteMpz(mpz_t const * A, long n)
{
  mpzMatrix_t * H;
  mpzMatrix_t const * A_mpz;

  A_mpz = mpzMatrix_initFromMpz(n, n, A);
  H = hermite(A_mpz);
  free((void*)A_mpz);

  return mpzMatrix_data(H);
}

mpzMatrix_t * hermiteWithTransform(mpzMatrix_t const * A, mpzMatrix_t ** U)
{
  return _hermite(A, U);
}
