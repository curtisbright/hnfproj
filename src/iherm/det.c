#include "iherm.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "arith_utils.h"
#include "imlsolve.h"
#include "highorder.h"
#include "reconstruct.h"
#include "dbg_print.h"
#include "timer.h"

static void mul_diag(mpz_t det, mpzMatrix_t const * A)
{
  long i;
  long n = A->nrows;
  mpz_set_ui(det, 1);
  for (i = 0; i < n; ++i) {
    mpz_mul(det, det, A->data[i*n + i]);
  }
}

static void mpzMatrix_realloc(mpzMatrix_t * M, long nrows, long ncols) /*TODO: don't realloc if (m,n) == (m,n) */
{
  long i;
  long nelms = nrows * ncols;
  for (i = 0; i < mpzMatrix_numElms(M); ++i) {
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

static void proj(mpzMatrix_t * X, mpz_t d, mpzMatrix_t const * B, rnsMatrix_t const * R, mpzMatrix_t const * v, rnsMatrix_t const * Binv)
{
  mpzMatrix_t * Rv = NULL;

  mpzMatrix_realloc(X, v->nrows, v->ncols);

  if (!R) {
    imlSolveInv(X, d, B, v, Binv);
  } else if (v->ncols != v->nrows) {
    Rv = mpzMatrix_init(v->nrows, v->ncols);
    mpzMatrix_rnsMpzMM(Rv, R, v);
    imlSolveInv(X, d, B, Rv, Binv);
    mpzMatrix_fini(Rv);
  } else {
    Rv = mpzMatrix_init(v->nrows, v->ncols);
    mpzMatrix_reconstruct(Rv, R);
    imlSolveInv(X, d, B, Rv, Binv);
    mpzMatrix_fini(Rv);
  }

}

static void chooseRHS(mpzMatrix_t * v, long n, long k)
{
  mpzMatrix_realloc(v, n, k);
  if (k == n) {
    mpzMatrix_identity(v);
  } else {
    mpzMatrix_rand(v, 10);
  }
}

static void hor(rnsMatrix_t ** R, rnsMatrix_t ** Binv, mpzMatrix_t const * B)
{
  if (*Binv) { basis_fini((*Binv)->basis); rnsMatrix_fini(*Binv); }
  if (*R) { basis_fini((*R)->basis); rnsMatrix_fini(*R); }
  *R = highOrderResidue_rns(B, Binv);
}

struct iherm_info {
  mpzMatrix_t const * A;
  mpzMatrix_t * B;

  pk_t * H;

  rnsMatrix_t * R;
  rnsMatrix_t * Binv;
};
typedef struct iherm_info iherm_info_t;

static int iter(iherm_info_t * s, int k)
{
  int ret = 0;
  long n = s->A->nrows;
  mpzMatrix_t * X = mpzMatrix_init(0, 0);
  mpzMatrix_t * v = mpzMatrix_init(0, 0);
  pk_t * Htmp = pkMatrix_init(n);
  pk_t * T = pkMatrix_init(n);
  pk_t tmp;

  mpz_t denom;
  mpz_init(denom);

  chooseRHS(v, n, k);

  TIMER("iml", 
  proj(X, denom, s->B, s->R, v, s->Binv);)

  TIMER("hermiteOfProj_balanced",
  hermiteOfProjection_balanced(T, X, denom);)

  /*TIMER("hermiteOfProj",
  hermiteOfProjection(T, X, denom);)*/

  TIMER("prod of hermite",
  pkMatrix_gemm_iter(Htmp, T, s->H);) /* dont use rnsgemm */
  tmp = *Htmp; *Htmp = *(s->H); *(s->H) = tmp; /*TODO: swap */

  TIMER("hnf of prod",
  pkMatrix_hnf(s->H);)

  if (k == n) { ret = 1; goto exit; }

  TIMER("applyHinv",
  pkMatrix_applyHinv(s->B, s->A, s->H);)

  TIMER("hor", 
  hor(&(s->R), &(s->Binv), s->B);)

  if (rnsMatrix_isZero(s->R)) { ret = 1; goto exit; }

  exit:
    mpz_clear(denom);
    mpzMatrix_fini(X);
    mpzMatrix_fini(v);
    pkMatrix_fini(Htmp);
    pkMatrix_fini(T);

  return ret;
}

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

static void iherm_fini(iherm_info_t * s)
{
  if(s->R) basis_fini(s->R->basis);
  if(s->Binv) basis_fini(s->Binv->basis); /*TODO: why (n=8)? */

  pkMatrix_fini(s->H);
  mpzMatrix_fini(s->B);
  rnsMatrix_fini(s->Binv);
  rnsMatrix_fini(s->R);
  free(s);
}

static void print_diagonal(FILE * stream, pk_t const * H)
{
  int j;
  for (j = 0; j < H->k; ++j) {
    mpz_out_str(stream, 10, pkgetc(H, H->idxs[j], j));
    fprintf(stream, ", ");
  }
}

int timer_depth;
mpzMatrix_t * myDet(mpzMatrix_t const * A, long k1, long k2)
{
  int i = 0;
  int done = 0;
  long kk [5] = {0};
  char buf[100];
  mpzMatrix_t * Hfull;
  iherm_info_t * s = iherm_init(A);

  kk[0] = k1 ? k1 : 8;
  kk[1] = k2 ? k2 : A->nrows/10;
  kk[2] = A->nrows;

  for(i = 0; !done; ++i) {
    dprintf(2, "\nITERATION %d (k = %ld; %ld/%ld)", i,  kk[i], s->H->k, A->ncols);
    sprintf(buf, "Iteration %d total", i);
    TIMER(buf,
    done = iter(s, kk[i]);)
  }
  dprintf(3, "H: \t");
  dprint2(3, print_diagonal, s->H);
  dprintf(3, "k: %ld/%ld", s->H->k, s->H->n);


  Hfull = pkMatrix_toFull(s->H);

  iherm_fini(s);

  return Hfull;
}

