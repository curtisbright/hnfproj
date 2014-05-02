#include "rns_conversion.h"

#include <assert.h>
#include <pthread.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "cblas.h"

#include "arith_utils.h"
#include "basis.h"
#include "reconstruct.h"
#include "rns_matrix.h"

/* TODO: switch between recon conversion and complicated recon */

static double * correctionFactor(rnsMatrix_t const * Ap)
{
  long i;
  basis_t const * Px = Ap->basis;
  long sz = Ap->nrows * Ap->ncols;
  long nP = Ap->nmod;
  double * r = malloc(sz*sizeof(double));
  memset(r, 0, sz*sizeof(double));

  for (i = 0; i < nP; ++i) {
    /* r += Ap[i]*W[i] */
    cblas_daxpy(sz, Px->W[i], Ap->residues[i]->_data, 1, r, 1);
  }
  for (i = 0; i < sz; ++i) {
    r[i] = round(r[i]);
  }

  return r;
}

static void lagrangeEval(residue_t * r, rnsMatrix_t const * Ap)
{
  /*       nP
         -----
          \    [                        ]
   r :=    )   [ Ap[i] * L[i] * invL[i] ]
          /    [                        ] q
         -----
         i = 1
   */
  long i;
  basis_t const * Px = Ap->basis;
  long nP = Px->nmod;
  long k;

  residue_zero(r);
  for (i = 0; i < nP; ++i) {
    k = modsMpzL(Px->C[i], r->mod.p);
    residue_add(r, Ap->residues[i], k);
  }
}

static void rnsMatrix_convertOne(residue_t * Aq_i, rnsMatrix_t const * Ap, double const * corr)
{
  /*       /    nP                               \
           |  -----                              |
           |   \    [                        ]   |    [           ]
   Aq_i := |    )   [ Ap[i] * L[i] * invL[i] ]   |  - [ PP * corr ]
           |   /    [                        ] q |    [           ]q
           |  -----                              |
           \  i = 1                              /
  */

  long sz = Ap->nrows * Ap->ncols;
  /* k := -(PP mod q[i]) */
  double k = -modsMpzL(Ap->basis->PP, Aq_i->mod.p);
  lagrangeEval(Aq_i, Ap);
  /* Aq_i := Aq_i + k * corr mod q[i] */
  cblas_daxpy(sz, k, corr, 1, Aq_i->_data, 1);
  residue_mods(Aq_i);
}

/** \brief Convert between rnsMatrix_t with different bases */
void rnsMatrix_convertFancy(rnsMatrix_t * Aq, rnsMatrix_t const * Ap)
{
  long j;
  double * corr;
  assert(Aq->nrows == Ap->nrows && Aq->ncols == Ap->ncols);

  corr = correctionFactor(Ap);
  for (j = 0; j < Aq->nmod; ++j) {
    rnsMatrix_convertOne(Aq->residues[j], Ap, corr);
  }
  free(corr);
}

#ifndef THREAD
void rnsMatrix_convertSimple(rnsMatrix_t * Aq, rnsMatrix_t const * Ap)
{
  static mpzMatrix_t * A = NULL;
  if(!A) { A = mpzMatrix_init(Aq->nrows, Aq->ncols); }
  else if(Aq->nrows != A->nrows || Aq->ncols != A->ncols) {
    mpzMatrix_fini(A);
    A = mpzMatrix_init(Aq->nrows, Aq->ncols);
  }
  assert(Aq->nrows == Ap->nrows && Aq->ncols == Ap->ncols);

  mpzMatrix_reconstruct(A, Ap);
  rnsMatrix_fromMpzMatrix(Aq, A);
}

#else
/* TODO: cleaner threaded code? qv rns_matrix.c */


struct conv_args {
  residue_t * Aq_i;
  rnsMatrix_t const * Ap;
  residue_t const * r;
};
typedef struct conv_args conv_args_t;

static void * rnsMatrix_convertOne_thread(void * arg)
{
  conv_args_t * args = (conv_args_t *)arg;
  rnsMatrix_convertOne(args->Aq_i, args->Ap, args->r);
  return NULL;
}

void rnsMatrix_convert(rnsMatrix_t * Aq, rnsMatrix_t const * Ap)
{
  conv_args_t args[32];
  pthread_t thread[32];
  long j;

  residue_t * r = correctionFactor(Ap);

  for (j = 0; j < Aq->nmod; ++j) {
    args[j].Aq_i = Aq->residues[j];
    args[j].Ap = Ap;
    args[j].r = r;

    pthread_create(&thread[j], NULL, rnsMatrix_convertOne_thread, &(args[j]));
  }
  for (j = 0; j < Aq->nmod; ++j) {
    pthread_join(thread[j], NULL);
  }
}

#endif



