#include "imlsolve.h"

#include <stdlib.h>

#include "iml.h"

void imlSolve(mpzMatrix_t * x, mpz_t d, mpzMatrix_t const * A, mpzMatrix_t const * b)
{
  nonsingSolvLlhsMM(RightSolu, b->nrows, b->ncols, A->data, b->data, x->data, d);
}

void imlSolveLeft(mpzMatrix_t * x, mpz_t d, mpzMatrix_t const * A, mpzMatrix_t const * b)
{
  nonsingSolvLlhsMM(LeftSolu, b->ncols, b->nrows, A->data, b->data, x->data, d);
}

void imlSolveInv(mpzMatrix_t * x, mpz_t d, mpzMatrix_t const * A, mpzMatrix_t const * b, rnsMatrix_t const * Ainv_rns)
{
  long i, nX;
  unsigned long *X;
  double **Ainv;

  if (!Ainv_rns) {
    nonsingSolvLlhsMM(RightSolu, b->nrows, b->ncols, A->data, b->data, x->data, d);
    return;
  }

  nX = Ainv_rns->nmod;
  Ainv = malloc(nX * sizeof(double*));
  X = malloc((1+nX) * sizeof(long));

  for (i = 0; i < nX; ++i) {
    Ainv[i] = Ainv_rns->residues[i]->_data;
    X[i] = basis_getPrime(Ainv_rns->basis, i);
  }
  X[nX] = 0;

  nonsingSolvLlhsMM_X(RightSolu, b->nrows, b->ncols, A->data, b->data, x->data, d, Ainv, X);
  free(X);
  free(Ainv);
}
