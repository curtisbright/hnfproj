#include "reconstruct.h"

#include <assert.h>
#include <stdlib.h>

#include "arith_utils.h"
#include "basis.h"


void mpzMatrix_gemmRnsMpz(mpzMatrix_t * dst, rnsMatrix_t const * A, mpzMatrix_t const * B)
{
  rnsMatrix_t * dst_rns = rnsMatrix_init(B->nrows, B->ncols, A->basis);
  rnsMatrix_t * B_rns = rnsMatrix_init(B->nrows, B->ncols, A->basis);

  rnsMatrix_fromMpzMatrix(B_rns, B);
  rnsMatrix_gemm(dst_rns, A, B_rns);
  mpzMatrix_reconstruct(dst, dst_rns);

  rnsMatrix_fini(dst_rns);
  rnsMatrix_fini(B_rns);
}

static void mpzMatrix_addResidue(mpzMatrix_t * A, residue_t const * B, mpz_t const k)
{
  /* A += k*B */
  long i, b;
  for (i = 0; i < mpzMatrix_numElems(A); ++i) {
    b = residue_getEntry(B, i);
    addmul_si(A->data[i], k, b);
  }
}

void mpzMatrix_reconstruct(mpzMatrix_t * A, rnsMatrix_t const * Ap)
{
  basis_t const * Px = Ap->basis;
  long nP = Px->nmod;
  long i;

  assert(Ap->nrows == A->nrows);
  assert(Ap->ncols == A->ncols);
  assert(Px->nmod == Ap->nmod);

  /* A:= sum( C_i * Ap_i ) mod PP*/
  mpzMatrix_zero(A);
  for (i = 0; i < nP; ++i) {
    mpzMatrix_addResidue(A, Ap->residues[i], Px->C[i]);
  }
  mpzMatrix_mods(A, Px->PP);
}

void rns_reconstruct(mpz_t ret, basis_t const * P, long * d)
{
  long j;
  mpz_set_ui(ret, 0);
  for (j = 0; j < P->nmod; ++j) {
    addmul_si(ret, P->C[j], d[j]);
  }
  modsMpzMpz(ret, P->PP);
}

void rnsMatrix_fromMpzMatrix(rnsMatrix_t * Ap, mpzMatrix_t const * A)
{
  long j;
  long nP = Ap->nmod;

  assert(Ap->nrows == A->nrows);
  assert(Ap->ncols == A->ncols);

  for (j = 0; j < nP; ++j) {
    residue_fromMpzMatrix(Ap->residues[j], A);
  }
}
