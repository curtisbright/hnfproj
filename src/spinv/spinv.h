#pragma once

#include "gmp.h"

#include "mpz_matrix.h"

/** \file spinv.h
 * \brief Sparse inverse expansion computation.
 */

struct spinv {
  long n;
  long k;
  mpzMatrix_t ** R;
  mpzMatrix_t ** M;
  mpzMatrix_t *C;
  mpz_t X;
};
typedef struct spinv spinv_t;

spinv_t initSpinv(long n, long k);
void finiSpinv(spinv_t * s);

spinv_t sparseInvLong(long const * A, long n, long k);
spinv_t sparseInvMpz(mpzMatrix_t const * A, long k);
