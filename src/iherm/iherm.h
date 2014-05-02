#pragma once

#include "mpz_matrix.h"

/** \file iherm.h
  *  Hermite normal form entry points.
  */

/** Square, nonsingular Hermite normal form.*/
mpzMatrix_t * hermite(mpzMatrix_t const * A);
mpzMatrix_t * hermiteWithTransform(mpzMatrix_t const * A, mpzMatrix_t ** U);

mpzMatrix_t * hermiteRect(mpzMatrix_t const * A);
