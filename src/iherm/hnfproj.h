#pragma once

#include "gmp.h"

#include "pk_matrix.h"
#include "mpz_matrix.h"

/** \file hnfproj.h
  * Hermite form of a single projection.
  *
  * Given projection (X/d) = (1/A).v, for a random v, compute the Hermite form
  * H corresponding to X/d.  X is destroyed in-place.
  */

/** Balanced algorithm. */
void hermiteOfProjection_balanced(pk_t * H, mpzMatrix_t * X, mpz_t const dd);

/** Simple, iterative algorithm. */
void hermiteOfProjection_simple(pk_t * H, mpzMatrix_t * X, mpz_t const dd);
