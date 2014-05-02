#pragma once

#include "mpz_matrix.h"

/** \file iherm_utils.h
  * Utility routines for HNF input, verification.
  */

int iherm_check(mpzMatrix_t const * A, mpzMatrix_t * H);
mpzMatrix_t * iherm_input(char const * type, long nrows, long ncols, long rank, long bitlen);
