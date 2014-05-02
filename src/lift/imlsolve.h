#pragma once

#include "gmp.h"
#include "iml.h"

#include "mpz_matrix.h"
#include "rns_matrix.h"

/** \file imlsolve.h
 * \brief IML interface for `mpzMatrix_t`.
 */

void imlSolve(mpzMatrix_t * x, mpz_t d, mpzMatrix_t const * A, mpzMatrix_t const * b);
void imlSolveLeft(mpzMatrix_t * x, mpz_t d, mpzMatrix_t const * A, mpzMatrix_t const * b);

void imlSolveInv(mpzMatrix_t * x, mpz_t d, mpzMatrix_t const * A, mpzMatrix_t const * b, rnsMatrix_t const * Ainv_rns);
