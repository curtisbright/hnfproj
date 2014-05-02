#pragma once

#include "mpz_matrix.h"
#include "libhnfproj.h"

/** \file unicert.h
 * \brief Unimodularity certification.
 */

int uniCert(mpzMatrix_t const * A);
