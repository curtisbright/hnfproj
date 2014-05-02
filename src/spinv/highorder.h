#pragma once

#include "rns_matrix.h"
#include "mpz_matrix.h"

/** \file highorder.h
 * \brief High order residue computation.
 */

mpzMatrix_t * highOrderResidue_mpz(mpzMatrix_t const * A, rnsMatrix_t ** Ainv);
rnsMatrix_t * highOrderResidue_rns(mpzMatrix_t const * A, rnsMatrix_t ** Ainv);

mpzMatrix_t * highOrderResidue(mpzMatrix_t const * A);


