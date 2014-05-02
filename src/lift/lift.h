#pragma once

#include <inttypes.h>

#include "gmp.h"
#include "mpz_matrix.h"
#include "rns_matrix.h"

/** \file lift.h
 * \brief Double-plus-one lifting.
 */

typedef struct lift_info_t lift_info_t;

lift_info_t * initLift(mpzMatrix_t const * A);
lift_info_t * lift(lift_info_t * info);
void finiLift(lift_info_t * info);

long numLiftIters(mpzMatrix_t const * A, mpz_t const XX);

mpz_srcptr liftInfo_modulus(lift_info_t const * info);
rnsMatrix_t * liftInfo_adoptInverse(lift_info_t * info);
rnsMatrix_t * liftInfo_adoptR(lift_info_t * info);
rnsMatrix_t const * liftInfo_inverse(lift_info_t const * info);
rnsMatrix_t const * liftInfo_R(lift_info_t const * info);
rnsMatrix_t const * liftInfo_M(lift_info_t const * info);

