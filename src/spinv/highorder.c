#include "libhnfproj.h"
#include "highorder.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <gmp.h>

#include "dbg_print.h"
#include "lift.h"
#include "reconstruct.h"

static lift_info_t * horBase(mpzMatrix_t const * A, rnsMatrix_t ** Ainv)
{
  int i;
  lift_info_t * info = initLift(A);
  long k = numLiftIters(A, liftInfo_modulus(info));

  for (i = 0; i <= k; ++i) {
    dprintf(3, "%d...", i);
    info = lift(info);
    if (rnsMatrix_isZero(liftInfo_R(info))) { break; }
  }

  if (Ainv) {
    *Ainv = liftInfo_adoptInverse(info);
  }

  return info;
}

mpzMatrix_t * highOrderResidue_mpz(mpzMatrix_t const * A, rnsMatrix_t ** Ainv)
{
  mpzMatrix_t * rslt;
  lift_info_t * info;
  info = horBase(A, Ainv);

  rslt = mpzMatrix_init(A->nrows, A->ncols);
  mpzMatrix_reconstruct(rslt, liftInfo_R(info));
  finiLift(info);

  return rslt;
}

rnsMatrix_t * highOrderResidue_rns(mpzMatrix_t const * A, rnsMatrix_t ** Ainv)
{
  lift_info_t * info = horBase(A, Ainv);
  rnsMatrix_t * R = liftInfo_adoptR(info);
  finiLift(info);

  return R;
}

mpzMatrix_t * highOrderResidue(mpzMatrix_t const * A)
{
  return highOrderResidue_mpz(A, NULL);
}

long * highOrderResidueLong(long const * A, long n)
{
  int i;
  long * rslt;
  mpzMatrix_t * A_mpz = mpzMatrix_initSet(n, n, A);
  mpzMatrix_t * tmp = highOrderResidue(A_mpz);

  rslt = malloc(mpzMatrix_numElems(tmp) * sizeof(long));

  for (i = 0; i < mpzMatrix_numElems(tmp); ++i) {
    assert(mpz_fits_slong_p(mpzMatrix_data(tmp)[i]));
    rslt[i] = mpz_get_si(mpzMatrix_data(tmp)[i]);
  }

  mpzMatrix_fini(A_mpz);
  mpzMatrix_fini(tmp);

  return rslt;
}

mpz_t * highOrderResidueMpz(mpz_t const * A, long n)
{
  mpzMatrix_t * rslt;
  mpzMatrix_t const * A_mpz;

  A_mpz = mpzMatrix_initFromMpz(n, n, A);
  rslt = highOrderResidue(A_mpz);
  free((void*)A_mpz);

  return mpzMatrix_data(rslt);
}
