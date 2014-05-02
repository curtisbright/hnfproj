#include "unicert.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "lift.h"
#include "rns_matrix.h"

int uniCert(mpzMatrix_t const * A)
{
  int i = 0;
  int rslt = 0;

  lift_info_t * info = initLift(A);
  long k = numLiftIters(A, liftInfo_modulus(info));
  for (i = 0; (i < k && !rslt); ++i) {
    info = lift(info);
    rslt = rnsMatrix_isZero(liftInfo_R(info));
  }
  finiLift(info);

  return rslt;
}

int uniCertLong(long const * A, long n)
{
  int rslt;
  mpzMatrix_t * A_mpz = mpzMatrix_initSet(n, n, A);

  rslt = uniCert(A_mpz);

  mpzMatrix_fini(A_mpz);
  return rslt;
}

int uniCertMpz(mpz_t const * A, long n)
{
  int rslt;
  mpzMatrix_t const * A_mpz;

  A_mpz = mpzMatrix_initFromMpz(n, n, A);
  rslt = uniCert(A_mpz);
  free((void*)A_mpz);

  return rslt;
}
