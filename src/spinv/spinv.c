#include "spinv.h"

#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "dbg_print.h"
#include "lift.h"
#include "reconstruct.h"
#include "timer.h"

spinv_t initSpinv(long n, long k)
{
  long i;
  spinv_t s;
  s.n = n;
  s.k = k;
  s.R = malloc(sizeof(mpzMatrix_t*) * (k+1));
  s.M = malloc(sizeof(mpzMatrix_t*) * (k));
  for (i = 0; i < k; ++i) {
    s.R[i] = mpzMatrix_init(n, n);
    s.M[i] = mpzMatrix_init(n, n);
  }
  s.R[k] = mpzMatrix_init(n, n);
  s.C = mpzMatrix_init(n, n);

  mpz_init_set_ui(s.X, 1);
  return s;
}

void finiSpinv(spinv_t * s)
{
  long i;
  for (i = 0; i < s->k; ++i) {
    mpzMatrix_fini(s->R[i]);
    mpzMatrix_fini(s->M[i]);
  }
  mpzMatrix_fini(s->R[s->k]);
  mpzMatrix_fini(s->C);
  free(s->R);
  free(s->M);
  mpz_clear(s->X);
}

static spinv_t sparseInv(mpzMatrix_t const * A, long k_req)
{
  long i;
  lift_info_t * info = initLift(A);
  long k_max = numLiftIters(A, liftInfo_modulus(info));
  long k = (k_req > k_max) ? k_max : k_req;
  spinv_t spinv;

  spinv = initSpinv(A->nrows, k);

  mpz_set(spinv.X, liftInfo_modulus(info));
  mpzMatrix_reconstruct(spinv.C, liftInfo_inverse(info));
  for (i = 0; i < k ; ++i) {
    dprintf(3, "%ld...", 1+i);
    info = lift(info);
    mpzMatrix_reconstruct(spinv.R[i], liftInfo_R(info));
    mpzMatrix_reconstruct(spinv.M[i], liftInfo_M(info));
  }

  info = lift(info);
  mpzMatrix_reconstruct(spinv.R[k], liftInfo_R(info));
  finiLift(info);

  return spinv;
}

spinv_t sparseInvLong(long const * A, long n, long k)
{
  spinv_t rslt;
  mpzMatrix_t * A_mpz = mpzMatrix_initSet(n, n, A);

  rslt = sparseInvMpz(A_mpz, k);
  mpzMatrix_fini(A_mpz);

  return rslt;
}

spinv_t sparseInvMpz(mpzMatrix_t const * A, long k)
{
  return sparseInv(A, k);
}
