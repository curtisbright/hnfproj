#include "lift.h"

#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cblas.h"
#include "iml.h"
#include "gmp.h"

#include "arith_utils.h"
#include "dbg_print.h"
#include "reconstruct.h"
#include "rns_conversion.h"
#include "rns_matrix.h"
#include "mpz_matrix.h"

struct lift_tmp {
  long * Xinv;
  rnsMatrix_t * Ap;

  rnsMatrix_t * Mp;
  rnsMatrix_t * Rx;
  rnsMatrix_t * Tp;
  rnsMatrix_t * Tx;
};
typedef struct lift_tmp lift_tmp_t;

struct lift_info_t {
  long n;
  basis_t * P;
  basis_t * X;

  lift_tmp_t tmp;

  rnsMatrix_t * Cx;
  rnsMatrix_t * Rp;
  rnsMatrix_t * Mx;
};

static void liftTmp_init(lift_tmp_t * tmp, long n, basis_t const * P, basis_t const * X)
{
  tmp->Xinv = malloc((P->nmod)*sizeof(long));
  tmp->Ap = rnsMatrix_init(n, n, P);
  tmp->Mp = rnsMatrix_init(n, n, P);
  tmp->Rx = rnsMatrix_init(n, n, X);
  tmp->Tp = rnsMatrix_init(n, n, P);
  tmp->Tx = rnsMatrix_init(n, n, X);
}

static void liftTmp_fini(lift_tmp_t * tmp)
{
  free(tmp->Xinv);
  rnsMatrix_fini(tmp->Ap);
  rnsMatrix_fini(tmp->Mp);
  rnsMatrix_fini(tmp->Rx);
  rnsMatrix_fini(tmp->Tp);
  rnsMatrix_fini(tmp->Tx);
}


/* lift_info_t accessors */
mpz_srcptr liftInfo_modulus(lift_info_t const * info)
{
  return info->X->PP;
}
rnsMatrix_t const * liftInfo_inverse(lift_info_t const * info)
{
  return info->Cx;
}
rnsMatrix_t * liftInfo_adoptInverse(lift_info_t * info)
{
  rnsMatrix_t * inv = info->Cx;
  info->Cx = NULL;
  return inv;
}
rnsMatrix_t const * liftInfo_R(lift_info_t const * info)
{
  return info->Rp;
}
rnsMatrix_t * liftInfo_adoptR(lift_info_t * info)
{
  rnsMatrix_t * R = info->Rp;
  info->Rp = NULL;
  return R;
}
rnsMatrix_t const * liftInfo_M(lift_info_t const * info)
{
  return info->Mx;
}


static void calcBounds(mpzMatrix_t const * A, mpz_t P, mpz_t X)
{
  long n = A->nrows;
  mpz_t S;

  /* S = ||A|| */
  mpz_init(S);
  mpzMatrix_max(S, A);

  /* P = 2*0.6001*n*S */
  mpz_mul_ui(P, S, n);
  mpz_mul_ui(P, P, 6001);
  mpz_cdiv_q_ui(P, P, 5000);

  /*X = 3.61*n*n*S */
  mpz_mul_ui(X, S, n);
  mpz_mul_ui(X, X, n);
  mpz_mul_ui(X, X, 361);
  mpz_cdiv_q_ui(X, X, 100);

  mpz_clear(S);
}


static void calcXinv(basis_t const * P, basis_t const * X, long * Xinv)
{
  long i;
  for (i = 0; i < P->nmod; ++i) {
    Xinv[i] = modInverseMpz(X->PP, basis_getPrime(P, i));
  }
}

static void setupBasis(lift_info_t * info, mpzMatrix_t const * A, long p0)
{
  mpz_t Pbound, Xbound;
  mpz_inits(Pbound, Xbound, 0);

  calcBounds(A, Pbound, Xbound);

  info->P = basis_initFromBound(p0, Pbound);
  info->X = basis_initFromBound(basis_getPrime(info->P,info->P->nmod-1), Xbound);

  dprintf(3, "X: ");
  dprint2(3, basis_print, info->X);
  dprintf(3, "P: ");
  dprint2(3, basis_print, info->P);

  mpz_clears(Pbound, Xbound, 0);
}

static int setupLift(lift_info_t * info, mpzMatrix_t const * A)
{
  int rc;
  long n = info->n;
  basis_t const * P = info->P;
  basis_t const * X = info->X;
  lift_tmp_t * tmp = &(info->tmp);

  info->Cx = rnsMatrix_init(n, n, X);
  info->Rp = rnsMatrix_init(n, n, P);
  info->Mx = rnsMatrix_init(n, n, X);
  liftTmp_init(&(info->tmp), n, P, X);

  /* Reduce A mod P, A mod X */
  rnsMatrix_fromMpzMatrix(tmp->Ap, A);
  rnsMatrix_fromMpzMatrix(info->Cx, A);

  /* C = A^-1 mod X */
  rc = rnsMatrix_inverse(info->Cx);
  if (!rc) { return 0; }

  /* Xinv = 1/X mod P */
  calcXinv(P, X, tmp->Xinv);

  /* R[0] = I */
  rnsMatrix_identity(info->Rp);

  /* M[0] = A^-1 mod X */
  rnsMatrix_copy(info->Mx, info->Cx);
  rnsMatrix_convert(tmp->Mp, info->Mx);
  return 1;
}


lift_info_t * initLift(mpzMatrix_t const * A)
{
  int rc;
  long p0;
  lift_info_t * info = malloc(sizeof(lift_info_t));
  info->n = A->nrows;

  p0 = pickStartModulus(A->nrows);
  do {
    setupBasis(info, A, p0);
    rc = setupLift(info, A);
    p0 = basis_getPrime(info->P, 0);
  } while (!rc);

  return info;
}

lift_info_t * lift(lift_info_t * info)
{
  lift_tmp_t * tmp = &(info->tmp);

  /* T = R[i].R[i] */
  rnsMatrix_gemm(tmp->Tp, info->Rp, info->Rp);

  /* R[i] = Xinv.(T - A.M[i-1]) mod P */
  rnsMatrix_quadLift(info->Rp, tmp->Tp, tmp->Ap, tmp->Mp, tmp->Xinv);
  rnsMatrix_convert(tmp->Rx, info->Rp);

  /* T = R[i].R[i] */
  rnsMatrix_gemm(tmp->Tx, tmp->Rx, tmp->Rx);

  /* M[i] = C.T mod X */
  rnsMatrix_gemm(info->Mx, info->Cx, tmp->Tx);
  rnsMatrix_convert(tmp->Mp, info->Mx);

  return info;
}

void finiLift(lift_info_t * info)
{
  basis_fini(info->P);
  basis_fini(info->X);

  rnsMatrix_fini(info->Cx);
  rnsMatrix_fini(info->Rp);
  rnsMatrix_fini(info->Mx);
  liftTmp_fini(&(info->tmp));
  free(info);
}

#define MAX(a, b) (((a) >= (b)) ? (a) : (b))

long numLiftIters(mpzMatrix_t const * A, mpz_t const XX)
{
  long n = A->nrows;
  long k = 1;
  mpz_t bound, S, X, Y;

  mpz_init(bound);
  mpz_init(S);
  mpz_init_set(X, XX);
  mpz_init_set_ui(Y, 1);

  /* bound := n^((n-1)/2)*||A||^(n-1)  /  n^2||A|| */
  mpzMatrix_max(S, A);
  mpz_pow_ui(S, S, MAX(n-2, 1));
  mpz_ui_pow_ui(bound, n, MAX(n/2-2, 1));
  mpz_mul(bound, bound, S);

  while(mpz_cmp(Y, bound) < 0) {
    mpz_mul(X, X, X);
    mpz_mul(Y, Y, X);
    ++k;
  }

  mpz_clear(bound);
  mpz_clear(S);
  mpz_clear(X);
  mpz_clear(Y);

  return k;
}
