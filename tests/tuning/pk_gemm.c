#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "gmp.h"

#include "highorder.h"
#include "mpz_matrix.h"
#include "pk_matrix.h"
#include "rns_matrix.h"
#include "timer.h"
#include "basis.h"



static int parse_args(int argc, char ** argv, long * n, long * j, long * k)
{
  char * endptr = NULL;
  int c;

  *n = 0;
  *j = 0;
  *k = 0;
  while ((c = getopt(argc, argv, "n:j:k:")) != -1) {
    switch (c) {
      case 'n':
        *n = strtol(optarg, &endptr, 10);
        break;
      case 'j':
        *j = strtol(optarg, &endptr, 10);
        break;
      case 'k':
        *k = strtol(optarg, &endptr, 10);
        break;
      default:
        break;
    }
  }
  return (*n != 0 && *j != 0 && *k != 0);
}

static pk_t * pkMatrix_rand(long n, long k)
{
  long i, j;
  mpzMatrix_t * A = mpzMatrix_init(n, n);
  pk_t * Z;
  mpz_t * data;
  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, rand());

  data = mpzMatrix_data(A);
  for (i = 0; i < A->nrows; ++i) {
    mpz_set_ui(data[i*A->nrows + i], 1);
  }

  for (i = 0; i < k; ++i) {
    long c = rand() % n;
    for (j = 0; j <= c; ++j) {
      mpz_urandomb(data[j*n+c], state, 10);
    }
  }

  gmp_randclear(state);

  Z = pkMatrix_fromFull(A);
  mpzMatrix_fini(A);

  return Z;
}

static int gemm_test(long n, long k1, long k2)
{
  pk_t *Z1, *Z2, *X, *Y;
  mpzMatrix_t *A, *B, *C1, *C2;

  printf("n: %ld\n", n);
  printf("k1: %ld\n", k1);
  printf("k2: %ld\n", k2);

  X = pkMatrix_rand(n, k1);
  Y = pkMatrix_rand(n, k2);
  Z1 = pkMatrix_init(n);
  Z2 = pkMatrix_init(n);

  /* compare straight iterative and blocked packed gemm */
  TIMER("gemm_iter",
  pkMatrix_gemm_iter(Z1, X, Y);)

  TIMER("gemm_block",
  pkMatrix_gemm_block(Z2, X, Y);)

  /* compare with standard dense gemm */
  A = pkMatrix_toFull(X);
  B = pkMatrix_toFull(Y);
  C1 = mpzMatrix_init(n, n);
  C2 = mpzMatrix_init(n, n);

  TIMER("mpz_gemm",
  mpzMatrix_gemm(C1, A, B);)

  TIMER("rns_gemm",
  mpzMatrix_rnsGemm(C2, A, B);)

  printf("%s\n", pkMatrix_equal(Z1, Z2) ? "Success." : "Fail!");

  pkMatrix_fini(X);
  pkMatrix_fini(Y);
  pkMatrix_fini(Z1);
  pkMatrix_fini(Z2);
  mpzMatrix_fini(A);
  mpzMatrix_fini(B);
  mpzMatrix_fini(C1);
  mpzMatrix_fini(C2);

  return 0;
}

int main(int argc, char ** argv)
{
  long n, j,k;

  g_timer_threshold = 1;
  g_timer_stream = stdout;

  srand(time(0));
  if (!parse_args(argc, argv, &n, &j, &k)) {
    return 0;
  }

  gemm_test(n, j, k);

  return 0;
}
