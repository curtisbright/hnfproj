#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "gmp.h"

#include "arith_utils.h"
#include "reconstruct.h"
#include "rns_conversion.h"
#include "timer.h"
#include "mpz_matrix.h"
#include "rns_matrix.h"

static int parse_args(int argc, char ** argv, long * n, long * l, long * k)
{
  char * endptr = NULL;
  int c;

  *n = 0;
  *l = 22;
  *k = 0;
  while ((c = getopt(argc, argv, "l:n:k:")) != -1) {
    switch (c) {
      case 'l':
        *l = strtol(optarg, &endptr, 10);
        break;
      case 'n':
        *n = strtol(optarg, &endptr, 10);
        break;
      case 'k':
        *k = strtol(optarg, &endptr, 10);
        break;
      default:
        break;
    }
  }
  return (*n != 0 && *k != 0);
}

static void rnsConv_test(long n, long k)
{
  long i, l;
  mpzMatrix_t *A;
  mpzMatrix_t *chk;
  rnsMatrix_t *Ap, *Aq, *Aq2;

  long * lP = malloc(k*sizeof(long));
  long * lQ = malloc(k*sizeof(long));
  basis_t *P, *Q;
  long start = pickStartModulus(n);

  lP[0] = start;
  for(i = 1; i < k; ++i) {
    lP[i] = prevprime(lP[i-1]);
  }
  lQ[0] = prevprime(lP[k-1]);
  for(i = 1; i < k; ++i) {
    lQ[i] = prevprime(lQ[i-1]);
  }

  P = basis_init(lP, k);
  Q = basis_init(lQ, k);
  l = mpz_sizeinbase(P->PP, 2) - 2;

  A = mpzMatrix_init(n, n);
  chk = mpzMatrix_init(n, n);

  Ap = rnsMatrix_init(n, n, P);
  Aq = rnsMatrix_init(n, n, Q);
  Aq2 = rnsMatrix_init(n, n, Q);

  mpzMatrix_rand(A, l);

  printf("n: %ld\n", n);
  printf("k: %ld\n", k);
  printf("l: %ld\n", l);
  fflush(stdout);

  mpzMatrix_mods(A, P->PP);

  TIMER("nmod",
  rnsMatrix_fromMpzMatrix(Ap, A);)

  TIMER("conv_fancy",
  rnsMatrix_convertFancy(Aq, Ap);)

  TIMER("conv_simple",
  rnsMatrix_convertSimple(Aq2, Ap);)

  TIMER("recon",
  mpzMatrix_reconstruct(chk, Aq);)

  if(mpzMatrix_equal(A, chk)) {
    printf("Success.\n");
  } else {
    printf("Fail!\n");
  }
  fflush(stdout);

  rnsMatrix_fini(Ap);
  rnsMatrix_fini(Aq);
  mpzMatrix_fini(A);
  mpzMatrix_fini(chk);

  basis_fini(P);
  basis_fini(Q);
  free(lP);
  free(lQ);
}



int main(int argc, char ** argv)
{
  long n, l,k;

  g_timer_threshold = 1;
  g_timer_stream = stdout;

  srand(time(0));
  if (!parse_args(argc, argv, &n, &l, &k)) {
    return 0;
  }

  rnsConv_test(n, k);

  return 0;

}
