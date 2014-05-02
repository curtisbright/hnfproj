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

static int parse_args(int argc, char ** argv, long * m, long * n, long * o, long * l1, long * l2)
{
  char * endptr = NULL;
  int c;
  long * arg = NULL;

  *m = 0;
  *n = 0;
  *o = 0;
  *l1 = 3;
  *l2 = 0;
  while ((c = getopt(argc, argv, "m:n:o:p:q:")) != -1) {
    switch (c) {
      case 'm': arg = m; break;
      case 'n': arg = n; break;
      case 'o': arg = o; break;
      case 'p': arg = l1; break;
      case 'q': arg = l2; break;
      default:
        break;
    }
    *arg = strtol(optarg, &endptr, 10);
  }
  if (*m == 0 && *o == 0) { *m = *n; *o = *n; }
  if (*l2 == 0) { *l2 = *l1; }

  return (*n != 0);
}

static int gemm_test(long m, long n, long o, long l1, long l2)
{
  mpzMatrix_t *C1, *A, *B;
  mpzMatrix_t *C2;

  printf("m: %ld\n", m);
  printf("n: %ld\n", n);
  printf("o: %ld\n", o);
  printf("l1: %ld\n", l1);
  printf("l2: %ld\n", l2);

  C1 = mpzMatrix_init(m ,o);
  C2 = mpzMatrix_init(m, o);
  A = mpzMatrix_init(m, n);
  B = mpzMatrix_init(n, o);
  mpzMatrix_rand(A, l1);
  mpzMatrix_rand(B, l2);

  TIMER("rns",
  mpzMatrix_rnsGemm(C2, A, B);)

  TIMER("mpz",
  mpzMatrix_gemm(C1, A, B);)

  printf("%s\n", mpzMatrix_equal(C1, C2) ? "Success." : "Fail!");

  mpzMatrix_fini(A);
  mpzMatrix_fini(B);
  mpzMatrix_fini(C1);
  mpzMatrix_fini(C2);

  return 0;
}

int main(int argc, char ** argv)
{
  long m, n, o, l1, l2;

  g_timer_stream = stdout;
  g_timer_threshold = 1;

  if (!parse_args(argc, argv, &m, &n, &o, &l1, &l2)) {
    return 0;
  }

  gemm_test(m, n, o, l1, l2);
  return 0;
}

