#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "gmp.h"

#include "horsolve.h"
#include "imlsolve.h"
#include "mpz_matrix.h"
#include "spinv.h"
#include "timer.h"

static void print_usage(void)
{
  printf("Usage: ./linsys_main -n <system dimension> -l <entry bit length>\n");
}

static int parse_args(int argc, char ** argv, long * n, long * l, long * v, long * m)
{
  char * endptr = NULL;
  int c;

  *m = 1;
  *n = 0;
  *l = 0;
  *v = 0;
  while ((c = getopt(argc, argv, "l:n:m:v")) != -1) {
    switch (c) {
      case 'l':
        *l = strtol(optarg, &endptr, 10);
        break;
      case 'n':
        *n = strtol(optarg, &endptr, 10);
        break;
      case 'm':
        *m = strtol(optarg, &endptr, 10);
        break;
      case 'v':
        *v = 1;
        break;
      default:
        break;
    }
  }
  return (*n != 0 && *l != 0);
}

int main(int argc, char ** argv)
{
  mpzMatrix_t *A, *b;
  mpzMatrix_t * x[4];
  mpz_t d[4];
  long i, n, m, l, v;
  int eq, Rzero;

  if (!parse_args(argc, argv, &n, &l, &v, &m)) {
    print_usage();
    return 0;
  }

  A = mpzMatrix_init(n, n);
  b = mpzMatrix_init(n, m);
  mpzMatrix_rand(A, l);
  mpzMatrix_rand(b, l);

  for (i = 0; i < 2; ++i) {
    mpz_init(d[i]);
  }

  printf("----------- %4ld, %4ld -----------\n", n, l);

  TIMER("w/ horSolve",
  x[0] = horSolveSpinv(d[0], &Rzero, A, b);)

  TIMER("w/ IML",
  x[1] = horSolveIML(d[1], &Rzero, A, b);)

  printf("=====\n");


  eq = mpzMatrix_equal(x[0], x[1]);
  if (!eq) {
    printf("Fail!\n");
  } else {
    printf("Success.\n");
  }

  if(v){
  }

  mpzMatrix_fini(A);
  mpzMatrix_fini(b);
  for (i = 0; i < 2; ++i) {
    mpzMatrix_fini(x[i]);
    mpz_clear(d[i]);
  }
  return 0;
}
