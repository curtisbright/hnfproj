#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "gmp.h"

#include "dbg_print.h"
#include "timer.h"
#include "highorder.h"
#include "unicert_utils.h"

static void print_usage(void)
{
  printf("Usage: ./hor_main -n <matrix dimension> -l <entry bit length>\n");
}

static int parse_args(int argc, char ** argv, long * n, long * l, int * v)
{
  char * endptr = NULL;
  int c;

  *n = 0;
  *l = 0;
  *v = 0;

  while ((c = getopt(argc, argv, "l:n:v")) != -1) {
    switch (c) {
      case 'l':
        *l = strtol(optarg, &endptr, 10);
        break;
      case 'n':
        *n = strtol(optarg, &endptr, 10);
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
  mpzMatrix_t * A, * R;
  mpz_t Amax, Rmax;
  int check, rc;
  long n, l;

  if (!parse_args(argc, argv, &n, &l, &check)) {
    print_usage();
    return 0;
  }
  g_print_stream = stdout;
  g_timer_threshold = 3;
  g_timer_stream = stdout;

  A = mpzMatrix_init(n, n);
  mpzMatrix_rand(A, l);
  mpzMatrix_print(stdout, A);

  TIMER("TOTAL",
  R = highOrderResidue(A);)
  mpzMatrix_print(stdout, R);

  mpz_inits(Amax, Rmax, 0);
  mpzMatrix_max(Amax, A);
  mpzMatrix_max(Rmax, R);
  printf("log2 Amax: %ld\n", (long)mpz_sizeinbase(Amax, 2));
  printf("log2 Rmax: %ld\n", (long)mpz_sizeinbase(Rmax, 2));


  if (check) {
    TIMER("check",
    rc = hor_check(A, R);)
    printf("%s\n", rc ? "Success." : "Fail!");
  }

  mpzMatrix_fini(A);
  mpzMatrix_fini(R);
  mpz_clear(Rmax);

  return 0;
}
