#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "gmp.h"

#include "dbg_print.h"
#include "iherm.h"
#include "iherm_utils.h"
#include "timer.h"

static void print_usage(void)
{
  printf("Usage: ./det_main -n <dimension> [-t 0-6] [-v] [-l bit length]\n");
}

static int parse_args(int argc, char ** argv, long * n, long * t, long * v, long * l)
{
  char * endptr = NULL;
  int c;

  *n = 0;
  *t = 1;
  *v = 0;
  *l = 3;
  while ((c = getopt(argc, argv, "n:t:vl:")) != -1) {
    switch (c) {
      case 'n':
        *n = strtol(optarg, &endptr, 10);
        break;
      case 't':
        *t = strtol(optarg, &endptr, 10);
        break;
      case 'v':
        *v = 1;
        break;
      case 'l':
        *l = strtol(optarg, &endptr, 10);
        break;
      default:
        break;
    }
  }
  return (*n != 0);
}

static void mul_diag(mpz_t det, mpzMatrix_t const * A)
{
  long i;
  long n = A->nrows;
  mpz_set_ui(det, 1);
  for (i = 0; i < n; ++i) {
    mpz_mul(det, det, A->data[i*n + i]);
  }
}

int g_print_level = 2;
FILE * g_print_stream;


int main(int argc, char ** argv)
{
  mpzMatrix_t *A, *H;
  mpz_t det, chk;
  long n, t, v, l;
  int rc;

  if (!parse_args(argc, argv, &n, &t, &v, &l)) {
    print_usage();
    return 0;
  }
  g_print_stream = stdout;

  srand(0);

  A = iherm_input(n, t, l);

  mpz_inits(det, chk, 0);
  TIMER("CRA determinant",
  mpzMatrix_determinant(det, A);)

  TIMER("HNF",
  H = hermite(A);)

  mul_diag(chk, H);
  rc = (0 == mpz_cmpabs(chk, det));

  printf("%s\n", rc ? "Success." : "Fail!");
  if (v) {
    mpz_out_str(stdout, 10, det); printf("\n");
  }

  mpzMatrix_fini(A);
  mpzMatrix_fini(H);
  mpz_clears(det, chk, 0);
  return 0;
}
