#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "gmp.h"

#include "linsys_utils.h"
#include "spinvsolve.h"
#include "mpz_matrix.h"
#include "timer.h"
#if 0
static void print_usage(void)
{
  printf("Usage: ./linsys_main -n <system dimension> -l <entry bit length>\n");
}

static int parse_args(int argc, char ** argv, long * n, long * l, long * v)
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


#include "dbg_print.h"
int g_timer_depth = 0;
int g_print_level = 2;
FILE * g_print_stream;

int main(int argc, char ** argv)
{
  mpzMatrix_t *A, *b, *x;
  mpz_t d;
  long n, l, v;
  int rc;
  mpz_init(d);

  if (!parse_args(argc, argv, &n, &l, &v)) {
    print_usage();
    return 0;
  }

  g_print_stream = stdout;

  x = mpzMatrix_init(n, 1);


  TIMER("INPUT",
  linsys_input(n, l, 0, &A, &b);)


  TIMER("SOLVE",
  x = solve_spinv(d, A, b);)

  TIMER("iml",
  imlSolve(x, d, A, b);)

  TIMER("CHECK",
  rc = linsys_check(A, b, x, d);)

  printf("%s\n", rc ? "Success." : "Fail!");

  mpzMatrix_fini(A);
  mpzMatrix_fini(b);
  mpzMatrix_fini(x);

  return 0;
}
#endif

int main(int argc, char ** argv){}
