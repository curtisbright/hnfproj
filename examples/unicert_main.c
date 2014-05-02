#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "gmp.h"

#include "dbg_print.h"
#include "timer.h"
#include "unicert.h"
#include "unicert_utils.h"

static void print_usage(void)
{
  printf("Usage: ./unicert_main -n <matrix dimension> -l <entry bit length> [-u]\n");
}

static int parse_args(int argc, char ** argv, long * n, long * l, int * u)
{
  char * endptr = NULL;
  int c;

  *n = 0;
  *l = 0;
  *u = 0;

  while ((c = getopt(argc, argv, "l:n:u")) != -1) {
    switch (c) {
      case 'l':
        *l = strtol(optarg, &endptr, 10);
        break;
      case 'n':
        *n = strtol(optarg, &endptr, 10);
        break;
      case 'u':
        *u = 1;
        break;
      default:
        break;
    }
  }
  return (*n != 0 && *l != 0);
}

int main(int argc, char ** argv)
{
  mpzMatrix_t * A = NULL;
  int uni, u, rc;
  long n, l;

  if (!parse_args(argc, argv, &n, &l, &u)) {
    print_usage();
    return 0;
  }
  g_print_stream = stdout;
  g_timer_threshold = 3;
  g_timer_stream = stdout;

  A = unicert_input(n, l, u);

  TIMER("Total",
  uni = uniCert(A);)
  printf("Unimodular?  %s\n", uni ? "Yes." : "No.");

  rc = unicert_check(A, uni);
  printf("%s\n", rc ? "Success." : "Fail!");

  mpzMatrix_fini(A);

  return 0;
}
