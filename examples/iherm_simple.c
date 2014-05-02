#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "gmp.h"

#include "libhnfproj.h"

static void print_usage(char const * a)
{
  printf("Usage: %s -n <matrix dimension> -l <entry bit length>\n", a);
}

static int parse_args(int argc, char ** argv, long * n, long * l)
{
  char * endptr = NULL;
  int c;

  *n = 0;
  *l = 3;

  while ((c = getopt(argc, argv, "l:n:")) != -1) {
    switch (c) {
      case 'l':
        *l = strtol(optarg, &endptr, 10);
        break;
      case 'n':
        *n = strtol(optarg, &endptr, 10);
        break;
      default:
        break;
    }
  }
  return (*n != 0);
}

static void print_diag(mpz_t * A, long n)
{
  long i;
  for (i = 0; i < n; ++i) {
    mpz_out_str(stdout, 10, A[i*n + i]);
    printf(", ");
  }
  printf("\n");
}

int main(int argc, char ** argv)
{
  long * A;
  mpz_t * H;
  long i, n, l;

  if (!parse_args(argc, argv, &n, &l)) {
    print_usage(argv[0]);
    return 0;
  }

  A = malloc(sizeof(long)*n*n);
  for (i = 0; i < n * n; ++i) {
    A[i] = rand() % (1L << l);
  }

  H = hermiteLong(A, n);
  print_diag(H, n);

  free(A);

  return 0;
}
