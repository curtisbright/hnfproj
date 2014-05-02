#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gmp.h"

#include "dbg_print.h"
#include "iherm.h"
#include "iherm_utils.h"
#include "timer.h"

/** \file iherm_main.c
  * Hermite normal form demonstration program.
  */

/** X-macro for commmand line arguments
  * ARG(short, long, type, default, flag, func)
  */
#define ALL_ARGS \
  ARG('n', nrows,                long,  0,        required_argument, do_num_arg) \
  ARG('m', ncols,                long,  0,        required_argument, do_num_arg) \
  ARG('r', rank,                long,  0,        required_argument, do_num_arg) \
  ARG('t', type,         char const *,  "random", required_argument, do_str_arg) \
  ARG('l', bitlen,               long,  8,        required_argument, do_num_arg) \
  ARG('d', print_diag,            int,  0,        no_argument,       do_flag_arg) \
  ARG('v', verify,                int,  0,        no_argument,       do_flag_arg) \
  ARG('h', timer_threshold,       int,  4,        required_argument, do_num_arg) \
  ARG('o', out_filename, char const *,  0,        required_argument, do_str_arg) \
  ARG('i', in_filename,  char const *,  0,        required_argument, do_str_arg) \
  ARG('S', seed,         unsigned int,  0,        required_argument, do_num_arg) \
  ARG('R', rect,                  int,  0,        no_argument, do_flag_arg) \

static void print_usage(void)
{
  printf("Usage: ./iherm_main -n <dimension>  [-l <bitlength>] [-t random|jaeger|steel]\n");
  printf("Options:\n");
  printf("\t-n, --nrows N\n");
  printf("\t      dimension of input matrix\n");
  printf("\t-l, --bitlen L\n");
  printf("\t      entry bit length\n");
  printf("\t-t, --type {random|jaeger|steel}\n");
  printf("\t      type of input matrix\n");
  printf("\t       - random, random entries\n");
  printf("\t       - jaeger, A_ij = i^j mod n\n");
  printf("\t       - steel,  apply random unimodular row ops to random diagonal matrix\n");
  printf("\t-d, --print_diag\n");
  printf("\t      print diagonal entries in HNF\n");
  printf("\t-v, --verify\n");
  printf("\t      run Monte Carlo verification of HNF\n");
  printf("\t-o, --out_filename FILE\n");
  printf("\t      write output to file in matrix market format\n");
  printf("\t-i, --in_filename FILE\n");
  printf("\t      read input from file in matrix market format; supersedes -n, -t, -l\n");
}

struct args {
  #define ARG(a, field, type, d, e, f) type field;
  ALL_ARGS
  #undef ARG
};
typedef struct args args_t;

static void set_empty_args(args_t * args)
{
  #define ARG(a, field, c, dflt_val, e, f) args->field = 0;
  ALL_ARGS
  #undef ARG
}
static void set_default_args(args_t * args)
{
  #define ARG(a, field, c, dflt_val, e, f) if(!args->field) {args->field = dflt_val;}
  ALL_ARGS
  #undef ARG
}

/* simple argument parsing callbacks */
static long do_num_arg(char const * arg)
{
  char * endptr = NULL;
  return strtol(arg, &endptr, 10);
}
static int do_flag_arg(char const * arg) { (void)arg; return 1; }
static char const * do_str_arg(char const * arg) { return arg; }


/* Fill arg structure from command line. */
static int parse_args(int argc, char ** argv, args_t * args)
{
  int c, indexptr;
  int rc = 1;

  #define ARG(sname, lname, c, d, flag, f) {#lname, flag, 0, sname}, 
  struct option longopts[] = { ALL_ARGS {0,0,0,0}};
  #undef ARG

  char shortops[100] = {0};
  #define ARG(name, b, c, d, flag, f) {\
    char buf[4] = {0}; \
    sprintf(buf, "%c%s", name, (flag==required_argument) ? ":" : ""); \
    strcat(shortops, buf); \
    }
  ALL_ARGS
  #undef ARG

  for(c = 0; c < argc; ++c) {
    fprintf(stderr, "%s ", argv[c]);
  }
  fprintf(stderr,"\n");


  /* fill arg structure with help from getopt */
  set_empty_args(args);
  while ((c = getopt_long(argc, argv, shortops, longopts, &indexptr)) != -1) {
    switch (c) {
      #define ARG(sname, lname, c, d, e, func) case sname: args->lname = func(optarg); break;
      ALL_ARGS
      #undef ARG
      default:
        assert(0);
        break;
    }
  }

  if (!args->ncols && !args->nrows && !args->in_filename) { rc = 0; }
  if (args->ncols && !args->nrows) { args->nrows = args->ncols; }
  if (args->nrows && !args->ncols) { args->ncols = args->nrows; }

  if (args->ncols != args->nrows) {
    if (args->type && 0 != strcmp(args->type, "random")) {
      rc = 0;
    }
  }

  if (args->rank) {
    if(args->type && 0 != strcmp(args->type, "rect")) {
      rc = 0;
    } else if (!args->type) {
      args->type = "rect";
    }
  }
  if (args->type && 0 == strcmp(args->type, "rect") && !args->rank) { rc = 0; }
  if(args->rank || (args->ncols != args->nrows) || args->in_filename) { args->rect = 1; }

  set_default_args(args);
  return rc;
}

static void print_diag(mpzMatrix_t const * A)
{
  long i;
  long n = A->nrows;
  for (i = 0; i < n; ++i) {
    mpz_out_str(stdout, 10, A->data[i*n + i]);
    printf(", ");
  }
  printf("\n");
}

int main(int argc, char ** argv)
{
  mpzMatrix_t *A, *H;
  int rc;

  args_t args;

  if (!parse_args(argc, argv, &args)) {
    print_usage();
    return 0;
  }
  g_print_stream = stderr;
  g_timer_stream = stderr;
  g_timer_threshold = args.timer_threshold;

  if (args.seed) {
    srand(args.seed);
  } else {
    srand(time(NULL));
  }

  if (args.in_filename) {
    A = mpzMatrix_initFromFile(args.in_filename);
    if (!A) {
      print_usage();
      printf("Nonexistent or invalid input file \"%s\".\n", args.in_filename);
      return 1;
    }
  } else {
    /* generate input as specified */
    A = iherm_input(args.type, args.nrows, args.ncols, args.rank, args.bitlen);

    if (!A) {
      print_usage();
      printf("Invalid input type \"%s\".\n", args.type);
      return 1;
    }
  }

  if(!args.rect) {
    TIMER("TOTAL", H = hermite(A);)
  } else {
    TIMER("TOTAL", H = hermiteRect(A);)
  }

  if (args.verify) {
    if(args.rect){
    mpzMatrix_writeToFile("A", A);
    mpzMatrix_writeToFile("H", H);
    } else{
    /* Monte Carlo verification*/

    rc = iherm_check(A, H);
    printf("%s\n", rc ? "Success." : "Fail!");}
  }

  if (args.out_filename) {
    mpzMatrix_writeToFile(args.out_filename, H);
  }

  if (args.print_diag) {
    print_diag(H);
  }

  mpzMatrix_fini(A);
  mpzMatrix_fini(H);
  return 0;
}

