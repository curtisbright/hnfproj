#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "gmp.h"

#include "iherm.h"
#include "iherm_utils.h"
#include "timer.h"

int g_timer_depth = 0;
int g_print_level = 2;
FILE * g_print_stream;

static int iherm_test(long n, long type, long k1, long k2)
{
  mpzMatrix_t *A, *H;
  int rc;
  char lbl[100] = {0};
  static clock_t t0, t1;
  static unsigned long long dt;

  snprintf(lbl, 100, "%s %ld %ld", "TOTAL", n, type);

  A = iherm_input(n, type);

  t0 = clock();
  H = myHermite(A, k1, k2);
  t1 = clock();
  dt = (unsigned long long)(t1 - t0) * 1000 / (unsigned long long)CLOCKS_PER_SEC;
  printf("%16s %llu.%.3llu\n", lbl, dt / 1000, dt % 1000);
  fflush(stdout);

  rc = iherm_check(A, H);
  if (!rc) { printf("Fail!"); }

  mpzMatrix_fini(A);
  mpzMatrix_fini(H);

  return 0;
}

int main(void)
{
  int i, j=1;
  g_print_stream = stderr;
  srand(0);

  iherm_test(100, 1, 0, 0);
  iherm_test(200, 1, 0, 0);
  iherm_test(400, 1, 0, 0);
  iherm_test(100, 6, 0, 0);
  iherm_test(200, 6, 0, 0);
  iherm_test(400, 6, 0, 0);
  /*for (j = 1; j < 50; ++j) {
    for (i = j/10+1; i < 400/10; ++i) {
      iherm_test(400, 6, j, 10*i);
    }
  }
  for (j = 1; j < 50; ++j) {
    for (i = j/10+1; i < 400/10; ++i) {
      iherm_test(400, 6, j, 10*i);
    }
  }*/


 


  return 0;
}

