#include <stdio.h>
#include <stdlib.h>

#include "AllTests.h"
#include "CuTest.h"

#include "mpz_matrix.h"
#include "pk_matrix.h"
#include "timer.h"

#include "iherm.c"
#include "hnfproj.c"

static void test_hcol(CuTest * tc)
{
  mpz_t d;
  long lw[] = {8, 104, 72, 20};
  mpzMatrix_t * w = mpzMatrix_initSet(4, 1, lw);


  long lH[] = { 1, 0,  1, 17,
                0, 1,  3,  5,
                0, 0,  5,  3,
                0, 0,  0, 21 };
  mpzMatrix_t * Hfull_chk = mpzMatrix_initSet(4, 4, lH);
  mpzMatrix_t * Hfull;
  pk_t * H = pkMatrix_init(4);

  mpz_init_set_ui(d, 105);

  hcol(H, w, d);

  Hfull = pkMatrix_toFull(H);

  CuAssertTrue(tc, mpzMatrix_equal(Hfull, Hfull_chk));
}

static void test_hermite(CuTest * tc)
{
  long lA[] = { 115, -85,  515, 1209,
                -57,  40, -242, -532,
                 33, -29,  176,  491,
                -54,  44, -267, -695};
  long lH[] = { 1, 0,  1, 17,
                0, 1,  3,  5,
                0, 0,  5,  3,
                0, 0,  0, 21 };

  mpzMatrix_t * A = mpzMatrix_initSet(4, 4, lA);
  mpzMatrix_t * H_chk = mpzMatrix_initSet(4, 4, lH);
  mpzMatrix_t * H;

  H = hermite(A);

  CuAssertTrue(tc, mpzMatrix_equal(H, H_chk));
}

static void test_hermiteLong(CuTest * tc)
{
  long i;
  long A[] = { 115, -85,  515, 1209,
                -57,  40, -242, -532,
                 33, -29,  176,  491,
                -54,  44, -267, -695};
  long lH[] = { 1, 0,  1, 17,
                0, 1,  3,  5,
                0, 0,  5,  3,
                0, 0,  0, 21 };
  mpz_t * H;

  H = hermiteLong(A, 4);

  for (i = 0; i < 4*4; ++i) {
    CuAssertIntEquals(tc, lH[i], mpz_get_si(H[i]));
  }
}

static void test_hermiteMpz(CuTest * tc)
{
  long i;
  long lA[] = { 115, -85,  515, 1209,
                -57,  40, -242, -532,
                 33, -29,  176,  491,
                -54,  44, -267, -695};
  long lH[] = { 1, 0,  1, 17,
                0, 1,  3,  5,
                0, 0,  5,  3,
                0, 0,  0, 21 };
  mpz_t * H;
  mpzMatrix_t * A = mpzMatrix_initSet(4, 4, lA);

  H = hermiteMpz((mpz_t const *)A->data, 4);

  for (i = 0; i < 4*4; ++i) {
    CuAssertIntEquals(tc, lH[i], mpz_get_si(H[i]));
  }
}

CuSuite * iherm_suite(void)
{
  CuSuite * suite = CuSuiteNew();

  SUITE_ADD_TEST(suite, test_hcol);
  SUITE_ADD_TEST(suite, test_hermite);
  SUITE_ADD_TEST(suite, test_hermiteLong);
  SUITE_ADD_TEST(suite, test_hermiteMpz);

  return suite;
}
