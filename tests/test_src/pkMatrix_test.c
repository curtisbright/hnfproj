#include <stdio.h>
#include <stdlib.h>

#include "AllTests.h"
#include "CuTest.h"

#include "mpz_matrix.h"
#include "pk_matrix.h"
#include "timer.h"


#include "pk_matrix.c"

#define AssertMpzMatEquals(tc, ex, ac)  AssertMpzMatEquals_LineMsg(tc, __FILE__, __LINE__, ex, ac)

static void AssertMpzMatEquals_LineMsg(CuTest * tc, char const * file, int line, mpzMatrix_t const * expected, mpzMatrix_t const * actual)
{
  long i;
  char msg[1024];
  if (expected->ncols != actual->ncols || expected->nrows != actual->nrows) {
    snprintf(msg, 1024, "dimensions differ: expected %ldx%ld; actual: %ldx%ld\n", expected->nrows, expected->ncols, actual->nrows, actual->ncols);
    CuFail_Line(tc, file, line, NULL, msg);
    return;
  }
  for (i = 0; i < expected->ncols*expected->nrows; ++i) {
    if (0 != mpz_cmp(expected->data[i], actual->data[i])) {
      snprintf(msg, 1024, "entries differ at position %ld; expected %ld actual %ld", i, mpz_get_ui(expected->data[i]), mpz_get_ui(actual->data[i]));
      CuFail_Line(tc, file, line, NULL, msg);
      return;
    }
  }
  return;
}

static void test_toFull(CuTest *tc)
{
  long lA[] = {1, 2, 0, 5, 0,
               0, 3, 0, 2, 0,
               0, 0, 1, 5, 0,
               0, 0, 0, 6, 0,
               0, 0, 0, 0, 1};

  mpzMatrix_t * A = mpzMatrix_initSet(5, 5, lA);
  pk_t * Z = pkMatrix_fromFull(A);
  mpzMatrix_t * B = pkMatrix_toFull(Z);

  AssertMpzMatEquals(tc, A, B);
}

static void test_getNonTrivialIdxs(CuTest *tc)
{
  long lA[] = {1, 2, 0, 5, 0,
               0, 3, 0, 2, 0,
               0, 0, 1, 5, 0,
               0, 0, 0, 6, 0,
               0, 0, 0, 0, 1};
  mpzMatrix_t * A = mpzMatrix_initSet(5, 5, lA);
  long k_ret;
  long * idxs = getNonTrivialIdxs(A, &k_ret);

  CuAssertIntEquals(tc, 2, k_ret);
  CuAssertIntEquals(tc, 1, idxs[0]);
  CuAssertIntEquals(tc, 3, idxs[1]);
}

static void test_makeRevIdxs(CuTest *tc)
{
  long idxs[] = {2, 4, 5, 8};
  long * revIdxs = makeRevIdxs(10, 4, idxs);

  CuAssertIntEquals(tc, 0, revIdxs[2]);
  CuAssertIntEquals(tc, 1, revIdxs[4]);
  CuAssertIntEquals(tc, 2, revIdxs[5]);
  CuAssertIntEquals(tc, 3, revIdxs[8]);
}

static void test_makeMinIdxs(CuTest *tc)
{
  long idxs[] = {2, 4, 5, 8};
  long * minIdxs = makeMinIdxs(10, 4, idxs);

  CuAssertIntEquals(tc, 0, minIdxs[0]);
  CuAssertIntEquals(tc, 0, minIdxs[1]);
  CuAssertIntEquals(tc, 0, minIdxs[2]);
  CuAssertIntEquals(tc, 1, minIdxs[3]);
  CuAssertIntEquals(tc, 1, minIdxs[4]);
  CuAssertIntEquals(tc, 2, minIdxs[5]);
  CuAssertIntEquals(tc, 3, minIdxs[6]);
  CuAssertIntEquals(tc, 3, minIdxs[7]);
  CuAssertIntEquals(tc, 3, minIdxs[8]);
  CuAssertIntEquals(tc, 4, minIdxs[9]);
}




static void test_gemm_block(CuTest *tc)
{
  mpzMatrix_t *A, *B, *C, *AB;
  pk_t *S, *T, *Z;
  long lA[] = {1, 2, 0, 5, 0,
               0, 3, 0, 2, 0,
               0, 0, 1, 5, 0,
               0, 0, 0, 6, 0,
               0, 0, 0, 0, 1};

  long lB[] = {1, 0, 0, 5, 6,
               0, 1, 0, 6, 7,
               0, 0, 1, 5, 5,
               0, 0, 0, 8, 1,
               0, 0, 0, 0, 9};

  long lC[] = {1, 2, 0, 57, 25,
               0, 3, 0, 34, 23,
               0, 0, 1, 45, 10,
               0, 0, 0, 48, 6,
               0, 0, 0, 0,  9};

  A = mpzMatrix_initSet(5, 5, lA);
  B = mpzMatrix_initSet(5, 5, lB);
  C = mpzMatrix_initSet(5, 5, lC);

  S = pkMatrix_fromFull(A);
  T = pkMatrix_fromFull(B);
  Z = pkMatrix_init(A->nrows);

  pkMatrix_gemm_block(Z, S, T);

  AB = pkMatrix_toFull(Z);
  AssertMpzMatEquals(tc, C, AB);
}

static void test_gemm_iter(CuTest *tc)
{
  mpzMatrix_t *A, *B, *C, *AB;
  pk_t *S, *T, *Z;
  long lA[] = {1, 2, 0, 5, 0,
               0, 3, 0, 2, 0,
               0, 0, 1, 5, 0,
               0, 0, 0, 6, 0,
               0, 0, 0, 0, 1};

  long lB[] = {1, 0, 0, 5, 6,
               0, 1, 0, 6, 7,
               0, 0, 1, 5, 5,
               0, 0, 0, 8, 1,
               0, 0, 0, 0, 9};

  long lC[] = {1, 2, 0, 57, 25,
               0, 3, 0, 34, 23,
               0, 0, 1, 45, 10,
               0, 0, 0, 48, 6,
               0, 0, 0, 0,  9};

  A = mpzMatrix_initSet(5, 5, lA);
  B = mpzMatrix_initSet(5, 5, lB);
  C = mpzMatrix_initSet(5, 5, lC);

  S = pkMatrix_fromFull(A);
  T = pkMatrix_fromFull(B);
  Z = pkMatrix_init(A->nrows);

  pkMatrix_gemm_iter(Z, S, T);

  AB = pkMatrix_toFull(Z);
  AssertMpzMatEquals(tc, C, AB);
}

static void test_hnf(CuTest *tc)
{
  mpzMatrix_t *A, *H, *ZH;
  pk_t * Z;
  long lA[] = {1, 19, 0, 57, 25,
               0, 3,  0, 34, 23,
               0, 0,  1, 45, 10,
               0, 0,  0, 19, 6,
               0, 0,  0, 0,  5};

  long lH[] = {1, 1, 0, 5,  0,
               0, 3, 0, 15, 2,
               0, 0, 1, 7,  3,
               0, 0, 0, 19, 1,
               0, 0, 0, 0,  5};

  A = mpzMatrix_initSet(5, 5, lA);
  H = mpzMatrix_initSet(5, 5, lH);

  Z = pkMatrix_fromFull(A);
  pkMatrix_hnf(Z);
  ZH = pkMatrix_toFull(Z);

  AssertMpzMatEquals(tc, H, ZH);
}


static void test_applyVector_simple(CuTest *tc)
{
  mpzMatrix_t *A, *x, *b;
  pk_t * Z;
  long lA[] = {1, 19, 0, 57, 25,
               0, 3,  0, 34, 23,
               0, 0,  1, 45, 10,
               0, 0,  0, 19, 6,
               0, 0,  0, 0,  5};
  long lx[] = {12,  29, 92, 43,
               51,  44, 68, 31,
               23,  92, 54, 96,
               54, -31, 78, 61,
               49,  67, 60, 25};
  long lb[] = {12,  773, 7330, 43,
               51,  619, 4236, 31,
               23, -633, 4164, 96,
               54, -187, 1842, 61,
               49,  335,  300, 25};

  A = mpzMatrix_initSet(5, 5, lA);
  x = mpzMatrix_initSet(5, 4, lx);
  b = mpzMatrix_initSet(5, 4, lb);

  Z = pkMatrix_fromFull(A);
  pkMatrix_applyVector_simple(x, Z, 1, 2);

  AssertMpzMatEquals(tc, b, x);
}

static void test_applyVector_block(CuTest *tc)
{
  mpzMatrix_t *A, *x, *b;
  pk_t * Z;
  long lA[] = {1, 19, 0, 57, 25,
               0, 3,  0, 34, 23,
               0, 0,  1, 45, 10,
               0, 0,  0, 19, 6,
               0, 0,  0, 0,  5};
  long lx[] = {12,  29, 92, 43,
               51,  44, 68, 31,
               23,  92, 54, 96,
               54, -31, 78, 61,
               49,  67, 60, 25};
  long lb[] = {12,  773, 7330, 43,
               51,  619, 4236, 31,
               23, -633, 4164, 96,
               54, -187, 1842, 61,
               49,  335,  300, 25};

  A = mpzMatrix_initSet(5, 5, lA);
  x = mpzMatrix_initSet(5, 4, lx);
  b = mpzMatrix_initSet(5, 4, lb);

  Z = pkMatrix_fromFull(A);
  pkMatrix_applyVector_block(x, Z, 1, 2);

  AssertMpzMatEquals(tc, b, x);
}

static void test_applyHinv(CuTest * tc)
{
  long lA[] = { -81,   -98,   -76,   -4,    29,
                -38,   -77,   -72,   27,    44,
                -18,    57,    -2,    8,    92,
                 87,    27,   -32,   69,   -31,
                 33,   -93,   -74,   99,    67};
  long lH[] = {  1, 19, 0, 57, 25,
                 0, 3,  0, 34, 23,
                 0, 0,  1, 45, 10,
                 0, 0,  0, 19, 6,
                 0, 0,  0, 0,  5};
  long lAH[] = { -81,    480,   -76,   -436,   -1122,
                 -38,    215,   -72,    -98,    -528,
                 -18,    133,    -2,   -178,    -285,
                  87,   -542,   -32,    788,    1170,
                  33,   -240,   -74,    510,     488};


  mpzMatrix_t *A, *Hfull, *AH, *AH_chk;
  pk_t * H;

  AH = mpzMatrix_init(5, 5);
  A = mpzMatrix_initSet(5, 5, lA);
  Hfull = mpzMatrix_initSet(5, 5, lH);
  AH_chk = mpzMatrix_initSet(5, 5, lAH);
  H = pkMatrix_fromFull(Hfull);

  pkMatrix_applyHinv(AH, A, H);

  AssertMpzMatEquals(tc, AH, AH_chk);
}

CuSuite * pkMatrix_suite(void)
{
  CuSuite * suite = CuSuiteNew();

  SUITE_ADD_TEST(suite, test_toFull);
  SUITE_ADD_TEST(suite, test_getNonTrivialIdxs);
  SUITE_ADD_TEST(suite, test_makeRevIdxs);
  SUITE_ADD_TEST(suite, test_makeMinIdxs);

  SUITE_ADD_TEST(suite, test_hnf);
  SUITE_ADD_TEST(suite, test_gemm_iter);
  SUITE_ADD_TEST(suite, test_gemm_block);
  SUITE_ADD_TEST(suite, test_applyHinv);

  SUITE_ADD_TEST(suite, test_applyVector_simple);
  SUITE_ADD_TEST(suite, test_applyVector_block);

  return suite;
}
