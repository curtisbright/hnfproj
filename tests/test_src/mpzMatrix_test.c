#include <stdio.h>
#include <stdlib.h>

#include "AllTests.h"
#include "CuTest.h"

#include "mpz_matrix.h"
#include "pk_matrix.h"
#include "timer.h"

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

static long lA[] = {3, 2, 7, -7, -9, 7};
static long lB[] = {-1, -4, 5, -3, -6, -6};
static long lC[] = {1, 2, 3, 4, 5, 6, 7, 8};
static mpzMatrix_t * A;
static mpzMatrix_t * B;
static mpzMatrix_t * C;

static void init(void)
{
  A = mpzMatrix_initSet(3, 2, lA);
  B = mpzMatrix_initSet(3, 2, lB);
  C = mpzMatrix_initSet(2, 4, lC);
}

static void fini(void)
{
  mpzMatrix_fini(A);
  mpzMatrix_fini(B);
  mpzMatrix_fini(C);
}


static void test_zero(CuTest* tc)
{
  mpzMatrix_t * R;
  init();
  R = mpzMatrix_init(A->nrows, A->ncols);

  mpzMatrix_zero(A);

  AssertMpzMatEquals(tc, R, A);

  fini();
}

static void test_set(CuTest* tc)
{
  mpzMatrix_t *R1, *R2;
  init();
  R1 = mpzMatrix_init(A->nrows, A->ncols);
  R2 = mpzMatrix_init(A->nrows, A->ncols);

  mpzMatrix_set(R1, A);
  mpzMatrix_set(R2, B);
  AssertMpzMatEquals(tc, R1, A);
  AssertMpzMatEquals(tc, R2, B);

  mpzMatrix_swap(R1, R2);

  AssertMpzMatEquals(tc, R1, B);
  AssertMpzMatEquals(tc, R2, A);

  fini();
}

static void test_swap(CuTest* tc)
{
  mpzMatrix_t * R;
  init();

  R = mpzMatrix_init(A->nrows, A->ncols);
  mpzMatrix_set(R, A);

  AssertMpzMatEquals(tc, R, A);

  fini();
}

static void test_scale(CuTest* tc)
{
  mpzMatrix_t * R1;
  mpz_t k;
  long lR1[] = {300, 200, 700, -700, -900, 700};

  init();

  R1 = mpzMatrix_initSet(A->nrows, A->ncols, lR1);
  mpz_init_set_ui(k, 100);

  mpzMatrix_scale(A, k);

  AssertMpzMatEquals(tc, R1, A);

  fini();
}

static void test_mods(CuTest* tc)
{
  long i;
  mpz_t k;
  const long n = 40;
  const long m = 50;
  long lM[2000];
  long lR[2000];
  long RR[] = {0, 1, 2, 3, -3, -2, -1};
  mpzMatrix_t * M;
  mpzMatrix_t * R;

  init();

  mpz_init_set_ui(k, 7);
  for (i = 0; i < n*m; ++i) {
    lM[i] = i- 7*((n*m/2)/7);
    lR[i] = RR[i%7];
  }
  M = mpzMatrix_initSet(m, n, lM);
  R = mpzMatrix_initSet(m, n, lR);

  mpzMatrix_mods(M, k);

  AssertMpzMatEquals(tc, R, M);

  fini();
}

static void test_add(CuTest* tc)
{
  mpzMatrix_t *R1, *R2;
  long lR1[] = {2, -2, 12, -10, -15, 1};
  long lR2[] = {6, 4, 14, -14, -18, 14};
  init();

  R1 = mpzMatrix_initSet(A->nrows, A->ncols, lR1);
  R2 = mpzMatrix_initSet(A->nrows, A->ncols, lR2);

  mpzMatrix_add(B, A);
  AssertMpzMatEquals(tc, R1, B);

  mpzMatrix_add(A, A);
  AssertMpzMatEquals(tc, R2, A);

  fini();
}

static void test_gemm(CuTest* tc)
{
  mpzMatrix_t *R1, *D;
  long lR1[] = {-21, -26, -31, -36, -10, -8, -6, -4, -36, -48, -60, -72};
  init();

  R1 = mpzMatrix_initSet(3, 4, lR1);
  D = mpzMatrix_init(3, 4);

  /* D = B.C */
  mpzMatrix_gemm(D, B, C);
  AssertMpzMatEquals(tc, R1, D);

  fini();
}

static void test_rnsGemm(CuTest* tc)
{
  mpzMatrix_t *R1, *D;
  long lR1[] = {-21, -26, -31, -36, -10, -8, -6, -4, -36, -48, -60, -72};
  init();

  R1 = mpzMatrix_initSet(3, 4, lR1);
  D = mpzMatrix_init(3, 4);

  /* D = B.C */
  mpzMatrix_rnsGemm(D, B, C);
  AssertMpzMatEquals(tc, R1, D);

  fini();
}

#if 0
static void test_gemm_alias(CuTest* tc)
{
  long lR1[] = {-21, -26, -31, -36, -10, -8, -6, -4, -36, -48, -60, -72};
  long lD[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  mpzMatrix_t * R1 = mpzMatrix_initSet(3, 3, lR1);
  mpzMatrix_t * D = mpzMatrix_initSet(3, 3, lD);
  init();

  /* D = B.C */
  mpzMatrix_print(D);
  mpzMatrix_gemm(D, D, D);
  mpzMatrix_print(D);
  AssertMpzMatEquals(tc, R1, D);

  fini();
}
#endif

static void test_max(CuTest* tc)
{
  mpzMatrix_t *M;
  mpz_t Mmax;
  long lM[] = {-21, -26, -31, -36, -10, -8, -6, -4, -36, -48, -60, -72};

  M = mpzMatrix_initSet(3, 4, lM);
  mpz_init(Mmax);

  mpzMatrix_max(Mmax, M);
  CuAssertIntEquals(tc, 72, mpz_get_si(Mmax));
}

static void test_isZero(CuTest* tc)
{
  mpzMatrix_t *M;
  long lM[] = {-21, -26, -31, -36, -10, -8, -6, -4, -36, -48, -60, -72};

  M = mpzMatrix_initSet(3, 4, lM);
  B = mpzMatrix_init(4, 4);

  CuAssertTrue(tc, !mpzMatrix_isZero(M));
  CuAssertTrue(tc, mpzMatrix_isZero(B));

  mpzMatrix_zero(M);
  CuAssertTrue(tc, mpzMatrix_isZero(M));
}

CuSuite * mpzMatrix_suite(void)
{
  CuSuite * suite = CuSuiteNew();

  SUITE_ADD_TEST(suite, test_zero);
  SUITE_ADD_TEST(suite, test_set);
  SUITE_ADD_TEST(suite, test_swap);

  SUITE_ADD_TEST(suite, test_scale);
  SUITE_ADD_TEST(suite, test_mods);
  SUITE_ADD_TEST(suite, test_add);
  SUITE_ADD_TEST(suite, test_gemm);
  SUITE_ADD_TEST(suite, test_rnsGemm);
  /*SUITE_ADD_TEST(suite, test_gemm_alias);*/

  SUITE_ADD_TEST(suite, test_max);
  SUITE_ADD_TEST(suite, test_isZero);

  return suite;
}
