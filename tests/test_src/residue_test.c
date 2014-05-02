#include <stdio.h>
#include <stdlib.h>

#include "AllTests.h"
#include "CuTest.h"

#include "residue.h"

#define AssertResidueEquals(tc, ex, ac)  AssertResidueEquals_LineMsg(tc, __FILE__, __LINE__, ex, ac)

static void AssertResidueEquals_LineMsg(CuTest * tc, char const * file, int line, residue_t const * expected, residue_t const * actual)
{
  long i;
  char msg[1024];
  if (expected->ncols != actual->ncols || expected->nrows != actual->nrows) {
    snprintf(msg, 1024, "dimensions differ: expected %ldx%ld; actual: %ldx%ld\n", expected->nrows, expected->ncols, actual->nrows, actual->ncols);
    CuFail_Line(tc, file, line, NULL, msg);
    return;
  }
  for (i = 0; i < expected->ncols*expected->nrows; ++i) {
    long exp = residue_getEntry(expected, i);
    long act = residue_getEntry(actual, i);
    if (exp != act) {
      snprintf(msg, 1024, "entries differ at position %ld; expected %ld actual %ld", i, (exp), (act));
      CuFail_Line(tc, file, line, NULL, msg);
      return;
    }
  }
  return;
}

static void set_residue(residue_t * A, long * lA)
{
  long i;
  for (i = 0; i < A->ncols*A->nrows; ++i) {residue_setEntry(A, i, lA[i]);}
}

static long lA[] = {3, 2, 7, -7, -9, 7};
static long lB[] = {-1, -4, 5, -3, -6, -6};
static long lC[] = {1, 2, 3, 4, 5, 6, 7, 8};
static residue_t * A;
static residue_t * B;
static residue_t * C;

static void init(void)
{
  A = residue_init(3, 2, 13);
  B = residue_init(3, 2, 13);
  C = residue_init(2, 4, 13);
  set_residue(A, lA);
  set_residue(B, lB);
  set_residue(C, lC);
}

static void fini(void)
{
  residue_fini(A);
  residue_fini(B);
  residue_fini(C);
}

static void test_copy(CuTest* tc)
{
  residue_t * R = residue_init(3, 2, 13);
  init();

  atlas_copy(B, A);
  atlas_copy(R, B);

  AssertResidueEquals(tc, A, B);
  AssertResidueEquals(tc, R, B);

  residue_zero(A);
  AssertResidueEquals(tc, R, B);

  fini();
}

static void test_identity(CuTest* tc)
{
  long lR[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  residue_t * M = residue_init(3, 3, 13);
  residue_t * R = residue_init(3, 3, 13);
  init();

  set_residue(R, lR);
  atlas_identity(M);
  AssertResidueEquals(tc, R, M);

  fini();
}

static void test_isZero(CuTest * tc)
{
  residue_t * M = residue_init(4, 4, 13);

  residue_zero(M);
  CuAssertTrue(tc, residue_isZero(M));
  residue_identity(M);
  CuAssertTrue(tc, !residue_isZero(M));
}


static void test_scale(CuTest* tc)
{
  long lR1[] = {-1,-5,2,-2,3,2};
  residue_t * R1 = residue_init(3, 2, 13);
  init();
  set_residue(R1, lR1);

  residue_scale(A, 4);
  AssertResidueEquals(tc, R1, A);

  fini();
}

static void test_add(CuTest* tc)
{
  residue_t *R1, *R2;
  long lR1[] = {2, -2, -1, 3, -2, 1};
  long lR2[] = {6, 4, 1, -1, -5, 1};
  init();
  R1 = residue_init(A->nrows, A->ncols, 13);
  R2 = residue_init(A->nrows, A->ncols, 13);
  set_residue(R1, lR1);
  set_residue(R2, lR2);

  residue_add(B, A, 1);
  AssertResidueEquals(tc, R1, B);

  residue_add(A, A, 1);
  AssertResidueEquals(tc, R2, A);

  fini();
}

static void test_gemm(CuTest* tc)
{
  long lR1[] = {5, 0, -5, 3, 3, 5, -6, -4, 3, 4, 5, 6};
  residue_t * R1 = residue_init(3, 4, 13);
  residue_t * D = residue_init(3, 4, 13);
  set_residue(R1, lR1);
  init();

  /* D = B.C */
  residue_gemm(D, B, C);
  AssertResidueEquals(tc, R1, D);

  fini();
}

static void test_mods(CuTest* tc)
{
  long i;
  const long p = 101;
  const long ncols = 40;
  const long nrows = 50;
  long lM[2000];
  long lR[2000];
  long RR[101];
  residue_t * M = residue_init(nrows, ncols, p);
  residue_t * R = residue_init(nrows, ncols, p);
  init();

  for (i = 0; i <= p/2; ++i) {
    RR[i] = i;
  }
  for (i = p/2+1; i < p; ++i) {
    RR[i] = i-p;
  }

  for (i = 0; i < ncols*nrows; ++i) {
    lM[i] = i- p*((ncols*nrows/2)/p);
    lR[i] = RR[i%p];
  }
  set_residue(M, lM);
  set_residue(R, lR);

  residue_mods(M);

  AssertResidueEquals(tc, R, M);

  fini();
}


CuSuite * residue_suite(void)
{
  CuSuite * suite = CuSuiteNew();

  SUITE_ADD_TEST(suite, test_copy);
  SUITE_ADD_TEST(suite, test_identity);
  SUITE_ADD_TEST(suite, test_isZero);
  SUITE_ADD_TEST(suite, test_scale);
  SUITE_ADD_TEST(suite, test_add);
  SUITE_ADD_TEST(suite, test_gemm);
  SUITE_ADD_TEST(suite, test_mods);

  return suite;
}
