
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "AllTests.h"
#include "CuTest.h"

#include "arith_utils.h"
#include "basis.h"
#include "mpz_matrix.h"
#include "reconstruct.h"
#include "rns_conversion.h"
#include "rns_matrix.h"

#include "basis.c"

static void test_initModulus(CuTest *tc)
{
  double delta = 0;
  long p = prevprime(pickStartModulus(20000));
  modulus_t M;
  modulus_init(&M, p);

  CuAssertIntEquals(tc, p, M.p);
  CuAssertIntEquals(tc, M.p, (long)M.pd);
  CuAssertDblEquals(tc, M.pd, (double)M.p, delta);
  CuAssertDblEquals(tc, 1.0, M.p*M.inv_p, delta);
}

static void test_calcBezout(CuTest *tc)
{
  mpz_t total;
  mpz_t * L;
  long * S;
  long primes[] = {11, 13, 23};
  int nprimes = 3;

  mpz_init_set_ui(total, 11*13*23);
  L = calcLagrange(total, primes, nprimes);
  S = calcBezout((mpz_t const *)L, primes, nprimes);

  CuAssertIntEquals(tc, -5, S[0]);
  CuAssertIntEquals(tc, -2, S[1]);
  CuAssertIntEquals(tc, -9, S[2]);

}

static void test_calcLagrange(CuTest *tc)
{
  mpz_t total;
  mpz_t * L;
  long primes[] = {11, 13, 23};
  int nprimes = 3;

  mpz_init_set_ui(total, 11*13*23);
  L = calcLagrange(total, primes, nprimes);

  CuAssertIntEquals(tc, 299, mpz_get_ui(L[0]));
  CuAssertIntEquals(tc, 253, mpz_get_ui(L[1]));
  CuAssertIntEquals(tc, 143, mpz_get_ui(L[2]));
}

static void test_calcWeights(CuTest *tc)
{
  mpz_t total;
  mpz_t * L;
  long * S;
  double * W;
  long primes[] = {3, 5, 7};
  int nprimes = 3;
  double delta = 0;

  mpz_init_set_ui(total, 3*5*7);
  L = calcLagrange(total, primes, nprimes);
  S = calcBezout((mpz_t const *)L, primes, nprimes);
  W = calcWeights(S, primes, nprimes);

  CuAssertDblEquals(tc, -1.0/3, W[0], delta);
  CuAssertDblEquals(tc, 1.0/5, W[1], delta);
  CuAssertDblEquals(tc, 1.0/7, W[2], delta);
}

static void test_calcRecon(CuTest *tc)
{
  mpz_t total;
  mpz_t * L;
  long * S;
  mpz_t * C;
  long primes[] = {11, 13, 23};
  int nprimes = 3;

  mpz_init_set_ui(total, 11*13*23);
  L = calcLagrange(total, primes, nprimes);
  S = calcBezout((mpz_t const *)L, primes, nprimes);
  C = calcRecon((mpz_t const *)L, S, nprimes);

  CuAssertIntEquals(tc, -1495, mpz_get_si(C[0]));
  CuAssertIntEquals(tc, -506, mpz_get_si(C[1]));
  CuAssertIntEquals(tc, -1287, mpz_get_si(C[2]));

}

/*static void test_calcMixedInv(CuTest * tc)
{
  long primes[] = {11, 13, 23, 31};
  int nprimes = 4;
  double * mixed_inv = calcMixedInv(primes, nprimes);
  double delta = 0;

  CuAssertDblEquals(tc, -1.0, mixed_inv[0], delta);
  CuAssertDblEquals(tc, -6.0, mixed_inv[1], delta);
  CuAssertDblEquals(tc, 9.0, mixed_inv[2], delta);
  CuAssertDblEquals(tc, 10.0, mixed_inv[3], delta);
}*/

static void test_pickStartModulus(CuTest * tc)
{
  long n = 10000;
  long p = pickStartModulus(n);
  long q = (1L << 26) * sqrt(2.0) / sqrt((double)n); /*approx bound*/

  CuAssertTrue(tc, p >= q);
  CuAssertTrue(tc, p - q < 10); /*approx bound should be close */
  CuAssertTrue(tc, checkStartModulus(n, q));
  CuAssertTrue(tc, checkStartModulus(n, p));
}

static void test_genPrimes(CuTest * tc)
{
  int n = 0;
  long * p;
  mpz_t bound;
  mpz_init_set_ui(bound, 11*13*17-1);

  p = genPrimes(18, bound, &n);

  CuAssertIntEquals(tc, 3, n);
  CuAssertIntEquals(tc, 17, p[0]);
  CuAssertIntEquals(tc, 13, p[1]);
  CuAssertIntEquals(tc, 11, p[2]);
}

static void test_initBasis(CuTest * tc)
{
  basis_t * P = 0;
  long primes[] = {11, 13, 23, 31};
  int nprimes = 4;

  P = basis_init(primes, nprimes);

  CuAssertPtrNotNull(tc, P);
  CuAssertIntEquals(tc, nprimes, P->nmod);
  CuAssertIntEquals(tc, 11*13*23*31, mpz_get_ui(P->PP));
  CuAssertPtrNotNull(tc, P->mod);
  CuAssertPtrNotNull(tc, P->L);
  CuAssertPtrNotNull(tc, P->inv_L);
  CuAssertPtrNotNull(tc, P->W);
}


static void test_reconstruct(CuTest * tc)
{
  rnsMatrix_t * Ap;
  mpzMatrix_t * A;
  long lP[4] = {11, 13, 17, 19};
  basis_t * P = basis_init(lP, 4);

  Ap = rnsMatrix_init(1, 1, P);
  A = mpzMatrix_init(1, 1);

  rnsMatrix_setEntry(Ap, 0, 0, 5);
  rnsMatrix_setEntry(Ap, 1, 0, 4);
  rnsMatrix_setEntry(Ap, 2, 0, 5);
  rnsMatrix_setEntry(Ap, 3, 0, 1);

  mpzMatrix_reconstruct(A, Ap);

  CuAssertIntEquals(tc, 20388, mpz_get_si(A->data[0]));
}

static void test_convResidue_small(CuTest * tc)
{
  long n;
  rnsMatrix_t *Ap, *Aq;
  long lP[] = {11, 13, 17};
  long lQ[] = {23, 29, 31};

  basis_t * P = basis_init(lP, 3);
  basis_t * Q = basis_init(lQ, 3);

  Ap = rnsMatrix_init(1, 1, P);
  Aq = rnsMatrix_init(1, 1, Q);

  n = 200;
  rnsMatrix_setEntry(Ap, 0, 0, modsL(n, lP[0]));
  rnsMatrix_setEntry(Ap, 1, 0, modsL(n, lP[1]));
  rnsMatrix_setEntry(Ap, 2, 0, modsL(n, lP[2]));

  rnsMatrix_convert(Aq, Ap);

  CuAssertDblEquals(tc, modsL(n, lQ[0]), rnsMatrix_getEntry(Aq, 0, 0), 0);
  CuAssertDblEquals(tc, modsL(n, lQ[1]), rnsMatrix_getEntry(Aq, 1, 0), 0);
  CuAssertDblEquals(tc, modsL(n, lQ[2]), rnsMatrix_getEntry(Aq, 2, 0), 0);
}


static void test_convResidue(CuTest * tc)
{
  long n, nn;
  rnsMatrix_t *Ap1, *Ap2, *Aq;
  long lP[] = {11, 13, 17};
  long lQ[] = {23, 29};
  basis_t * P = basis_init(lP, 3);
  basis_t * Q = basis_init(lQ, 2);

  Ap1 = rnsMatrix_init(1, 1, P);
  Ap2 = rnsMatrix_init(1, 1, P);
  Aq = rnsMatrix_init(1, 1, Q);

  for (nn = 0;   nn < lQ[0]*lQ[1]*2; ++nn) {
    n = modsL(nn, lQ[0]*lQ[1]);
    rnsMatrix_setEntry(Aq, 0, 0, modsL(n, lQ[0]));
    rnsMatrix_setEntry(Aq, 1, 0, modsL(n, lQ[1]));

    rnsMatrix_setEntry(Ap1, 0, 0, modsL(n, lP[0]));
    rnsMatrix_setEntry(Ap1, 1, 0, modsL(n, lP[1]));
    rnsMatrix_setEntry(Ap1, 2, 0, modsL(n, lP[2]));

    rnsMatrix_convert(Ap2, Aq);


    CuAssertDblEquals(tc, rnsMatrix_getEntry(Ap1, 0, 0), rnsMatrix_getEntry(Ap2, 0, 0), 0);
    CuAssertDblEquals(tc, rnsMatrix_getEntry(Ap1, 1, 0), rnsMatrix_getEntry(Ap2, 1, 0), 0);
    CuAssertDblEquals(tc, rnsMatrix_getEntry(Ap1, 2, 0), rnsMatrix_getEntry(Ap2, 2, 0), 0);
  }
}

static void test_convResidue2(CuTest * tc)
{
  long n, nn;
  rnsMatrix_t *Ap1, *Ap2, *Aq;
  long lP[] = {23, 29};
  long lQ[] = {11, 13, 17};
  basis_t * P = basis_init(lP, 2);
  basis_t * Q = basis_init(lQ, 3);

  Ap1 = rnsMatrix_init(1, 1, P);
  Ap2 = rnsMatrix_init(1, 1, P);
  Aq = rnsMatrix_init(1, 1, Q);

  for (nn = 0; nn < lQ[0]*lQ[1]*lQ[2]*2; ++nn) {
    n = modsL(nn, lQ[0]*lQ[1]*lQ[2]);
    rnsMatrix_setEntry(Aq, 0, 0, modsL(n, lQ[0]));
    rnsMatrix_setEntry(Aq, 1, 0, modsL(n, lQ[1]));
    rnsMatrix_setEntry(Aq, 2, 0, modsL(n, lQ[2]));

    rnsMatrix_setEntry(Ap1, 0, 0, modsL(n, lP[0]));
    rnsMatrix_setEntry(Ap1, 1, 0, modsL(n, lP[1]));

    rnsMatrix_convert(Ap2, Aq);

    CuAssertDblEquals(tc, rnsMatrix_getEntry(Ap1, 0, 0), rnsMatrix_getEntry(Ap2, 0, 0), 0);
    CuAssertDblEquals(tc, rnsMatrix_getEntry(Ap1, 1, 0), rnsMatrix_getEntry(Ap2, 1, 0), 0);
  }

}

static void test_convResidue3(CuTest * tc)
{
  rnsMatrix_t *Ap1, *Ap2, *Aq;
  long lQ[] = {6710887, 6710897, 6710909};
  long lP[] = {6710927, 6710929, 6710947};
  basis_t * P = basis_init(lP, 3);
  basis_t * Q = basis_init(lQ, 3);

  Ap1 = rnsMatrix_init(1, 1, P);
  Ap2 = rnsMatrix_init(1, 1, P);
  Aq = rnsMatrix_init(1, 1, Q);

  rnsMatrix_setEntry(Aq, 0, 0, 2175439);
  rnsMatrix_setEntry(Aq, 1, 0, -3264792);
  rnsMatrix_setEntry(Aq, 2, 0, 944438);

  rnsMatrix_setEntry(Ap1, 0, 0, 558336);
  rnsMatrix_setEntry(Ap1, 1, 0, -3211407);
  rnsMatrix_setEntry(Ap1, 2, 0, 3146388);

  rnsMatrix_convert(Ap2, Aq);

  CuAssertDblEquals(tc, rnsMatrix_getEntry(Ap1, 0, 0), rnsMatrix_getEntry(Ap2, 0, 0), 0);
  CuAssertDblEquals(tc, rnsMatrix_getEntry(Ap1, 1, 0), rnsMatrix_getEntry(Ap2, 1, 0), 0);
  CuAssertDblEquals(tc, rnsMatrix_getEntry(Ap1, 2, 0), rnsMatrix_getEntry(Ap2, 2, 0), 0);
}

CuSuite * rns_suite(void)
{
  CuSuite * suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, test_initModulus);

  SUITE_ADD_TEST(suite, test_calcLagrange);
  SUITE_ADD_TEST(suite, test_calcBezout);
  /*SUITE_ADD_TEST(suite, test_calcMixedInv);*/
  SUITE_ADD_TEST(suite, test_calcWeights);
  SUITE_ADD_TEST(suite, test_calcRecon);
  SUITE_ADD_TEST(suite, test_pickStartModulus);

  SUITE_ADD_TEST(suite, test_initBasis);

  SUITE_ADD_TEST(suite, test_convResidue_small);
  SUITE_ADD_TEST(suite, test_convResidue);
  SUITE_ADD_TEST(suite, test_convResidue2);
  SUITE_ADD_TEST(suite, test_convResidue3);

  SUITE_ADD_TEST(suite, test_reconstruct);
  SUITE_ADD_TEST(suite, test_genPrimes);

  return suite;
}
