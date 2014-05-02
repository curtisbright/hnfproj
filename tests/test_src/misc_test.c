#include <stdlib.h>

#include "AllTests.h"
#include "CuTest.h"

#include "arith_utils.h"
#include "basis.h"

static void test_nextprime(CuTest * tc)
{
  CuAssertIntEquals(tc, 97, nextprime(90));
  CuAssertIntEquals(tc, 123456791, nextprime(123456789));
}
static void test_prevprime(CuTest * tc)
{
  CuAssertIntEquals(tc, 89, prevprime(97));
  CuAssertIntEquals(tc, 123456761, prevprime(123456791));
}

static void test_modsMpz(CuTest * tc)
{
  mpz_t nn;
  int n;
  int p = 5;
  int r1[] = {0, 1, 2, -2, -1};
  int r2[] = {0, -1, -2, 2, 1};

  mpz_init(nn);

  for(n = 0; n <= 100; ++n) {
    mpz_set_si(nn, n);
    CuAssertIntEquals(tc, r1[n % p], modsMpzL(nn, p));
  }
  for(n = 0; n >= -100; --n) {
    mpz_set_si(nn, n);
    CuAssertIntEquals(tc, r2[(-n) % p], modsMpzL(nn, p));
  }
}

static void test_modsL(CuTest * tc)
{
  int n;
  int p = 5;
  int r1[] = {0, 1, 2, -2, -1};
  int r2[] = {0, -1, -2, 2, 1};
  for(n = 0; n <= 100; ++n) {
    CuAssertIntEquals(tc, r1[n % p], modsL(n, p));
  }
  for(n = 0; n >= -100; --n) {
    CuAssertIntEquals(tc, r2[(-n) % p], modsL(n, p));
  }
}

static void test_modsL1(CuTest * tc)
{
  long n = 1216;
  long p = 23;
  long r = modsL(n, p);
  CuAssertIntEquals(tc, -3, r);
}


#include "mods.c"
static void test_mods(CuTest* tc)
{
  long ip, ia, np, na;
  modulus_t M;

  np = 10;
  na = 10;

  for (ip = 0; ip < np; ++ip) {
    long p = nextprime(rand() % (1L << 23));
    modulus_init(&M, p);
    for (ia = 0; ia < na; ++ia) {
      long s = (rand() % 2) ? 1 : -1;
      long a = s * ((((unsigned long long)rand() << 32) | rand()) & (((unsigned long long)1 << 53)-1));
      double r1 = double_mods(a, M);
      double r2 = mods_check(a, M);
      if (r1 != r2) {
        printf("%ld %% %ld --> %f (%f)\n", a, p, r1, r2);
      }
      CuAssertIntEquals(tc, (long)r1, (long)r2);
    }
  }
}

static void test_gcdex(CuTest * tc)
{
  mpz_t a, b, g, u, v, s, t;
  mpz_init(a);
  mpz_init(b);
  mpz_init(g);
  mpz_init(u);
  mpz_init(v);
  mpz_init(s);
  mpz_init(t);
  mpz_set_ui(a, 162);
  mpz_set_ui(b, 63);

  gcdex(g, s, t, u, v, a, b);

  CuAssertIntEquals(tc, 9, mpz_get_si(g));
  CuAssertIntEquals(tc, -5, mpz_get_si(s));
  CuAssertIntEquals(tc, 13, mpz_get_si(t));
  CuAssertIntEquals(tc, -7, mpz_get_si(u));
  CuAssertIntEquals(tc, 18, mpz_get_si(v));

  mpz_set_ui(a, 105);
  mpz_set_ui(b, 20);

  gcdex(g, s, t, u, v, a, b);

  CuAssertIntEquals(tc, 5, mpz_get_si(g));
  CuAssertIntEquals(tc, -3, mpz_get_si(s));
  CuAssertIntEquals(tc, 16, mpz_get_si(t));
  CuAssertIntEquals(tc, -4, mpz_get_si(u));
  CuAssertIntEquals(tc, 21, mpz_get_si(v));

}


CuSuite * misc_suite(void)
{
  CuSuite * suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, test_nextprime);
  SUITE_ADD_TEST(suite, test_prevprime);
  SUITE_ADD_TEST(suite, test_modsL);
  SUITE_ADD_TEST(suite, test_modsL1);
  SUITE_ADD_TEST(suite, test_modsMpz);
  SUITE_ADD_TEST(suite, test_mods);
  SUITE_ADD_TEST(suite, test_gcdex);
  return suite;
}
