#include "basis.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "gmp.h"

#include "arith_utils.h"

/** \brief Fill a new `modulus_t`. */
void modulus_init(modulus_t * M, long p)
{
  M->p = p;
  M->pd = (double)p;
  M->inv_p = 1.0/M->pd;
  M->half_p = M->pd/2.0;
}

static
long * calcBezout(mpz_t const * L, long const * primes, int nprimes)
{
  /* S_i = [(p_1..p_n)/p_i]^(-1) mod p_i */
  long i;
  long * S = (long *)malloc(nprimes * sizeof(long));
  mpz_t s, p;
  mpz_init(s);
  mpz_init(p);
  for (i = 0; i < nprimes; ++i) {
    mpz_set_ui(p, primes[i]);
    mpz_invert(s, L[i], p);
    S[i] = modsMpzL(s, primes[i]);
  }
  mpz_clear(s);
  mpz_clear(p);

  return S;
}

static
mpz_t * calcLagrange(mpz_t const total, long const * primes, int nprimes)
{
/* L_i = (p_1...p_n)/p_i */
  long i, p;
  mpz_t * L = (mpz_t *)malloc(nprimes * sizeof(mpz_t));
  for (i = 0; i < nprimes; ++i) {
    p = primes[i];
    mpz_init(L[i]);
    mpz_divexact_ui(L[i], total, p);
  }
  return L;
}

static
double * calcWeights(long const * S, long const * primes, int nprimes)
{
  long i;
  double * W = malloc(nprimes * sizeof(double));
  for (i = 0; i < nprimes; ++i) {
    W[i] = (double)(S[i])/primes[i];
  }
  return W;
}

static
mpz_t * calcRecon(mpz_t const * L, long const * S, int nprimes)
{
  long i;
  mpz_t * C = malloc(nprimes * sizeof(mpz_t));
  for (i = 0; i < nprimes; ++i) {
    mpz_init(C[i]);
    mpz_mul_si(C[i], L[i], S[i]);
  }
  return C;
}


#if 0
static
double * calcMixedInv(long const * primes, int nprimes)
{
  long i;
  mpz_t total;
  double * mixed_inv = (double *)malloc(nprimes * sizeof(double));

  mpz_init_set_ui(total, 1);
  for (i = 0; i < nprimes; ++i) {
    mixed_inv[i] = -modsL(modInverseMpz(total, primes[i]), primes[i]);
    mpz_mul_ui(total, total, primes[i]);
  }
  return mixed_inv;
}
#endif

#ifdef DEBUG
static int checkStartModulus(long n, long p)
{
  /* Check: n(p-1)^2 + (p-1) < 2^53 - 1. */
  int rc;
  mpz_t lhs, rhs;
  mpz_inits(lhs, rhs, 0);
  mpz_set_ui(lhs, n);
  mpz_mul_ui(lhs, lhs, p-1);
  mpz_mul_ui(lhs, lhs, p-1);
  mpz_add_ui(lhs, lhs, p-1);
  mpz_ui_pow_ui(rhs, 2, 53);
  mpz_sub_ui(rhs, rhs, 1);
  rc = (mpz_cmp(lhs, rhs) < 0);
  mpz_clears(lhs, rhs, 0);
  return rc;
}
#endif

static long * genPrimes(int start, mpz_t const bound, int * nprimes)
{
  long i = 0;
  mpz_t total;
  long ret_alloc = 8;
  long * ret = malloc(ret_alloc * sizeof(modulus_t));

  mpz_init_set_ui(total, 1);

  while (mpz_cmp(total, bound) < 0) {
    if (i == ret_alloc - 1) {
      ret_alloc *= 2;
      ret = realloc(ret, ret_alloc * sizeof(modulus_t));
    }
    start = prevprime(start);
    ret[i] = start;
    mpz_mul_ui(total, total, start);
    ++i;
  }

  mpz_clear(total);

  *nprimes = i;
  return ret;
}

/** \brief Largest modulus compatible with dimension `n` matrix.
 *
 * Returns the largest value `p` such that a residue number system backed with
 * BLAS operations having all moduli less than `p` can represent a matrix of
 * dimension `n`. As the BLAS operations use double-precision floating point
 * numbers, the result of a dot product between vectors of dimension `n` must
 * be exactly representable in the 53 bits of a `double`'s mantissa.
 *
 * Viz. \f$n(p-1)^2 + (p-1) < 2^{53} - 1.\f$; equivalently, \f$p < \left(\sqrt{4n(2^{53}-1)+1} + (2n-1)\right) / 2n \f$
 */
long pickStartModulus(long n)
{
  /* Choose p such that n(p-1)^2 + (p-1) < 2^{53} - 1.
   * i.e, p < [sqrt(4n(2^53-1)+1) + (2n-1)] / 2n
   * Alternate, looser bound: p < [ (1 << 26) * sqrt(2.0) / sqrt((double)n) ];
   */
  long p;
  mpz_t x;

  /* x = [sqrt(4n(2^53-1)+1) + (2n-1)] / 2n */
  mpz_init(x);
  mpz_ui_pow_ui(x, 2, 53);
  mpz_sub_ui(x, x, 1);
  mpz_mul_ui(x, x, 4*n);
  mpz_add_ui(x, x, 1);
  mpz_sqrt(x, x);
  mpz_add_ui(x, x, 2*n-1);
  mpz_fdiv_q_ui(x, x, 2*n);

  p = mpz_get_si(x);

  assert(mpz_fits_slong_p(x));
  assert(checkStartModulus(n, p));

  mpz_clear(x);
  return p;
}

/** Initialize a new basis_t.
 * \param primes array of coprime moduli
 * \param nprimes number of moduli
 *
 * Initializes an RNS basis with precomputed quantities useful in common
 * operations (e.g., Lagrange interpolants for reconstruction by CRA).  The
 * supplied moduli are assumed to be pairwise coprime.  The caller is
 * responsible for freeing the new basis with the `basis_fini` method.
 */
basis_t * basis_init(long const * primes, int nprimes)
{
  long i;
  basis_t * P = malloc(sizeof(basis_t));
  modulus_t * M = malloc(nprimes * sizeof(modulus_t));

  for (i = 0; i < nprimes; ++i) {
    modulus_init(&M[i], primes[i]);
  }

  mpz_init_set_ui(P->PP, 1);
  for (i = 0; i < nprimes; ++i) {
    mpz_mul_ui(P->PP, P->PP, primes[i]);
  }

  P->mod = M;
  P->nmod = nprimes;
  P->L = calcLagrange(P->PP, primes, nprimes);
  P->inv_L = calcBezout((mpz_t const *)P->L, primes, nprimes);
  P->W = calcWeights(P->inv_L, primes, nprimes);
  P->C = calcRecon((mpz_t const *)P->L, P->inv_L, nprimes);

  return P;
}

basis_t * basis_initFromBound(long start, mpz_t const bound)
{
  int nprimes;
  long * primes = genPrimes(start, bound, &nprimes);
  basis_t * P = basis_init(primes, nprimes);
  free(primes);
  return P;
}



/** \brief Free a basis_t.
*/
void basis_fini(basis_t * P)
{
  long i;
  free(P->mod);
  for (i = 0; i < P->nmod; ++i) {
    mpz_clear(P->L[i]);
    mpz_clear(P->C[i]);
  }
  free(P->C);
  free(P->L);
  free(P->inv_L);
  free(P->W);
  mpz_clear(P->PP);
  free(P);
}

/** \brief Deep basis copy */
basis_t * basis_copy(basis_t const * src)
{
  long i;
  basis_t * dst;
  long nprimes = src->nmod;
  long * primes = malloc(nprimes*sizeof(long));

  for (i = 0; i < nprimes; ++i) {
    primes[i] = basis_getPrime(src, i);
  }
  dst = basis_init(primes, nprimes);
  free(primes);

  return dst;
}

/** \brief Print in a simple human-readable format. */
void basis_print(FILE * stream, basis_t const * P)
{
  long i;
  for (i = 0; i < P->nmod; ++i) {
    fprintf(stream, "%ld ", P->mod[i].p);
  }
  fprintf(stream, "\n");
}

/** \brief Return `i`th basis modulus. */
long basis_getPrime(basis_t const * P, int i)
{
  assert(i >= 0 && i < P->nmod);
  return P->mod[i].p;
}

/** \brief Print in a detailed human-readable format. */
void basis_longPrint(FILE * stream, basis_t const * P)
{
  long i;
  long nP = P->nmod;

  fprintf(stream, "p:\t");
  for (i = 0; i < nP; ++i) { fprintf(stream, "%ld ", P->mod[i].p); }
  fprintf(stream, "\n");

  fprintf(stream, "L:\t");
  for (i = 0; i < nP; ++i) {
    mpz_out_str(stream, 10, P->L[i]);
    fprintf(stream, " ");
  }
  fprintf(stream, "\n");

  fprintf(stream, "inv_L:\t");
  for (i = 0; i < nP; ++i) { fprintf(stream, "%ld ", P->inv_L[i]); }
  fprintf(stream, "\n");

  fprintf(stream, "PP:\t");
  mpz_out_str(stream, 10, P->PP);
  fprintf(stream, "\n");

}


