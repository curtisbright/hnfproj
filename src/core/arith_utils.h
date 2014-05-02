#pragma once

#include "gmp.h"

/** \file arith_utils.h
 * Arithmetic helpers.
 *
 * Convenience wrappers to GMP number theoretic and arithmetic functions for
 * word-sized arguments.
 */

/** Scaled addtion: ret := ret + a*b */
void addmul_si(mpz_t ret, mpz_t const a, long b);

/** Extended gcd.
  * Return g = gcd(a,b) and s, t, u, v such that   \verbatim
   [ g ] = [ s  t ][ a ]
   [ 0 ]   [ u  v ][ b ]  \endverbatim
   v > 0 and 0 <= t < v.
*/
void gcdex(mpz_t g, mpz_t s, mpz_t t, mpz_t u, mpz_t v, mpz_t const a, mpz_t const b);

/** Modular inverse (long) */
long modInverseL (long nn, long pp);

/** Modular inverse (mpz) */
long modInverseMpz (mpz_t const n, long pp);

/** Symmetric modular reduction (long) */
long modsL(long r, long p);

/** Symmetric modular reduction (mpz) */
long modsMpzL(mpz_t const r, long p);

/** Symmetric modular reduction (mpz modulus) */
void modsMpzMpz(mpz_t r, mpz_t const p);

/** Next largest prime */
long nextprime(long nn);

/** Next smallest prime */
long prevprime(long nn);

