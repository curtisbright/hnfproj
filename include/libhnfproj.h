#pragma once

#include "gmp.h"

/* High order residue
 * Inputs:
 *      A - nxn integer matrix in row-major order
 *      n - matrix dimension
 * Returns:
 *      R - nxn matrix with A^{-1} = [*] + A^{-1}R and ||A^{-1}R|| small
 */
long * highOrderResidueLong(long const * A, long n);
mpz_t * highOrderResidueMpz(mpz_t const * A, long n);

/* Unimodularity certification
 * Inputs:
 *      A - nxn integer matrix in row-major order
 *      n - matrix dimension
 * Returns:
 *      1 if det(A) = +/- 1
 *      0 otherwise
 */
int uniCertLong(long const * A, long n);
int uniCertMpz(mpz_t const * A, long n);

/* Hermite normal form (square, non-singular).
 * Inputs:
 *      A - nxn integer matrix in row-major order
 *      n - matrix dimension
 * Returns:
 *      H - nxn matrix of GMP integers in Hermite normal form
 */
mpz_t * hermiteLong(long const * A, long n);
mpz_t * hermiteMpz(mpz_t const * A, long n);
