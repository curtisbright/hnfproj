#pragma once

#include <stdio.h>

#include "gmp.h"

/** \file basis.h
 * \brief Precomputed basis for RNS
 *
 */

struct modulus_t
{
  long p;
  double pd;
  double inv_p;
  double half_p;
};
typedef struct modulus_t modulus_t;

struct basis_t
{
  /**<\brief Modulus count */
  long nmod;

  /**<\brief Array of moduli */
  modulus_t * mod;

 /** \brief Product of moduli
  * \details \f$PP := \prod_{i=0}^{\mathrm{nmod}-1}p_i\f$
  */
  mpz_t PP;

 /** \brief Lagrange interpolants
  * \details \f$L_i := PP/p_i\f$ */
  mpz_t * L;

 /** \brief Inverse Lagrange interpolants
  * \details \f$\mathrm{inv\_L}_i := 1/L_i \bmod p_i \f$ */
  long * inv_L;

  /** \brief RNS weights
   * \details \f$W_i := \mathrm{inv\_L}_i/p_i \f$ */
  double * W;

  mpz_t * C;
};
typedef struct basis_t basis_t;

void modulus_init(modulus_t * M, long p);

basis_t * basis_init(long const * primes, int nprimes);

/** Initialize a new basis_t with moduli exceeding a given bound.
 * \param primes array of coprime moduli
 * \param nprimes number of moduli

 */
basis_t * basis_initFromBound(long start, mpz_t const bound);
void basis_fini(basis_t * P);

long basis_getPrime(basis_t const * P, int n);

basis_t * basis_copy(basis_t const * src);

void basis_longPrint(FILE * stream, basis_t const * P);
void basis_print(FILE * stream, basis_t const * P);

long pickStartModulus(long n);
