#pragma once

#include "rns_matrix.h"
#include "mpz_matrix.h"

/** \file reconstruct.h
 * RNS/mpz matrix conversion.
 *
 * Logically, these methods belong with those in mpz_matrix.c or rns_matrix.c;
 * placing them there, however, so creates a circular dependency between the
 * corresponding headers.
 */

/** Matrix multiplication: rnsMatrix_t by mpzMatrix_t */
void mpzMatrix_gemmRnsMpz(mpzMatrix_t * dst, rnsMatrix_t const * A, mpzMatrix_t const * B);

/** Reconstruct matrix from residue number system. */
void mpzMatrix_reconstruct(mpzMatrix_t * A, rnsMatrix_t const * Ap);

/** Reconstruct arbitrary precision integer from residues. */
void rns_reconstruct(mpz_t ret, basis_t const * P, long * d);

/** Reduce matrix to residue number system. */
void rnsMatrix_fromMpzMatrix(rnsMatrix_t * Ap, mpzMatrix_t const * A);

