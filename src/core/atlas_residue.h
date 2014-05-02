#pragma once

#include <inttypes.h>

#include "basis.h"
#include "mpz_matrix.h"

/** \file atlas_residue.h
 * BLAS-compatible matrix
 *
 * An `atlasResidue_t` is a matrix of double-precision floating point numbers
 * representing an integer matrix reduced in the symmetric range with respect
 * to an appropriately-sized modulus. The operations defined here are
 * implemented in terms of ATLAS/BLAS routines.  The associated modulus is
 * selected such that the results of underlying BLAS operations can always be
 * represented exactly in the mantissa (53 bits) of a `double` (see basis.h).
 */
struct atlasResidue_t {
  long nrows; /**< row dimension */
  long ncols; /**< column dimension */
  double * _data; /**< data array */
  modulus_t mod; /**< corresponding modulus */
};
typedef struct atlasResidue_t atlasResidue_t;

/** Allocate a new `atlasResidue_t`.
 *
 * The caller is responsible for freeing the new matrix via the `atlas_fini` method.
 * \param nrows row dimension
 * \param ncols column dimension
 * \param p modulus
 */
atlasResidue_t * atlas_init(long nrows, long ncols, long p);

/** Free an `atlasResidue_t` */
void atlas_fini(atlasResidue_t * A);

/* Entry access */
long atlas_getEntry(atlasResidue_t const * A, long idx);
void atlas_setEntry(atlasResidue_t * A, long idx, long val);

/* Setters */
/** Deep matrix copy.
 * Copies the underlying array of doubles from `src` to `dst`.  Source and
 * destination matrices must be of equal dimension.
 */
void atlas_copy(atlasResidue_t * dst, atlasResidue_t const * src);

/** Reduction from arbitrary precision.
 *  Set matrix `dst` to `src` reduced in the symmetric range with respect to `dst`'s modulus.  Matrices are assumed to be of equal dimension.
 */
void atlas_fromMpzMatrix(atlasResidue_t * dst, mpzMatrix_t const * src);

/** Set target to the identity matrix.
 * The target matrix is assumed to be square.
 */
void atlas_identity(atlasResidue_t * A);

/** Zero the target matrix. */
void atlas_zero(atlasResidue_t * A);

/* Query operations*/

/** Check for the zero matrix.
 * \return 1 if matrix A is the zero matrix; 0 otherwise.
 */
int atlas_isZero(atlasResidue_t const * A);

/** Print matrix in a human-readable format.
 */
void atlas_print(FILE * stream, atlasResidue_t const * A);

/* Mutating operations */
/** Scaled matrix addition.
 * Set matrix `dst` to `dst`+ `k`*`src`.
 */
void atlas_add(atlasResidue_t * dst, atlasResidue_t const * src, long k);

/** In-place modular determinant.
 * The input matrix is destroyed in place.
 */
long atlas_determinant(atlasResidue_t * A);

/** In-place modular inverse.
 * Compute the modular inverse of matrix `A` in place.
 * \return
 *  - 1, if the inverse exists
 *  - 0, if the inverse does not exist
 */
int atlas_inverse(atlasResidue_t * A);

/** Matrix multiplication.
 * Set matrix `dst` to `A`*`B`.  Arguments cannot be aliased. Input matrices
 * are assumed to be of compatible dimension.
 */
void atlas_gemm(atlasResidue_t * dst, atlasResidue_t const * A, atlasResidue_t const * B);

/** Symmetric modular reduction.
 * Reduce (in place) the entries of `dst` in the symmetric range.
 */
void atlas_mods(atlasResidue_t * dst);

/** Convenience quadratic lifting routine.
 * `dst := Xinv.(T - A.M)`
 */
void atlas_quadLift(atlasResidue_t * dst, atlasResidue_t const * T, atlasResidue_t const * A, atlasResidue_t const * M, long xinv);

/** Scalar multiplication.
 * Scale matrix `dst` in place by scalar `k`.
 */
void atlas_scale(atlasResidue_t * dst, long k);

