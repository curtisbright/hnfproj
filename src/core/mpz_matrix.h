#pragma once

#include <stdio.h>

#include "gmp.h"

/** \file mpz_matrix.h
 * Arbitrary precision integer matrix.
 *
 * An `mpzMatrix_t` is a matrix of arbitrary precision integers.  The
 * operations defined here are implemented (mostly naively) in terms of the
 * appropriate GMP functions.
 */

struct mpzMatrix_t {
  mpz_t * data;
  long nrows;
  long ncols;
};
typedef struct mpzMatrix_t mpzMatrix_t;

/** Allocate a new `mpzMatrix_t`.
 *
 * The new matrix has each element set to zero. The caller is responsible for
 * freeing the new matrix via the `mpzMatrix_fini` method.
 *
 * \param nrows row dimension
 * \param ncols column dimension
 */
mpzMatrix_t * mpzMatrix_init(long nrows, long ncols);

/** Free an `mpzMatrix_t` */
void mpzMatrix_fini(mpzMatrix_t * M);

/** Initialize a new mpzMatrix_t from an array of mpz_t integers.
 * The supplied array is assumed to be in row-major order with nrows*ncols elements.
 */
mpzMatrix_t const * mpzMatrix_initFromMpz(long nrows, long ncols, mpz_t const * data);

/** Initialize a new mpzMatrix_t from an array of long integers.
 * The supplied array is assumed to be in row-major order with nrows*ncols elements.
 */
mpzMatrix_t * mpzMatrix_initSet(long nrows, long ncols, long const * M_long);

/** Initialize a new mpzMatrix_t from a file in MatrixMarket format. */
mpzMatrix_t * mpzMatrix_initFromFile(char const * filename);
/** Write matrix to a MatrixMarket file. */
void mpzMatrix_writeToFile(char const * filename, mpzMatrix_t const * M);


/** Return a constant pointer to the underlying array of mpz_t integers */
mpz_t const * mpzMatrix_constData(mpzMatrix_t const * M);

/** Return a mutable pointer to the underlying array of mpz_t integers*/
mpz_t * mpzMatrix_data(mpzMatrix_t * M);

/** Entry access (mutable). */
mpz_ptr mmget(mpzMatrix_t * A, int row, int col);
/** Entry access (constant). */
mpz_srcptr mmgetc(mpzMatrix_t const * A, int row, int col);

/** Return the number of matrix entries */
long mpzMatrix_numElems(mpzMatrix_t const * M);


/** Set target to the identity matrix.
 *
 * The target matrix is assumed to be square.
 */
void mpzMatrix_identity(mpzMatrix_t * M);

/** Fill matrix with random entries.
 *
 * \param M Destination matrix
 * \param l Bit length of random entries
 */
void mpzMatrix_rand(mpzMatrix_t * M, long l);

/** Matrix assignment.
 *
 * Set entries of dst from those of src. Destination and source matrices are
 * assumed to be of the same dimensions.
 */
void mpzMatrix_set(mpzMatrix_t * dst, mpzMatrix_t const * src);

/** Efficient matrix swap.
 *
 * Efficiently swap - i.e., without copying - the contents of A and B.  The
 * swap operation only requires an exchange of the underlying data pointers.
 */
void mpzMatrix_swap(mpzMatrix_t * A, mpzMatrix_t * B);

/** Zero the target matrix. */
void mpzMatrix_zero(mpzMatrix_t * M);


/** Matrix equality.
 *
 * \return 1 if matrices A and B are equal; 0 otherwise.
 */
int mpzMatrix_equal(mpzMatrix_t const * A, mpzMatrix_t const * B);

/** Check for the zero matrix.
 *
 * \return 1 if matrix A is the zero matrix; 0 otherwise.
 */
int mpzMatrix_isZero(mpzMatrix_t const * A);

/** Return maximum matrix entry.
 *
 * Return (as 'ret') the largest entry in absolute value of A.
 */
void mpzMatrix_max(mpz_t ret, mpzMatrix_t const * A);

/** Print matrix in a human-readable format.
 */
void mpzMatrix_print(FILE * stream, mpzMatrix_t const * A);


/** Matrix addition.
 *
 * Set matrix `dst` to `dst`+ `src`.
 */
void mpzMatrix_add(mpzMatrix_t * dst, mpzMatrix_t const * src);

/** Scaled matrix addition.
 *
 * Set matrix `dst` to `dst`+ `c`*`src`.
 */
void mpzMatrix_addmul(mpzMatrix_t * dst, mpzMatrix_t const * src, mpz_t const c);

/** Simple matrix determinant. */
void mpzMatrix_determinant(mpz_t det, mpzMatrix_t const * M);

/** Matrix multiplication.
 *
 * Set matrix `dst` to `A`*`B`.  Arguments cannot be aliased. Input matrices
 * are assumed to be of compatible dimension.
 */
void mpzMatrix_gemm(mpzMatrix_t * dst, mpzMatrix_t const * A, mpzMatrix_t const * B);
/** Matrix multiplication via residue number system. */
void mpzMatrix_rnsGemm(mpzMatrix_t * C, mpzMatrix_t * const A, mpzMatrix_t * const B);

/** Symmetric modular reduction.
 *
 * Reduce (in place) the entries of `dst` in the symmetric range modulo `p`.
 */
void mpzMatrix_mods(mpzMatrix_t * dst, mpz_t const p);

/** Scalar multiplication.
 *
 * Scale matrix `dst` (in place) by scalar `c`.
 */
void mpzMatrix_scale(mpzMatrix_t * dst, mpz_t const c);

