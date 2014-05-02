#pragma once

#include <inttypes.h>

#include "basis.h"
#include "residue.h"

/** \file rns_matrix.h
 * \brief Residue number system matrix.
 *
 * An rnsMatrix_t is a represents an integer matrix in a residue number system
 * as a collection of (residue_t) matrices with word-sized entries and the RNS'
 * associated basis (a basis_t).
 *
 */

struct rnsMatrix_t {
  long nrows;             /**<\brief row dimension */
  long ncols;             /**<\brief column dimension */
  long nmod;              /**<\brief residue count */
  residue_t ** residues;  /**<\brief array of residues */
  basis_t * basis;  /**<\brief corresponding RNS basis */
};

typedef struct rnsMatrix_t rnsMatrix_t;

/* Memory management */
rnsMatrix_t * rnsMatrix_init(long nrows, long ncols, basis_t const * P);
void rnsMatrix_fini(rnsMatrix_t * M);

/* Entry access */
double rnsMatrix_getEntry(rnsMatrix_t const * M, long res, long idx);
void rnsMatrix_setEntry(rnsMatrix_t * M, long res, long idx, double entry);

/* Setters */
void rnsMatrix_copy(rnsMatrix_t * dst, rnsMatrix_t const * src);
void rnsMatrix_identity(rnsMatrix_t * M);

/* Query operations */
int rnsMatrix_isZero(rnsMatrix_t const * M);
void rnsMatrix_print(rnsMatrix_t const * A);

/* Mutating operations */
void rnsMatrix_determinant(mpz_t det, rnsMatrix_t * A);
void rnsMatrix_gemm(rnsMatrix_t * dst, rnsMatrix_t const * A,
                    rnsMatrix_t const * B);
int rnsMatrix_inverse(rnsMatrix_t * M);
void rnsMatrix_quadLift(rnsMatrix_t * dst, rnsMatrix_t const * T,
                        rnsMatrix_t const * A, rnsMatrix_t const * M,
                        long const * Xinv);

