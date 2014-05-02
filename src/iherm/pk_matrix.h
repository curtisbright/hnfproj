#pragma once

#include "mpz_matrix.h"

/** \file pk_matrix.h
  * Packed triangular matrix.
  *
  * Specialized representation/routines for upper triangular matrices with few
  * non-trivial columns (non-identity columns) as arising in the Hermite normal
  * form computation.
  */

struct pk_t {
  mpz_t * data;
  long n;
  long k;
  long * idxs;
  long * rev_idxs;
  long * min_idxs;
};
typedef struct pk_t pk_t;

void pkMatrix_print(FILE * stream, pk_t const * T);

pk_t * pkMatrix_init(long n);
pk_t * pkMatrix_initFromIdxs(long * idxs, long k, long n);
void pkMatrix_reinit(pk_t * Z, long n);
void pkMatrix_realloc(pk_t * Z, long * idxs, long k, long n);
void pkMatrix_fini(pk_t * Z);

int pkMatrix_equal(pk_t const * A, pk_t const * B);

pk_t * pkMatrix_fromFull(mpzMatrix_t const * A);
mpzMatrix_t * pkMatrix_toFull(pk_t const * Z);

void pkMatrix_gemm(pk_t * Z, pk_t const * S, pk_t const * T);
void pkMatrix_gemm_iter(pk_t * Z, pk_t const * S, pk_t const * T);
void pkMatrix_gemm_block(pk_t * Z, pk_t const * S, pk_t const * T);

void pkMatrix_hnf(pk_t * Z);

void pkMatrix_applyVector_simple(mpzMatrix_t * X, pk_t const * Z, long start, long count);
void pkMatrix_applyVector_block(mpzMatrix_t * X, pk_t const * Z, long start, long count);

void pkMatrix_applyHinv(mpzMatrix_t * B, mpzMatrix_t const * A, pk_t const * H);

mpz_ptr pkget(pk_t * A, int row, int col);
mpz_srcptr pkgetc(pk_t const * A, int row, int col);
