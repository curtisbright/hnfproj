#include "rns_matrix.h"

#include <assert.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "iml.h"

#include "residue.h"
#include "reconstruct.h"

/** \brief Allocate a new rnsMatrix_t.
 *
 * The new rnsMatrix_t has each element of each constituent matrix set to zero. The
 * caller is responsible for freeing the new matrix via the rnsMatrix_fini
 * method.
 *
 * \param nrows row dimension
 * \param ncols column dimension
 * \param P RNS basis
 */

rnsMatrix_t * rnsMatrix_init(long nrows, long ncols, basis_t const * P)
{
  long j;
  rnsMatrix_t * A = malloc(sizeof(rnsMatrix_t));
  A->nrows = nrows;
  A->ncols = ncols;
  A->basis = basis_copy(P);
  A->nmod = P->nmod;
  A->residues = malloc(A->nmod * sizeof(residue_t*));

  for (j = 0; j < A->nmod; ++j) {
    A->residues[j] = residue_init(nrows, ncols, basis_getPrime(P, j));
  }

  return A;

  #if 0
  /*double * block;
  size_t sz = n * n * sizeof(double);*/
  residues = malloc(P->nmod*sizeof(residue_t));
  /* Allocate a single block for all residue matrices.  Initialize A to
   * point into the single block appropriately.
   */
  block = (double*)malloc(M->nmod * sz);
  memset(block, 0, M->nmod * sz);

  for (j = 0; j < M->nmod; ++j) {
    residue_t * Mi = &(residues[j]);
    Mi->data = block;
    Mi->m = m;
    Mi->n = n;
    Mi->mod = P->mod[j];
    block += m*n;
  }
  M->residues = residues;
  #endif

}

/** \brief Free an rnsMatrix_t.
 */
void rnsMatrix_fini(rnsMatrix_t * A)
{
  long j;
  if (!A) { return; }
  for (j = 0; j < A->nmod; ++j) {
    residue_fini(A->residues[j]);
  }
  basis_fini(A->basis);
  free(A->residues);
  free(A);
}

double rnsMatrix_getEntry(rnsMatrix_t const * A, long res, long idx)
{
  assert(A->residues);
  assert(res < A->nmod);
  assert(idx < A->nrows*A->ncols);
  return residue_getEntry(A->residues[res], idx);
}

void rnsMatrix_setEntry(rnsMatrix_t * A, long res, long idx, double entry)
{
  assert(A->residues);
  assert(res < A->nmod);
  assert(idx < A->nrows*A->ncols);
  residue_setEntry(A->residues[res], idx, entry);
}

/** \brief Print matrix in a human-readable format.
 */
void rnsMatrix_print(rnsMatrix_t const * A)
{
  long nrows = A->nrows;
  long ncols = A->ncols;
  long nmod = A->nmod;

  long i, j, k;
  for (i = 0; i < nrows; i++) /* rows */
    {
      fprintf(stdout, "  ");
      fprintf(stdout, "  |  ");
      for (k = 0; k < nmod; k++) { /* residues */
        for (j = 0; j < ncols; j++) { /* columns */
          fprintf(stdout, "%4.0f  ", rnsMatrix_getEntry(A, k, i*ncols+j));
        }
      fprintf(stdout, "  |  ");
      }
      fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
}

/**************************/

/** \brief Deep matrix copy
 *
 * Copies the underlying array of residue_t matrices and their contents from
 * `src` to `dst`.  Source and destination matrices must be of equal dimension.
 */
void rnsMatrix_copy(rnsMatrix_t * dst, rnsMatrix_t const * src)
{
  long j;
  long nmod = dst->nmod;

  assert(dst->nrows == src->nrows);
  assert(dst->ncols == src->ncols);
  assert(dst->nmod == src->nmod);

  for (j = 0; j < nmod; ++j) {
    residue_copy(dst->residues[j], src->residues[j]);
  }
}

void rnsMatrix_determinant(mpz_t det, rnsMatrix_t * A)
{
  long j;
  long nmod = A->nmod;
  long * d = malloc(sizeof(long) * nmod);
  assert(A->nrows == A->ncols);

  for (j = 0; j < nmod; ++j) {
    d[j] = residue_determinant(A->residues[j]);
  }
  rns_reconstruct(det, A->basis, d);
}


/** \brief In-place inverse
 *
 * Compute the modular inverse of matrix `A` in place.
 *
 * \return
 *  - 1, if the inverse exists
 *  - 0, if the inverse does not exist
 */
int rnsMatrix_inverse(rnsMatrix_t * A)
{
  int rc;
  long j;
  long nmod = A->nmod;
  assert(A->nrows == A->ncols);

  for (j = 0; j < nmod; ++j) {
    rc = residue_inverse(A->residues[j]);
    if (!rc) { return 0; }
  }
  return 1;
}

/** \brief Set target to the identity matrix.
 *
 * The target matrix is assumed to be square.
 */
void rnsMatrix_identity(rnsMatrix_t * A)
{
  long j;
  long nmod = A->nmod;

  assert(A->nrows == A->ncols);

  for (j = 0; j < nmod; ++j) {
    residue_identity(A->residues[j]);
  }
}

/** \brief Check for the zero matrix.
 *
 * \return 1 if matrix A is the zero matrix; 0 otherwise.
 */
int rnsMatrix_isZero(rnsMatrix_t const * A)
{
  long j;
  for (j = 0; j < A->nmod; ++j) {
    if (!residue_isZero(A->residues[j])) { return 0; }
  }
  return 1;
}

#ifndef THREAD

/** \brief Matrix multiplication
 *
 * Set matrix `dst` to `A`*`B`.  Arguments cannot be aliased. Input matrices
 * are assumed to be of compatible dimension.
 */
void rnsMatrix_gemm(rnsMatrix_t * dst, rnsMatrix_t const * A, rnsMatrix_t const * B)
{
  long j;
  long nmod = dst->nmod;
  assert(dst->nmod == A->nmod);
  assert(A->nmod == B->nmod);
  assert(dst->nrows == A->nrows);
  assert(dst->ncols == B->ncols);
  assert(A->ncols == B->nrows);

  for (j = 0; j < nmod; ++j) {
    residue_gemm(dst->residues[j], A->residues[j], B->residues[j]);
  }
}

/** \brief Convenience quadratic lifting routine.
 *
 * `dst := Xinv.(T - A.M)`
 */
void rnsMatrix_quadLift(rnsMatrix_t * dst, rnsMatrix_t const * T,
                        rnsMatrix_t const * A, rnsMatrix_t const * M,
                        long const * Xinv)
{
  /* R = Xinv.(T - AM)*/
  long j;
  long nmod = dst->nmod;
  for (j = 0; j < nmod; ++j) {
    residue_quadLift(dst->residues[j], T->residues[j],
                     A->residues[j], M->residues[j], Xinv[j]);
  }
}

#else
struct gemm_args {
  residue_t * dst;
  residue_t const * A;
  residue_t const * B;
};
typedef struct gemm_args gemm_args_t;

static void * residue_gemm_thread(void * arg)
{
  gemm_args_t * args = (gemm_args_t *)arg;
  residue_gemm(args->dst, args->A, args->B);
  return NULL;
}

void rnsMatrix_gemm(rnsMatrix_t * dst, rnsMatrix_t const * A, rnsMatrix_t const * B)
{
  gemm_args_t args[32];
  pthread_t thread[32];

  long j;
  long nmod = dst->nmod;
  for (j = 0; j < nmod; ++j) {
    args[j].dst = dst->residues[j];
    args[j].A = A->residues[j];
    args[j].B = B->residues[j];

    pthread_create(&thread[j], NULL, residue_gemm_thread, &(args[j]));
  }
  for (j = 0; j < nmod; ++j) {
    pthread_join(thread[j], NULL);
  }
}

struct quadLift_args {
  residue_t * dst;
  residue_t const * T;
  residue_t const * A;
  residue_t const * M;
  long const * Xinv;
};
typedef struct quadLift_args quadLift_args_t;

static void * residue_quadLift_thread(void * arg)
{
  quadLift_args_t * args = (quadLift_args_t *)arg;
  residue_quadLift(args->dst, args->T, args->A, args->M, *args->Xinv);
  return NULL;
}

void rnsMatrix_quadLift(rnsMatrix_t * dst, rnsMatrix_t const * T,
                        rnsMatrix_t const * A, rnsMatrix_t const * M,
                        long const * Xinv)
{
  quadLift_args_t args[32];
  pthread_t thread[32];

  long j;
  long nmod = dst->nmod;
  for (j = 0; j < nmod; ++j) {
    args[j].dst = dst->residues[j];
    args[j].T = T->residues[j];
    args[j].A = A->residues[j];
    args[j].M = M->residues[j];
    args[j].Xinv = &Xinv[j];

    pthread_create(&thread[j], NULL, residue_quadLift_thread, &(args[j]));
  }
  for (j = 0; j < nmod; ++j) {
    pthread_join(thread[j], NULL);
  }
}
#endif


