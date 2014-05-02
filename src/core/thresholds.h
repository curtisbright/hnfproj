#pragma once

/** Column dimension above which packed matrix multiplication switches from
 * iterative to blocked method (see src/iherm/pk_matrix.h). */
extern long PK_GEMM_BLOCK_THRESHOLD;

/** Column dimension above which applying a packed matrix to a dense block of
 * vectors switches from iterative to blocked method.*/
extern long PK_APPLY_VECTOR_BLOCK_THRESHOLD;

/** Number of basis elements in residue number system below which conversion
 * switches to ``fancy'' method from naive reconstruct/reduce method (see
 * src/lift/rns_conversion.c). */
extern long RNS_CONV_FANCY_THRESHOLD;
