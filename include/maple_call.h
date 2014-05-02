#pragma once

#include "gmp.h"
#include "maplec.h"

/*
 * Args:  1
 *        A - nxn RTable
 */
ALGEB M_DECL highOrderResidue_maple(MKernelVector kv, ALGEB args);
/*
 * Args:  1
 *        A - nxn RTable
 */
ALGEB M_DECL unicert_maple(MKernelVector kv, ALGEB args);

/*
 * Args:  2
 *        A - nxn RTable
 *        b - 1xn RTable
 */
ALGEB M_DECL horSolveIML_maple(MKernelVector kv, ALGEB args);
ALGEB M_DECL horSolveSpinv_maple(MKernelVector kv, ALGEB args);
ALGEB M_DECL imlSolve_maple(MKernelVector kv, ALGEB args);

ALGEB M_DECL maple_passing(MKernelVector kv, ALGEB args);

ALGEB M_DECL iherm_maple(MKernelVector kv, ALGEB args);
