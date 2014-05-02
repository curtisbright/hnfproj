
#include "maple_call.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "mpz_matrix.h"
#include "highorder.h"
#include "iherm.h"
#include "linsys.h"
#include "horsolve.h"
#include "imlsolve.h"
#include "timer.h"
#include "unicert.h"

/*static void whattype( MKernelVector kv, ALGEB s )
{
  if( IsMapleAssignedName(kv, s) )        MaplePrintf(kv, "AssignedName\n");
  if( IsMapleComplexNumeric(kv, s) )      MaplePrintf(kv, "ComplexNumeric\n");
  if( IsMapleComplex64(kv, s) )           MaplePrintf(kv, "Complex64\n");
  if( IsMapleNumeric(kv, s) )             MaplePrintf(kv, "Numeric\n");
  if( IsMapleFloat64(kv, s) )             MaplePrintf(kv, "Float64\n");
  if( IsMapleInteger(kv, s) )             MaplePrintf(kv, "Integer\n");
  if( IsMapleInteger8(kv, s) )            MaplePrintf(kv, "Integer8\n");
  if( IsMapleInteger16(kv, s) )           MaplePrintf(kv, "Integer16\n");
  if( IsMapleInteger32(kv, s) )           MaplePrintf(kv, "Integer32\n");
  if( IsMapleInteger64(kv, s) )           MaplePrintf(kv, "Integer64\n");
  if( IsMapleList(kv, s) )                MaplePrintf(kv, "List\n");
  if( IsMapleExpressionSequence(kv, s) )  MaplePrintf(kv, "ExprSeq\n");
  if( IsMapleName(kv, s) )                MaplePrintf(kv, "Name\n");
  if( IsMapleNULL(kv, s) )                MaplePrintf(kv, "Null\n");
  if( IsMaplePointer(kv, s) )             MaplePrintf(kv, "Pointer\n");
  if( IsMaplePointerNULL(kv, s) )         MaplePrintf(kv, "PointerNull\n");
  if( IsMapleProcedure(kv, s) )           MaplePrintf(kv, "Proc\n");
  if( IsMapleRTable(kv, s) )              MaplePrintf(kv, "RTable\n");
  if( IsMapleSet(kv, s) )                 MaplePrintf(kv, "Set\n");
  if( IsMapleStop(kv, s) )                MaplePrintf(kv, "Stop\n");
  if( IsMapleString(kv, s) )              MaplePrintf(kv, "String\n");
  if( IsMapleTable(kv, s) )               MaplePrintf(kv, "Table\n");
  if( IsMapleUnassignedName(kv, s) )      MaplePrintf(kv, "Unassignedname\n");
  if( IsMapleUnnamedZero(kv, s) )         MaplePrintf(kv, "UnnamedZero\n");
}*/

static mpzMatrix_t * mpzMatrix_fromRTable( MKernelVector kv, ALGEB rt )
{
  long src_idx, dst_idx = 0, nrows, ncols;
  mpzMatrix_t * A;
  ALGEB * A_mpl;
  mpz_srcptr tmp;
  RTableSettings rts;

  assert(IsMapleRTable(kv, rt));

  RTableGetSettings(kv, &rts, rt);
  assert(rts.subtype == RTABLE_MATRIX);
  assert(rts.num_dimensions == 2);

  nrows = RTableUpperBound(kv, rt, 1);
  ncols = RTableUpperBound(kv, rt, 2);
  assert(nrows*ncols == RTableNumElements(kv, rt));

  A_mpl = (ALGEB*)RTableDataBlock(kv, rt);

  A = mpzMatrix_init(nrows, ncols);

  for(src_idx = 0; src_idx < nrows*ncols; ++src_idx) {
      tmp = MapleToGMPInteger(kv, A_mpl[src_idx]);
      switch (rts.order) {
        case RTABLE_C:
          dst_idx = src_idx;
          break;
        case RTABLE_FORTRAN:
          dst_idx = (src_idx/nrows)+(src_idx%nrows)*ncols;
          break;
        default:
          assert(0);
          break;
      }
      mpz_set(A->data[dst_idx], tmp);
  }
  return A;
}

static ALGEB mpzMatrix_toRTable( MKernelVector kv, mpzMatrix_t const * A )
{
  long i;
  RTableSettings rts;
  M_INT bounds[4];
  ALGEB * dst;
  mpz_t const * src;
  ALGEB rt;

  RTableGetDefaults(kv, &rts);
  rts.order = RTABLE_C;
  rts.num_dimensions = 2;
  rts.subtype = RTABLE_MATRIX;
  rts.data_type = RTABLE_DAG;

  bounds[0] = 1;
  bounds[1] = A->m;
  bounds[2] = 1;
  bounds[3] = A->n;

  rt = RTableCreate(kv,&rts, NULL, bounds);
  dst = (ALGEB*)RTableDataBlock(kv, rt);
  src = mpzMatrix_constData(A);

  for (i = 0; i < A->m * A->n; ++i) {
    dst[i] = GMPIntegerToMaple(kv, (mpz_ptr)src[i]);
  }
  return rt;
}

static void * alloc_func(size_t sz) { return malloc(sz); }
static void * realloc_func(void * ptr, size_t old_sz, size_t new_sz) { (void)old_sz; return realloc(ptr, new_sz); }
static void free_func(void * ptr, size_t sz) { (void)sz; free(ptr); }

ALGEB M_DECL highOrderResidue_maple( MKernelVector kv, ALGEB args )
{
  int argc;
  mpzMatrix_t * A;
  ALGEB A_rt;
  mpzMatrix_t * R;
  ALGEB R_rt;

  argc = MapleNumArgs(kv,args);
  if (argc != 1 || !IsMapleRTable(kv, (ALGEB)args[1])) {
    MapleRaiseError(kv,"one argument expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }

  R = highOrderResidue(A);
  R_rt = mpzMatrix_toRTable(kv, R);

  mpzMatrix_fini(A);
  mpzMatrix_fini(R);

  MaplePopGMPAllocators(kv);

  return R_rt;
}

ALGEB M_DECL iherm_maple( MKernelVector kv, ALGEB args )
{
  int argc;
  mpzMatrix_t * A;
  ALGEB A_rt;
  mpzMatrix_t * H;
  ALGEB H_rt;

  argc = MapleNumArgs(kv,args);
  if (argc != 1 || !IsMapleRTable(kv, (ALGEB)args[1])) {
    MapleRaiseError(kv,"one argument expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }

  H = pkMatrix_toFull(myHermite(A));
  H_rt = mpzMatrix_toRTable(kv, H);

  mpzMatrix_fini(A);
  mpzMatrix_fini(H);

  MaplePopGMPAllocators(kv);

  return H_rt;
}

ALGEB M_DECL unicert_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  mpzMatrix_t * A;
  ALGEB A_rt;
  int rslt;

  argc = MapleNumArgs(kv,args);
  if (argc != 1 || !IsMapleRTable(kv, (ALGEB)args[1])) {
    MapleRaiseError(kv,"one argument expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }

  rslt = uniCert(A);
  mpzMatrix_fini(A);

  MaplePopGMPAllocators(kv);

  return ToMapleBoolean(kv, rslt);
}

ALGEB M_DECL horSolveIML_maple(MKernelVector kv, ALGEB args)
{
  int argc, Rzero;
  ALGEB A_rt, b_rt, y_rt, d_rt, Rzero_maple;
  mpzMatrix_t * A, * b;
  mpz_t d;
  mpzMatrix_t * y;

  argc = MapleNumArgs(kv,args);
  if (argc != 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  b_rt = (ALGEB)args[2];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }
  b = mpzMatrix_fromRTable(kv, b_rt);
  if (!b) { return NULL; }
  mpz_init(d);

  y = horSolveIML(d, &Rzero, A, b);

  y_rt = mpzMatrix_toRTable(kv, y);
  d_rt = GMPIntegerToMaple(kv, d);
  Rzero_maple = ToMapleBoolean(kv, Rzero);

  mpzMatrix_fini(A);
  mpzMatrix_fini(b);
  mpzMatrix_fini(y);
  mpz_clear(d);

  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 3, Rzero_maple, y_rt, d_rt);
}

ALGEB M_DECL horSolveSpinv_maple(MKernelVector kv, ALGEB args)
{
  int argc, Rzero;
  ALGEB A_rt, b_rt, y_rt, d_rt, Rzero_maple;
  mpzMatrix_t * A, * b;
  mpz_t d;
  mpzMatrix_t * y;

  argc = MapleNumArgs(kv,args);
  if (argc != 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  b_rt = (ALGEB)args[2];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }
  b = mpzMatrix_fromRTable(kv, b_rt);
  if (!b) { return NULL; }

  mpz_init(d);

  y = horSolveSpinv(d, &Rzero, A, b);

  y_rt = mpzMatrix_toRTable(kv, y);
  d_rt = GMPIntegerToMaple(kv, d);
  Rzero_maple = ToMapleBoolean(kv, Rzero);

  mpzMatrix_fini(A);
  mpzMatrix_fini(b);
  mpzMatrix_fini(y);
  mpz_clear(d);

  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 3, Rzero_maple, y_rt, d_rt);
}

ALGEB M_DECL imlSolve_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB A_rt, b_rt, y_rt, d_rt;
  mpzMatrix_t * A, * b;
  mpz_t d;
  mpzMatrix_t * y;

  argc = MapleNumArgs(kv,args);
  if (argc != 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  b_rt = (ALGEB)args[2];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }
  b = mpzMatrix_fromRTable(kv, b_rt);
  if (!b) { return NULL; }
  mpz_init(d);
  y = mpzMatrix_init(b->m, b->n);

  imlSolve(y, d, A, b);

  y_rt = mpzMatrix_toRTable(kv, y);
  d_rt = GMPIntegerToMaple(kv, d);

  mpzMatrix_fini(A);
  mpzMatrix_fini(b);
  mpzMatrix_fini(y);
  mpz_clear(d);

  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 2, y_rt, d_rt);
}

ALGEB M_DECL maple_passing(MKernelVector kv, ALGEB args)
{
  int argc;
  mpzMatrix_t * A;
  mpzMatrix_t * B;
  ALGEB A_rt, B_rt;

  argc = MapleNumArgs(kv,args);
  if (argc != 1 || !IsMapleRTable(kv, (ALGEB)args[1])) {
    MapleRaiseError(kv,"one argument expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  B = mpzMatrix_init(A->m, A->n);
  if (!A) { return NULL; }

  TIMER("mpz_set", mpzMatrix_set(B, A));

  B_rt = mpzMatrix_toRTable(kv, B);

  mpzMatrix_fini(A);
  mpzMatrix_fini(B);

  MaplePopGMPAllocators(kv);

  return B_rt;
}
