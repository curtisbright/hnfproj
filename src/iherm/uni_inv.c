#include "uni_inv.h"

#include <assert.h>
#include <timer.h>

#include "imlsolve.h"

static void rankOneUpdate(mpzMatrix_t const * U, mpz_t d, mpzMatrix_t ** px, mpzMatrix_t ** py)
{
  mpz_t denom;
  long n = U->ncols;
  mpzMatrix_t * u = mpzMatrix_init(n, 1);
  mpzMatrix_t * v = mpzMatrix_init(1, n);
  mpzMatrix_t * x = mpzMatrix_init(n, 1);
  mpzMatrix_t * y = mpzMatrix_init(1, n);
  mpzMatrix_t * dd = mpzMatrix_init(1,1);

  mpz_init(denom);
  mpzMatrix_rand(u, 1);
  mpzMatrix_rand(v, 1);

  /* x := U^-1u */
  imlSolve(x, denom, U, u);
  assert(0 == mpz_cmp_ui(denom, 1));

  /* y := vU^-1 */
  imlSolveLeft(y, denom, U, v);
  assert(0 == mpz_cmp_ui(denom, 1));

  /* d:= 1 + vU^-1u = 1+vx */
  mpzMatrix_gemm(dd, v, x);
  mpz_set(d, mmgetc(dd, 0, 0));
  mpz_add_ui(d, d, 1);

  *px = x;
  *py = y;

  mpz_clear(denom);
  mpzMatrix_fini(dd);
  mpzMatrix_fini(u);
  mpzMatrix_fini(v);
}


static void mpzMatrix_iquos(mpzMatrix_t * A, mpz_t const dd)
{
  int i;
  mpzMatrix_t * B = mpzMatrix_init(A->nrows, A->ncols);
  mpz_t c, d;
  mpz_inits(c, d, 0);
  mpz_set_si(c, -1);
  mpz_abs(d, dd);

  /* A := A - (A mod d) */
  mpzMatrix_set(B, A);
  mpzMatrix_mods(B, d);
  mpzMatrix_scale(B, c);
  mpzMatrix_add(A, B);

  /* A:= A / d */
  for (i = 0; i < mpzMatrix_numElems(A); ++i) {
    mpz_divexact(A->data[i], A->data[i], dd);
  }

  mpzMatrix_fini(B);
  mpz_clears(c, d, 0);
}

static mpzMatrix_t * approxSoln(mpz_t const d, mpzMatrix_t const * x, mpzMatrix_t const * y, mpzMatrix_t const *b)
{
  mpzMatrix_t * t = mpzMatrix_init(x->nrows, y->ncols);
  mpzMatrix_t * z = mpzMatrix_init(x->nrows, b->ncols);

  /* z := (x.y quo d).b */
  TIMER("x.y",
  mpzMatrix_gemm(t, x, y);)
  TIMER("iquos",
  mpzMatrix_iquos(t, d);)
  TIMER("t.b",
  mpzMatrix_gemm(z, t, b);)

  mpzMatrix_fini(t);

  return z;
}

static mpzMatrix_t * checkSoln(mpzMatrix_t const * U, mpzMatrix_t const * z, mpzMatrix_t const *b)
{
  mpz_t c;
  mpzMatrix_t * chk = mpzMatrix_init(b->nrows, b->ncols);

  /* chk := b - Uz */
  mpz_init_set_si(c, -1);
  mpzMatrix_gemm(chk, U, z);
  mpzMatrix_scale(chk, c);
  mpzMatrix_add(chk, b);

  mpz_clear(c);
  return chk;
}


mpzMatrix_t * uniSolve(mpzMatrix_t const * U, mpzMatrix_t const * b)
{
  mpzMatrix_t *e = NULL, *x, *y, *z, *chk;
  mpz_t d, denom;
  mpz_inits(d, denom, 0);
  TIMER("rankOne",
  rankOneUpdate(U, d, &x, &y);)

  TIMER("approxSoln",
  z = approxSoln(d, x, y, b);)

  TIMER("checkSoln",
  chk = checkSoln(U, z, b);)

  printf("%s\n", mpzMatrix_isZero(chk) ? "True." : "False.");

  if(!mpzMatrix_isZero(chk)) {
    /* U(z+e) = b; Ue = chk */
    e = mpzMatrix_init(chk->nrows, chk->ncols);
    TIMER("e",
    imlSolve(e, denom, U, chk);)
    assert(0 == mpz_cmp_ui(denom, 1));
    mpzMatrix_add(z, e);
  }

  if(e) { mpzMatrix_fini(e); }
  mpzMatrix_fini(x);
  mpzMatrix_fini(y);
  mpzMatrix_fini(chk);
  mpz_clears(d, denom, 0);

  return z;
}
