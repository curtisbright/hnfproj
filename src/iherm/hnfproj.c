#include "hnfproj.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "arith_utils.h"
#include "thresholds.h"
#include "timer.h"


static long hcol_diagonal(mpz_t ** H, mpz_t ** t, mpzMatrix_t const * x, mpz_t const d)
{
  long n = x->nrows;
  long i, k = 0;
  mpz_t g, g_tmp, s, v;
  mpz_inits(g, g_tmp, s, v, 0);
  mpz_set(g, d);

  for (i = n-1; i >= 0; --i) {
    /* g:= gcdex(g, x[i]) */
    gcdex(g_tmp, s, (*t)[i], v, (*H)[i], g, x->data[i]);
    mpz_set(g, g_tmp);

    if (mpz_cmp_ui((*H)[i], 1) != 0) { ++k; }
  }

  mpz_clears(g, g_tmp, s, v, 0);

  return k;
}

static long * hcol_idxsFromDiagonal(mpz_t const * H, long k, long n)
{
  long i, j = 0;
  long * idxs = malloc(k * sizeof(long));
  for (i = 0; i < n; ++i) {
    if (mpz_cmp_ui(H[i], 1) != 0) {
      idxs[j] = i;
      ++j;
    }
  }
  assert(j == k);
  return idxs;
}

static void hcol_buildMatrix(pk_t * Z, mpzMatrix_t * x, mpz_t const dd, mpz_t const * H, mpz_t const * t)
{
  long n = Z->n;
  long i, j, q;
  mpz_t d;
  mpz_init_set(d, dd);
  for (j = 0; j < Z->k; ++j) {
    q = Z->idxs[j];

    /* d := d / H[q] */
    mpz_tdiv_q(d, d, H[q]);

    /* Z[q, q] := H[q]; */
    mpz_set(pkget(Z, q, j), H[q]);

    for (i = 0; i < q; ++i) {
      /* Z[i,q] := -t[q]*x[i] mod H[q] */
      mpz_mul(pkget(Z, i, j), t[q], x->data[i]);
      mpz_neg(pkget(Z, i, j), pkget(Z, i, j));
      mpz_mod(pkget(Z, i, j), pkget(Z, i, j), H[q]);

      /* x[i] := (x[i] + Z[i,q]*x[q]) / H[q] mod d*/
      mpz_addmul(x->data[i], pkget(Z, i, j), x->data[q]);
      mpz_tdiv_q(x->data[i], x->data[i], H[q]);
      mpz_mod(x->data[i], x->data[i], d);
    }
    for (i = q+1; i < n; ++i) {
      /* x[i] := x[i] / H[q] */
      mpz_tdiv_q(x->data[i], x->data[i], H[q]);
    }
  }
  mpz_clear(d);
}

static void hcol(pk_t * Z, mpzMatrix_t const * w, mpz_t const d)
{
  long i, k, n = w->nrows;
  long * idxs;
  mpz_t * H, * t;
  mpzMatrix_t * x = mpzMatrix_init(w->nrows, w->ncols);
  mpzMatrix_set(x, w);

  assert(w->ncols == 1);

  t = malloc(n * sizeof(mpz_t));
  H = malloc(n * sizeof(mpz_t));
  for (i = 0; i < n; ++i) {
    mpz_init(t[i]);
    mpz_init(H[i]);
  }

  k = hcol_diagonal(&H, &t, x, d);
  idxs = hcol_idxsFromDiagonal((mpz_t const *)H, k, n);

  pkMatrix_realloc(Z, idxs, k, n);
  hcol_buildMatrix(Z, x, d, (mpz_t const *)H, (mpz_t const *)t);

  for (i = 0; i < n; ++i) {
    mpz_clear(t[i]);
    mpz_clear(H[i]);
  }
  free(t);
  free(H);
  mpzMatrix_fini(x);
}

static void overshoot(mpz_t d, mpzMatrix_t * x)
{
  long i;
  mpz_t g;
  mpz_init_set(g, d);

  /* g:= gcd(d, x_1, x_2, ..., x_n)*/
  for (i = 0; i < x->nrows; ++i) {
    mpz_gcd(g, g, x->data[i]);
  }
  /* d := (d/g) */
  mpz_divexact(d, d, g);

  /* x:= (x/g) */
  for (i = 0; i < x->nrows; ++i) {
    mpz_divexact(x->data[i], x->data[i], g);
  }

  mpz_clear(g);
}

struct hSlices {
    pk_t ** L;
    long * count;
    long len;
    long nrows;
    long nalloc;
};
typedef struct hSlices hSlices_t;

static hSlices_t * hSlices_init(long nrows, long ncols)
{
  long c;
  hSlices_t * s = malloc(sizeof(hSlices_t));
  s->L = malloc(ncols * sizeof(pk_t*));
  s->count = malloc(ncols * sizeof(long));
  memset(s->count, 0, ncols * sizeof(long));
  s->len = 0;
  s->nalloc = ncols;
  s->nrows = nrows;
  for (c = 0; c < ncols; ++c) {
    s->L[c] = pkMatrix_init(nrows);
  }
  return s;
}

static void hSlices_fini(hSlices_t * s)
{
  long i;
  for (i = 0; i < s->nalloc; ++i) {
    pkMatrix_fini(s->L[i]);
  }
  free(s->count);
  free(s->L);
  free(s);
}

static void hSlices_collapse(hSlices_t * s)
{
  pk_t * last, * next;
  pk_t * Ltmp = pkMatrix_init(s->nrows);
  pk_t tmp;
  long * last_count, * next_count;

  while (s->len >= 2) {
    last_count = &(s->count[s->len - 1]);
    next_count = &(s->count[s->len - 2]);
    last = s->L[s->len - 1];
    next = s->L[s->len - 2];

    if (*last_count != *next_count) { break; }

    pkMatrix_gemm(Ltmp, last, next);
    tmp = *Ltmp; *Ltmp = *next; *next = tmp;

    *next_count += *last_count;
    *last_count = 0;
    s->len -= 1;
  }
  pkMatrix_fini(Ltmp);
}

static void hSlices_collapseAll(pk_t * H, hSlices_t * s)
{
  long i;
  pk_t tmp;
  pk_t * Htmp = pkMatrix_init(s->nrows);
  pkMatrix_reinit(H, s->nrows);

  for (i = 0; i < s->len; ++i) {
    pkMatrix_gemm(Htmp, s->L[i], H);
    tmp = *Htmp; *Htmp = *H; *H = tmp;
  }
  pkMatrix_fini(Htmp);
}

static void applyVector(hSlices_t * s, mpzMatrix_t * X, long start)
{
  long count;
  pk_t * Z;
  if (s->len == 0) { return; }

  Z = s->L[s->len-1];

  count = s->count[s->len-1];
  if (X->ncols - start < count) { count = X->ncols - start;}

  if( Z->k == 0 || count == 0) { return; }

  if (count >= PK_APPLY_VECTOR_BLOCK_THRESHOLD) {
    pkMatrix_applyVector_block(X, Z, start, count);
  } else {
    pkMatrix_applyVector_simple(X, Z, start, count);
  }
}

static void set_column(mpzMatrix_t * dst, mpzMatrix_t const * src, long col)
{
  long i;
  assert(dst->nrows == src->nrows);
  assert(dst->ncols == 1);
  assert(col < src->ncols);
  for (i = 0; i < dst->nrows; ++i) {
    mpz_set(mmget(dst, i, 0), mmgetc(src, i, col));
  }
}

void hermiteOfProjection_balanced(pk_t * H, mpzMatrix_t * X, mpz_t const dd)
{
  TIMER_ACC_INIT(Hx);
  TIMER_ACC_INIT(hcol);
  TIMER_ACC_INIT(gemm);

  mpzMatrix_t * y = mpzMatrix_init(X->nrows, 1);
  long col;
  mpz_t d;
  hSlices_t * s = hSlices_init(X->nrows, X->ncols);

  mpz_init(d);

  for (col = 0; col < X->ncols; ++col) {
    mpz_set(d, dd);

    TIMER_ACC(Hx,
    applyVector(s, X, col);
    set_column(y, X, col);
    mpzMatrix_mods(y, d);)
    overshoot(d, y);

    /*if(mpz_cmp_ui(d, 1) == 0) { continue; }*/
    TIMER_ACC(hcol,
    hcol(s->L[s->len], y, d);)
    s->count[s->len] = 1;
    s->len += 1;

    TIMER_ACC(gemm,
    hSlices_collapse(s);)
  }
  TIMER_ACC_FINI(Hx);
  TIMER_ACC_FINI(hcol);
  TIMER_ACC_FINI(gemm);

  TIMER("gemm2",
  hSlices_collapseAll(H, s);
  )

  TIMER("hnf",
  pkMatrix_hnf(H);)

  hSlices_fini(s);
  mpzMatrix_fini(y);
  mpz_clear(d);
}

void hermiteOfProjection_simple(pk_t * H, mpzMatrix_t * X, mpz_t const dd)
{
  TIMER_ACC_INIT(Hx);
  TIMER_ACC_INIT(hcol);
  TIMER_ACC_INIT(gemm);

  long col;
  mpzMatrix_t * y = mpzMatrix_init(X->nrows, 1);
  pk_t * Htmp = pkMatrix_init(X->nrows);
  pk_t * T = pkMatrix_init(X->nrows);
  pk_t tmp;
  mpz_t d;
  mpz_init(d);

  /* H := I_n */
  pkMatrix_reinit(H, X->nrows);
  for (col = 0; col < X->ncols; ++col) {

    mpz_set(d, dd);
    /* y := H.x mod d */
    TIMER_ACC(Hx,
    pkMatrix_applyVector_simple(X, H, col, 1);
    set_column(y, X, col);
    mpzMatrix_mods(y, d);)

    /* Remove gcd of entries in y.*/
    overshoot(d, y);
    if(mpz_cmp_ui(d, 1) == 0) { continue; }

    /* Minimal triangular denominator T. */
    TIMER_ACC(hcol,
    hcol(T, y, d);)

    /* H := T.H */
    TIMER_ACC(gemm,
    pkMatrix_gemm(Htmp, T, H);)
    tmp = *Htmp; *Htmp = *H; *H = tmp;
  }

  TIMER_ACC_FINI(Hx);
  TIMER_ACC_FINI(hcol);
  TIMER_ACC_FINI(gemm);

  TIMER("hnf",
  pkMatrix_hnf(H);)

  mpzMatrix_fini(y);
  pkMatrix_fini(Htmp);
  pkMatrix_fini(T);
  mpz_clear(d);
}
