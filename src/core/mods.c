
#include <assert.h>
#include <math.h>

#ifdef DEBUG
static double mods_check(double n, modulus_t M)
{
  double r = fmod(n, M.pd);
  if (r < -(M.pd/2.0)) {
    r += M.pd;
  } else if (r > M.pd/2.0) {
    r -= M.pd;
  }
  return r;
}
#endif

static double double_mods(double n, modulus_t M)
{
  long q = n*(M.inv_p);
  double r = n - (M.pd)*q;
  if (r < -M.half_p) {
    r += M.pd;
  } else if (r > M.half_p) {
    r -= M.pd;
  }
  assert(r == mods_check(n, M));

  return r;
}

#if 0
static double mods2(double nn, modulus_t M)
{
  double n = fabs(nn);
  long q = n*(M.inv_p);
  double r = n - (M.pd)*q;
  if (r > M.half_p) {
    r -= M.pd;
  }
  return r;
}

static double double_mods(double n, modulus_t M)
{
  double qd = (n*(M.inv_p));
  double r;
  long q;
  asm("cvtsd2siq %1, %0"
    :"=r"(q)
    :"xm"(qd));
  r = n - (M.pd)*q;

  assert(r == mods_check(n, M));

  return r;
}

static double double_mods(double n, modulus_t M)
{
  /* Compute quotient by multiplication by reciprocal; compute remainder
   * by subsequent subtraction.
   * Cast from double to long truncates towards 0.  Adding a shift allows
   * cheap rounding.
   */
  double s = (n > 0) ? 0.5 : -0.5;
  double qd = n*M.inv_p;
  long q = (qd + s);
  double r = n - (M.pd)*q;

  assert(r == mods_check(n, M));

  return r;
}
#endif

#if 0
static double double_modp(double n, modulus_t M)
{
  /* Compute quotient by multiplication by reciprocal; compute remainder
   * by subsequent subtraction.  Adjust to get 0 <= r < p.
   */

  double qd = n*M.inv_p;
  long q = qd;
  double r = n - (M.pd)*q;
  r = (r < 0) ? r + M.pd : r;

  return r;
}
#endif
