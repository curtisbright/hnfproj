#pragma once

#include <stdio.h>
#include <time.h>

/** \file timer.h
 * \brief Timing convenience macros.
 */

extern int g_timer_depth;
extern int g_timer_threshold;
extern FILE * g_timer_stream;

#ifdef NOTIMER

#define TIMER(lbl, stmt) stmt
#define REAL_TIMER(lbl, stmt) stmt
#define TIMER_ACC_INIT(lbl)
#define TIMER_ACC(lbl, stmt) stmt
#define TIMER_ACC_FINI(lbl)

#else

#define TIMER(lbl, stmt) \
{ \
  static clock_t t0, t1; \
  static unsigned long long dt; \
  g_timer_depth += 1; \
  t0 = clock(); \
  { stmt } \
  t1 = clock(); \
  g_timer_depth -= 1; \
  dt = (unsigned long long)(t1 - t0) * 1000 / (unsigned long long)CLOCKS_PER_SEC; \
  timer_print(g_timer_depth, lbl, dt); \
}

#define TIMERF(lbl, stmt, ...) \
{ \
  char buf[100]; \
  snprintf(buf, 100, lbl, __VA_ARGS__); \
  TIMER(buf, stmt); \
}

#define REAL_TIMER(lbl, stmt) \
{ \
  static time_t start_real, end_real; \
  static long diff_real; \
  time(&start_real); \
  stmt; \
  time(&end_real); \
  diff_real = (long)difftime(end_real, start_real); \
  printf("%16s: %ld seconds (real time)\n", lbl, diff_real); \
}

#define TIMER_ACC_INIT(lbl) \
unsigned long long t##lbl = 0

#define TIMER_ACC(lbl, stmt) \
{ \
  static clock_t t0, t1; \
  static unsigned long long dt; \
  t0 = clock(); \
  stmt; \
  t1 = clock(); \
  dt = (unsigned long long)(t1 - t0) * 1000 / (unsigned long long)CLOCKS_PER_SEC; \
  t##lbl += dt; \
}

#define TIMER_ACC_FINI(lbl) \
  timer_print(g_timer_depth, #lbl, t##lbl);
#endif

void timer_print(int depth, char const * lbl, unsigned long long dt);
