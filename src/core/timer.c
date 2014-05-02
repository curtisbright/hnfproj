#include "timer.h"

#include <stdio.h>
#include <time.h>

int g_timer_threshold = 0;
int g_timer_depth = 0;
FILE * g_timer_stream;

void timer_print(int depth, char const * lbl, unsigned long long dt)
{
  int i;
  if (depth >= g_timer_threshold) { return; }
  for (i = 0; i < depth; ++i) {
    fprintf(g_timer_stream, "\t");
  }
  fprintf(g_timer_stream, "%16s:", lbl);
  fprintf(g_timer_stream, "\t\t\t");
  fprintf(g_timer_stream, "%llu.%.3llu s\n", dt / 1000, dt % 1000);
  fflush(g_timer_stream);
}
