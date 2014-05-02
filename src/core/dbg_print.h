#pragma once

#include <stdio.h>

extern int g_print_level;
extern FILE * g_print_stream;

#ifdef NOPRINT
  #define dprintf(...)
  #define dprint2(...)

#else

  #define dprintf(lvl, ...) \
  do { \
    if (!g_print_stream) { g_print_stream = stderr; } \
    if (lvl <= g_print_level) { \
      fprintf(g_print_stream, __VA_ARGS__); \
      fprintf(g_print_stream, "\n"); \
      fflush(g_print_stream); \
    } \
  } while(0)


  #define dprint2(lvl, func, ...) \
  do { \
    if (!g_print_stream) { g_print_stream = stderr; } \
    if (lvl <= g_print_level) { \
      func(g_print_stream,  __VA_ARGS__); \
    } \
  } while(0)

#endif
