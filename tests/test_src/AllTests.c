#include "AllTests.h"

#include <stdio.h>

static void RunAllTests(void)
{
  CuString *output = CuStringNew();
  CuSuite* suite = CuSuiteNew();

  CuSuiteAddSuite(suite, residue_suite());
  CuSuiteAddSuite(suite, rns_suite());

  /*CuSuiteAddSuite(suite, lift_suite());
  CuSuiteAddSuite(suite, applications_suite());*/

  CuSuiteAddSuite(suite, misc_suite());

  CuSuiteAddSuite(suite, mpzMatrix_suite());
  CuSuiteAddSuite(suite, pkMatrix_suite());
  CuSuiteAddSuite(suite, iherm_suite());

  CuSuiteRun(suite);
  CuSuiteSummary(suite, output);
  CuSuiteDetails(suite, output);
  printf("%s\n", output->buffer);
}

#include "timer.h"
#include "dbg_print.h"

int main(void)
{
  g_print_stream = stderr;
  g_print_level = 0;
  g_timer_stream = stderr;
  RunAllTests();
  return 0;
}
