#include <stdio.h>

#include "cblas.h"
#include "gmp.h"

void ATL_buildinfo(void);

int main(void)
{
  printf("---------- ATLAS ----------\n");
  ATL_buildinfo();
  printf("\n");
  printf("---------- GMP ----------\n");
  printf("%s\n", gmp_version);

  return 0;
}

