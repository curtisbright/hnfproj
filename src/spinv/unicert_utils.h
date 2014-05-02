#pragma once

#include "mpz_matrix.h"

mpzMatrix_t * unicert_input(long n, long l, int u);
int unicert_check(mpzMatrix_t const * A, int uni);

int hor_check(mpzMatrix_t const * A, mpzMatrix_t const * R);
