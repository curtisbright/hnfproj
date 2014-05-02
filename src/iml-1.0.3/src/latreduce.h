/* ---------------------------------------------------------------------
 *
 * -- Integer Matrix Library (IML)
 *    (C) Copyright 2004, 2006 All Rights Reserved
 *
 * -- IML routines -- Version 1.0.1 -- November, 2006
 *
 * Author         : Zhuliang Chen
 * Contributor(s) : Arne Storjohann, Cory Fletcher
 * University of Waterloo -- School of Computer Science
 * Waterloo, Ontario, N2L3G1 Canada
 *
 * ---------------------------------------------------------------------
 *
 * -- Copyright notice and Licensing terms:
 *
 *  Redistribution  and  use in  source and binary forms, with or without
 *  modification, are  permitted provided  that the following  conditions
 *  are met:
 *
 * 1. Redistributions  of  source  code  must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce  the above copyright
 *    notice,  this list of conditions, and the  following disclaimer in
 *    the documentation and/or other materials provided with the distri-
 *    bution.
 * 3. The name of the University,  the IML group,  or the names of its
 *    contributors  may not be used to endorse or promote products deri-
 *    ved from this software without specific written permission.
 *
 * -- Disclaimer:
 *
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,  INDIRECT, INCIDENTAL, SPE-
 * CIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO,  PROCUREMENT  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEO-
 * RY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  (IN-
 * CLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */



#ifndef LATREDUCE_H
#define LATREDUCE_H 1

#include "gmp.h"
#include "common.h"

#define A(i,j) *(AA+(i)*m+j)
#define B(i,j) *(BB+(i)*n+j)
#define d(i)   *(dd+i)

#define	NTMP	20
#define	mpz_next_tmp()	(&mpz_t_tmp[mpz_t_ntmp++])
#define	mpz_free_tmp(n)	(mpz_t_ntmp -= (n))


void LLL(mpz_t *AA, int n, int m);

void ired(mpz_t *A, long n, long m, long n1);

void mpz_initall_tmp(void);

void mpz_freeall_tmp(void);

void mpz_mods(mpz_t r, mpz_t n, mpz_t d);

void mpz_div_round(mpz_t r, mpz_t n, mpz_t d);

void SubtractRow(mpz_t *A, mpz_t *T, long n, long k, long r, mpz_t q);

void SwitchRow(mpz_t *A, mpz_t *T, long n, long k);

void ModSubtractRow(mpz_t *A, mpz_t *T, long n, mpz_t M, mpz_t *dd, long k, \
                    long r, mpz_t q);

void ModSwitchRow(mpz_t *A, mpz_t *T, long n, mpz_t M, mpz_t *dd, long k);

void GetNextU(mpz_t *U[], mpz_t t00, mpz_t t11, mpz_t t12, mpz_t t22);

void TwoReduce(mpz_t *A, mpz_t *T, long n, mpz_t M, mpz_t *dd, long k);

void UpdateMdd(mpz_t M, mpz_t *dd, long n, mpz_t *T);

#endif /* !LATREDUCE_H */
