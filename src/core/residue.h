#pragma once

/*#ifdef ATLAS*/

#include "atlas_residue.h"

typedef atlasResidue_t residue_t;

#define residue_init atlas_init
#define residue_fini atlas_fini

#define residue_getEntry atlas_getEntry
#define residue_setEntry atlas_setEntry


#define residue_copy atlas_copy
#define residue_fromMpzMatrix atlas_fromMpzMatrix
#define residue_identity atlas_identity
#define residue_zero atlas_zero

#define residue_isZero atlas_isZero
#define residue_print atlas_print

#define residue_add atlas_add
#define residue_determinant atlas_determinant
#define residue_inverse atlas_inverse
#define residue_gemm atlas_gemm
#define residue_mods atlas_mods
#define residue_quadLift atlas_quadLift
#define residue_scale atlas_scale

/*#endif*/
