with(ExternalCalling):

LIBNAME := "libhnfproj.so":

# High order residue
# Input:
#      A - square integer matrix
# Returns:
#      R - high order residue R
#          matrix R with A^{-1} = [*] + A^{-1}R and ||A^{-1}R|| small
#
highOrderResidue := DefineExternal("highOrderResidueMaple", LIBNAME):


# Unimodularity certification
# Input:
#      A - square integer matrix
# Returns:
#      true if det(A) = +/- 1
#      false otherwise
#
unicert := DefineExternal("uniCertMaple", LIBNAME):

