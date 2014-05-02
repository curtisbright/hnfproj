randomize():
with(LinearAlgebra):
with(ExternalCalling):

n:=8:

A:=LinearAlgebra[RandomMatrix](n,n):
H:=LinearAlgebra[HermiteForm](A):
U1:=Matrix(A.(1/H), order=C_order):
U2:=Matrix(H.(1/A), order=C_order):

unicert := DefineExternal("unicert_maple", "../lib/libhnfproj.so"):

print(A);
Determinant(A), unicert(A);
Determinant(U1), unicert(U1);
Determinant(U2), unicert(U2);

