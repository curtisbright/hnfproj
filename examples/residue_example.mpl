randomize():
with(LinearAlgebra):
with(ExternalCalling):

highOrderResidue_Reference := proc(A)
     local barA,X,k,n,B0,i,R,RR,M;

     n := RowDimension(A);
     barA := max(op(map(abs,map(op,convert(A,listlist)))));
     X := nextprime(max(10000,ceil(3.61*n^2*barA^1)));
     k := 1;
     while (X^(2^(k+1)-2))^2 * (n^2*barA)^2 < (n^(n-1)*barA^(2*n-2)) do 
          k := k+1 
     od; 
     B0 := LinearAlgebra:-Modular:-Mod(X,A,integer);
     B0 := LinearAlgebra:-Modular:-Inverse(X,B0);
     B0 := Matrix(B0);
     B0 := map(mods,B0,X);
     R := (1/X)*(IdentityMatrix(n)-A.B0);
     for i to k do
          RR := R.R;
          M := map(mods,B0.RR,X);
          R := (1/X)*(RR-A.M);
     od;
     return R;
end:

n:=40:
l:=100:

A := RandomMatrix(n,n, generator=-2^l..2^l);
Ainv := MatrixInverse(A):
C := RandomMatrix(n,n, generator=-2^l..2^l);

highOrderResidue_Native := DefineExternal("highOrderResidue_maple", "../lib/libhnfproj.so"):

R1 := highOrderResidue_Native(A):
R2 := highOrderResidue_Reference(A):

#print(A);
printf("||A||      : %a\n", Norm(A)):
#print(R);
printf("||C||      : %a\n", Norm(C)):
printf("||R1||      : %a\n", Norm(R1)):
printf("||R2||      : %a\n", Norm(R2)):
#print(Ainv.R);
printf("||A^-1.R1||      : %f\n", Norm(Ainv.R1)):
printf("||A^-1.R2||      : %f\n", Norm(Ainv.R2)):
printf("||A^-1.C||      : %f\n", Norm(Ainv.C)):
