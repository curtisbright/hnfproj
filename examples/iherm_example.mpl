with(LinearAlgebra):
with(ExternalCalling):

randHerm := proc(n)
     local A,f,g,h,i,r1,r2,s;

     A := Matrix(n,n):
     f := rand(1..floor(n/10)):
     g := rand(1..n):
     h := rand(0..1):
     for i to n do A[i,i] := f() od:
     #A[1..floor(n/2),1..floor(n/2)] := RandomMatrix(floor(n/2),floor(n/2)):
     for i to 10*n do
          r1 := g(); r2 := g(); while r2 = r1 do r2 := g() od; s := 1-2*h();
          A[r1,1..-1] := A[r1,1..-1]+s*A[r2,1..-1];
          r1 := g(); r2 := g(); while r2 = r1 do r2 := g() od; s := 1-2*h();
          A[1..-1,r1] := A[1..-1,r1]+s*A[1..-1,r2];
     od:

    return A;

end:


read "../../isdecomp/trimul/trimul.mpl":
QQ:=PackedTriangular:
read "../../isdecomp/trimul/iherm.mpl":

n:=100:

A := randHerm(n);

iherm_native := DefineExternal("iherm_maple", "../lib/libhnfproj.so"):

H1 := iherm_native(A):
printf("\n------------------------------\n");
H2 := myHermiteInstrumented(A,n):

printf("\n------------------------------\n");

print(ArrayTools[IsEqual](H1, H2));

for HH in [H1, H2] do
  for i from 1 to n do
    printf("%d  ", HH[i,i]);
  od;
  printf("\n\n");
od;
