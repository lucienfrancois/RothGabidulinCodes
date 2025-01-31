//=============================================================================
//=============================================================================
//          Fibre decoding algorithm for RothGabidulin codes
//=============================================================================
//=============================================================================



//===================Choice of parameters=======================
q := 3;                                            //===========
n := 5;                                            //==========   
//==============================================================


//===
// Introduction of objects
//===
Fq := GF(q);
Fqn := ext<Fq | n>;
g := Generator(Fqn);
AssignNames(~Fqn,["g"]);
FqnVs, Fqn2V := VectorSpace(Fqn,Fq);
alphazero := NormalElement(Fqn,Fq);
alpha := [alphazero^(q^i) : i in [0..n-1]];
RXY<X,Y> := PolynomialRing(Fqn,2);
RZ<Z> := PolynomialRing(Fqn);


//====================
// Necessary functions
//====================

function FqRank(v)
    //Returns the Fq rank of a vector over Fqn.
    return Dimension(sub<FqnVs | Fqn2V(v) >);
end function;

function DualBasis(beta)
   //Given a Fq-basis beta of Fqn as a sequence, returns the dual basis, as a sequence.
   M := ZeroMatrix(Fqn,n,n);
   for i in [1..n], j in [1..n] do
     M[i,j] := beta[j]^(q^(i-1)) ; 
   end for;
   N := M^(-1);
   return [N[i,1] : i in [1..n]];
end function;
procedure TestDualBasis()
   B := [Random(FqnVs) : i in [1..n]];
   while Dimension(sub<FqnVs | B>) ne n do
      B := [Random(FqnVs) : i in [1..n]];
   end while;
   beta := Inverse(Fqn2V)(B);
   printf "Basis : %o\n",beta;
   betadual := DualBasis(beta);
   printf "Tested dual : %o\n",betadual;
   TrMat := ZeroMatrix(Fqn,n,n);
   for i in [1..n], j in [1..n] do
    TrMat[i,j] := &+ [(beta[i]* betadual[j])^(q^L) : L in [0..n-1]]; 
   end for;
   printf "Trace matrix :\n";
   TrMat;
end procedure;

alphadual := DualBasis(alpha);
alphazerodual := alphadual[1];

function qdeg(a)
   //Returns the qdegree of a non-zero lin poly, and -1 for 0.
   if a eq RZ!0 then
    return -1;
   end if;
   return Floor(Log(q,Degree(a)));
end function;


function RDiv(a,b)
   //RightEuclideanDivision a(Z) = q(b(Z)) + r(Z)
   if qdeg(a) lt qdeg(b) then 
     return [RZ!0,a];
   end if;
   da := qdeg(a);
   db := qdeg(b);
   aa := a - LeadingCoefficient(a) *(LeadingCoefficient(b))^(-q^(da - db)) * (b)^(q^(da - db));
   RES := RDiv(aa,b);
   return [ LeadingCoefficient(a) *(LeadingCoefficient(b))^(-q^(da - db)) * Z^(q^(da - db))+ RES[1],RES[2]];
end function;
procedure TestRDiv();
   a := &+[Random(Fqn)*Z^(q^t) : t in [0..Random([0..12])]];
   b := &+[Random(Fqn)*Z^(q^t) : t in [0..Random([0..8])]];
   RRRR:= RDiv(a,b);
   a;
   b;
   RRRR;
   qq:=RRRR[1];
   rr:=RRRR[2];
   Evaluate(qq,b) + rr eq a;
end procedure;

function LDiv(a,b)
   //LeftEuclideanDivision a(Z) = b(q(Z)) + r(Z)
   if qdeg(a) lt qdeg(b) then 
     return [RZ!0,a];
   end if;
   da := qdeg(a);
   db := qdeg(b);
   aa := a - Evaluate(b, Root(LeadingCoefficient(a) *(LeadingCoefficient(b))^(-1),q^(db)) * (Z)^(q^(da - db)));
   RES := LDiv(aa,b);
   return [ Root(LeadingCoefficient(a) *(LeadingCoefficient(b))^(-1),q^(db)) * Z^(q^(da - db))+ RES[1],RES[2]];
end function;
procedure TestLDiv();
   a := &+[Random(Fqn)*Z^(q^t) : t in [0..Random([0..12])]];
   b := &+[Random(Fqn)*Z^(q^t) : t in [0..Random([0..8])]];
   RRRR:= LDiv(a,b);
   a;
   b;
   RRRR;
   qq:=RRRR[1];
   rr:=RRRR[2];
   Evaluate(b,qq) + rr eq a;
end procedure;

function LEEA(a,b,ds)
   //Linearized Extended Euglidean Algorithm
   r := [a,b];
   i := 1;
   u := [RZ!0,Z];
   v := [Z,RZ!0];
   
   while qdeg(r[i+1]) ge ds do
     //printf "r:\n"; r;
     //printf "u:\n"; u;
     //printf "v:\n"; v;
     RRRR := RDiv(r[i],r[i+1]);
     //RRRR;
     qq := RRRR[1];
     r := r cat [RRRR[2]];
     u := u cat [u[i] - Evaluate(qq,u[i+1])];
     v := v cat [v[i] - Evaluate(qq,v[i+1])];
     i := i+1;
     //printf "===\n";
   end while;
   
   return [r[i+1],u[i+1],v[i+1]];
end function;
procedure TestLEEA();
   a := &+[Random(Fqn)*Z^(q^t) : t in [0..Random([1..10])]];
   b := &+[Random(Fqn)*Z^(q^t) : t in [0..Random([0..5])]];
   while qdeg(a) lt qdeg(b) do
      a := &+[Random(Fqn)*Z^(q^t) : t in [0..Random([1..10])]];
      b := &+[Random(Fqn)*Z^(q^t) : t in [0..Random([0..5])]];
   end while;
   RRRR:= RDiv(a,b);
   a;
   b;
   RRRR;
   RRRR := LEEA(a,b,1);
   RRRR;
   RRRR[1] eq Evaluate(RRRR[3],a) + Evaluate(RRRR[2],b);
end procedure;


function DecodeGabidulinWZAS(r,k)
   //Given r = y + e where f is in Gab[alpha](n,k) 
   // and e of rank at most [(n-k)/2],
   // returns y.
   //[Algorithm WachterZeh-Afanassiev-Sidorenko]
   RR := [ &+[ r[i+1] * alphazerodual^(q^(i+j))  : i in [0..n-1]]  : j in [0..n-1]];
   R := &+ [RR[j] * Z^(q^(j-1)) : j in [1..n]];
   RES := LEEA(Z^(q^n)-Z,R,Floor((n+k)/2));
   rbis := RES[1];
   u := RES[2];
   v := RES[3];
   RES := LDiv(rbis,u);
   if RES[2] eq RZ!0 then
      return [Evaluate(RES[1],alpha[i]) : i in [1..n]];
   else
      //Decoding Failure;
      return -1;
   end if;
end function;
procedure TestDecodeGabidulinWZAS()
   k := Random([1..n-1]);
   tt := Floor((n-k)/2);
   printf "\n\n== Decoding in Gab(%o,%o) with radius %o.==\n",n,k,tt;
   f := &+([Random(Fqn) * Z^(q^j) : j in [0..k-1]] cat [RZ!0]);
   y := [Evaluate(f,alpha[i]) : i in [1..n]]; // equiv to qTransformInv(fcoefs).
   printf "Codeword :\n%o\nFrom:\n%o\n",y,f;

   SpanE := sub<FqnVs | {Random(FqnVs) : t in [1..tt]}>;
   e := [Inverse(Fqn2V)(Random(SpanE)) : i in [1..n]];
   printf "Error :\n%o\n",e;

   r :=[y[i] + e[i] : i in [1..n]];
   printf "Recieved message :\n%o\n",r;
   
   printf "Result :\n";
   DecodeGabidulinWZAS(r,k);
end procedure;

function RandomErrorOfFqRank(tt)
    //returns a random error in Fqn^n of FqRank at most tt.
    SpanE := sub<FqnVs | {Random(FqnVs) : t in [1..tt]}>;
    return [Inverse(Fqn2V)(Random(SpanE)) : i in [1..n]];
end function;

//====================
// DecodingAlgorithms
//====================

function ColumnWiseDecoder(R,mu1,mu2)
   //Given R = M + E where M in C([0..mu1]x[0..mu2]) 
   // and E an error with E[:,i2] of rank at most [(n-mu1-1)/2]
   // returns M.
   M := ZeroMatrix(Fqn,n,n);
   for i2 in [1..n] do
      RES := DecodeGabidulinWZAS([R[i1,i2]:i1 in [1..n]],mu1+1);
      if Type(RES) eq RngIntElt then
         printf "Known decoding failure on column %o\n",i2;
      else 
         for i1 in [1..n] do
            M[i1,i2] := RES[i1];
         end for;
      end if;
   end for;
   return M;
end function;
procedure TestColumnWiseDecoder()
   mu1 := Random([1..n-1]);
   mu2 := Random([1..n-1]);
   tt := Floor((n-mu1-1)/2);
   printf "\n\n== ColumnWiseDecoding in C([0..%o]x[0..%o]) in (Fq^%o)^(ox3) with radius %o.==\n",mu1,mu2,n,tt;

   f := &+([Random(Fqn)*X^(q^i)*Y^(q^j) : i in [0..mu1],j in [0..mu2]] cat [RXY!0]);
   M := ZeroMatrix(Fqn,n,n);
   for i1 in [1..n], i2 in [1..n] do
     M[i1,i2] := Evaluate(f,[alpha[i1],alpha[i2]]);
   end for; 
   printf "  --Codeword:\n%o\n  --From:\n%o\n",M,f;

   E := ZeroMatrix(Fqn,n,n);
   for i2 in [1..n] do 
      e := RandomErrorOfFqRank(tt);
      for i1 in [1..n] do
          E[i1,i2] := e[i1];
      end for;
   end for;
   printf "  --Error:\n%o\n",E;

   R := M + E;
   printf "  --Recieved:\n%o\n",R;
   MM := ColumnWiseDecoder(R,mu1,mu2);
   printf "  --Result:\n%o   ie %o\n",MM,MM eq M;
end procedure;

function FibreWiseDecoder(R,mu1,mu2)
   //Given R = M + E where M in C([0..mu1]x[0..mu2]) 
   // and E an error with E[:,i2] of rank at most [(n-mu1-1)/2] on at least [(n+mu2+1)/2] cols
   // returns M.
   MM := ZeroMatrix(Fqn,n,n);
   M := ZeroMatrix(Fqn,n,n);
   for i2 in [1..n] do
      RES := DecodeGabidulinWZAS([R[i1,i2]:i1 in [1..n]],mu1+1);
      if Type(RES) eq RngIntElt then
         printf "Known decoding failure on column %o\n",i2;
      else 
         for i1 in [1..n] do
            MM[i1,i2] := RES[i1];
         end for;
      end if;
   end for;
   for i1 in [1..n] do
      RES := DecodeGabidulinWZAS([MM[i1,i2]:i2 in [1..n]],mu2+1);
      if Type(RES) eq RngIntElt then
         printf "Known decoding failure on row %o\n",i1;
      else 
         for i2 in [1..n] do
            M[i1,i2] := RES[i2];
         end for;
      end if;
   end for;   
   return M;
end function;
procedure TestFibreWiseDecoder()
   mu1 := Random([1..n-1]);
   mu2 := Random([1..n-1]);
   tt := Floor((n-mu1-1)/2);
   jj := Ceiling((n+mu2+1)/2);
   printf "\n\n== FibreW.D. in C([0..%o]x[0..%o]) in (Fq^%o)^(ox3) with radius %o on %o cols ==\n",mu1,mu2,n,tt,jj;

   f := &+([Random(Fqn)*X^(q^i)*Y^(q^j) : i in [0..mu1],j in [0..mu2]] cat [RXY!0]);
   M := ZeroMatrix(Fqn,n,n);
   for i1 in [1..n], i2 in [1..n] do
     M[i1,i2] := Evaluate(f,[alpha[i1],alpha[i2]]);
   end for; 
   printf "  --Codeword:\n%o\n  --From:\n%o\n",M,f;

   E := ZeroMatrix(Fqn,n,n);
   JJ := Random(Subsets({1..n},jj));
   for i2 in JJ do 
      e := RandomErrorOfFqRank(tt);
      for i1 in [1..n] do
          E[i1,i2] := e[i1];
      end for;
   end for;
   for i1 in [1..n], i2 in {1..n} diff JJ do
       E[i1,i2] := Random(Fqn);
   end for;
   printf "  --Error:\n%o\nJ=%o\n",E,JJ;
   
   R := M + E;
   printf "  --Recieved:\n%o\n",R;
   MM := FibreWiseDecoder(R,mu1,mu2);
   printf "  --Result:\n%o   ie %o\n",MM,MM eq M;
end procedure;

//================
// Examples in paper 
//================
procedure TestColumnWiseDecoderExample1()
   if n ne 5 or q ne 3 then
    error "The parameters are not q=3 and n =5";
   end if;
   mu1 := 2;
   mu2 := 2;
   tt := Floor((n-mu1-1)/2);
   printf "\n\n== ColumnWiseDecoding in C([0..%o]x[0..%o]) in (Fq^%o)^(ox3) with radius %o.==\n",mu1,mu2,n,tt;

   f := RXY!0;
   M := ZeroMatrix(Fqn,n,n);
   printf "  --Codeword:\n%o\n  --From:\n%o\n",M,f;

   E := ZeroMatrix(Fqn,n,n);
   for j in [1..n] do
     E[1,j] := alpha[j];
     for i in [2..n] do
        if j ne 3 and j ne 4 then
          E[i,j] := alpha[j];
        end if;
     end for;
   end for;
   printf "  --Error:\n%o\n",E;

   R := M + E;
   printf "  --Recieved:\n%o\n",R;
   MM := ColumnWiseDecoder(R,mu1,mu2);
   printf "  --Result:\n%o   ie %o\n",MM,MM eq M;
end procedure;
procedure TestColumnWiseDecoderExample2()
   if n ne 5 or q ne 3 then
    error "The parameters are not q=3 and n =5";
   end if;
   mu1 := 2;
   mu2 := 2;
   tt := Floor((n-mu1-1)/2);
   printf "\n\n== ColumnWiseDecoding in C([0..%o]x[0..%o]) in (Fq^%o)^(ox3) with radius %o.==\n",mu1,mu2,n,tt;

   f := RXY!0;
   M := ZeroMatrix(Fqn,n,n);
   printf "  --Codeword:\n%o\n  --From:\n%o\n",M,f;

   E := ZeroMatrix(Fqn,n,n);
   for j in [1..n] do
     E[1,j] := alpha[j];
     for i in [2..n] do
        if j ne 3 and j ne 4 and j ne 5 then
          E[i,j] := alpha[j];
        elif j eq 5 then
          E[i,j] := alpha[6-i];
        end if;
     end for;
   end for;
   printf "  --Error:\n%o\n",E;

   R := M + E;
   printf "  --Recieved:\n%o\n",R;
   MM := ColumnWiseDecoder(R,mu1,mu2);
   printf "  --Result:\n%o   ie %o\n",MM,MM eq M;
end procedure;
procedure TestFibreWiseDecoderExample2()
   if n ne 5 or q ne 3 then
    error "The parameters are not q=3 and n =5";
   end if;
   mu1 := 2;
   mu2 := 2;
   tt := Floor((n-mu1-1)/2);
   printf "\n\n== ColumnWiseDecoding in C([0..%o]x[0..%o]) in (Fq^%o)^(ox3) with radius %o.==\n",mu1,mu2,n,tt;

   f := RXY!0;
   M := ZeroMatrix(Fqn,n,n);
   printf "  --Codeword:\n%o\n  --From:\n%o\n",M,f;

   E := ZeroMatrix(Fqn,n,n);
   for j in [1..n] do
     E[1,j] := alpha[j];
     for i in [2..n] do
        if j ne 3 and j ne 4 and j ne 5 then
          E[i,j] := alpha[j];
        elif j eq 5 then
          E[i,j] := alpha[6-i];
        end if;
     end for;
   end for;
   JJ := {1,2,3,4};
   printf "  --Error:\n%o\nJ=%o\n",E,JJ;
   
   R := M + E;
   printf "  --Recieved:\n%o\n",R;
   MM := FibreWiseDecoder(R,mu1,mu2);
   printf "  --Result:\n%o   ie %o\n",MM,MM eq M;
end procedure;
