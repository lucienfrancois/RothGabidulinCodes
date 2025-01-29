//=============================================================================
//=============================================================================
//          Radical decoding algorithm for RothGabidulin codes
//=============================================================================
//=============================================================================


//===================Choice of parameters=======================
q := 4;                                            //===========
n := 8;                                            //==========
mthd := "Fibre";    //"Fibre" or "Rad"              //========== 
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
MatSpFq := KMatrixSpace(Fq,n,n);
MMM := KMatrixSpace(Fqn,n,n);



//====================
// Necessary functions
//====================

function FqRank(v)
    //Returns the Fq rank of a vector over Fqn.
    return Dimension(sub<FqnVs | Fqn2V(v) >);
end function;

function qdeg(a)
   //Returns the qdegree of a non-zero lin poly, and -1 for 0.
   if a eq RZ!0 then
    return -1;
   end if;
   return Floor(Log(q,Degree(a)));
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

alphadual := DualBasis(alpha);
alphazerodual := alphadual[1];

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


function FactoringOnTheLeft(N,V,S)
  //Given N(X,Y) a bilinearized polynomial,
  // given V(Z) a non-zero linearized polynomial,
  // and S a subset of [0..n-1]x[0..n-1],
  // if Supp(N) is in S + [0..qdegV]
  // and if N is of the form N(X,Y) = V(f(X,Y)),
  // returns V.
  M := Max({s[1] : s in S} join {s[2] : s in S});
  theta := qdeg(V);
  nu := Min({L : L in [0..theta] | Coefficient(V,q^L) ne Fqn!0});
  if N eq RXY!0 then
     return RXY!0;
  end if;
  f := RXY!0;
  for delta in [-M..M] do
     //Getting the coefficients of f on the diagonal X^(q^(delta + tau)) Y^(q^(tau)) for tau = 0..M.
     for tau in [0..M] do
       s := [delta + tau,tau];
       if s in S then
          SIGMA := &+ ([Fqn!0] cat [ Coefficient(V,q^L) * ( MonomialCoefficient(f,X^(q^(s[1] + nu -L))*Y^(q^(s[2] + nu - L))) )^(q^L)  
                                                                        : L in [nu+1..theta] | [s[1] + nu -L,s[2] + nu - L] in S ]);
          c := Root( (Coefficient(V,q^nu))^(-1) * ( MonomialCoefficient(N,X^(q^(s[1] + nu))*Y^(q^(s[2] + nu)))    - SIGMA ) , q^nu);
          f := f + c*X^(q^(s[1]))*Y^(q^(s[2]));
       end if;
     end for;
  end for;
  return f;
end function;

//Bijection [1..n]² -> [1..n²].
Tpl2Int := map< CartesianPower([1..n],2)-> [1..n^2] | x :-> (x[1]-1)*n  + x[2] >;
//Inverse function.
Int2Tpl := map< [1..n^2] -> CartesianPower([1..n],2) | y :-> < (y-1) div n +1 , (y-1) mod n +1> >;
procedure TestTpl2Int()
  for i in [1..n],j in [1..n] do
    printf "Tpl2Int(%2o,%2o) = %3o",i,j,Tpl2Int(i,j);
    printf "      --Int2Tpl-->  %10o   ie %o\n",Int2Tpl(Tpl2Int(i,j)), Int2Tpl(Tpl2Int(i,j)) eq <i,j> ;
  end for;
end procedure;

function RandomErrorTRANK(T)
   //Given an integer T,
   //returns a random (non-uniform) error matrix in Fqn^(n x n) of tensor-rank at most T
   E := MMM!0;
   EE := MMM!0;
   for t in [1..T] do
      u := [Random(Fq) : i in [1..n]];
      v := [Random(Fq) : i in [1..n]];
      w := Random(Fqn);
      for i in [1..n],j in [1..n] do
          EE[i,j] := u[i] * v[j] * w;
      end for;
      E := EE + E;
   end for;
   return E;
end function;


function RadicalDecoder(R,mu)
  //Given R = M + E where M in C([0..mu]²) 
  // and E an error with wSigma(E) <= n-mu-1
  // returns M.

  Sf := {[s1 ,s2 ]  : s1 in [0..mu],s2 in [0..mu]};

  //Find by bisection method the smallest t such that the system has a non-trivial solution.
  ta := -1; //always no non-trivial solution
  tb := n; //always a solution.
  while tb-ta gt 1 do
     ttest := Floor((tb + ta)/2);
     //printf "Step ta = %3o, tb = %3o, ttest = %3o\n",ta,tb,ttest;
     SnSet := {[s1 + L,s2 + L]  : s1 in [0..mu],s2 in [0..mu], L in [0..ttest]};
     Sn := SetToSequence(SnSet);
     M := ZeroMatrix(Fqn,n^2,#Sn + ttest + 1);
     for nbline in [1..n^2] do 
        i := Int2Tpl(nbline); // equation <i[1],i[2]>;
        for dV in [1..ttest+1] do
           M[nbline,dV] := R[i[1],i[2]]^(q^(dV-1));
        end for;
        for nbS in [1..#Sn] do
           M[nbline,ttest+1+nbS] := - alpha[i[1]]^(q^(Sn[nbS][1])) * alpha[i[2]]^(q^(Sn[nbS][2]));
        end for;
     end for;
     if Rank(M) eq #Sn + ttest + 1 then //ie. no nontrivial solution
        ta := ttest;
     else
        tb := ttest;
        Mb := M;
        Snb := Sn;
     end if;
  end while;

  if ta eq n-1 then
     printf "Decoding failure. The error does not have that shape.\n";
     return -1;
  end if;
  //Since ta diffrent than n-1, it means that there exists a system that worked.
  K := Kernel(Transpose(Mb));
  Solution := K.1;

  //Reconstruct
  V := &+ [Solution[i] * Z^(q^(i-1)) : i in [1..tb+1]];
  N := &+ [Solution[tb+1+nbS]* X^(q^(Snb[nbS][1])) * Y^(q^(Snb[nbS][2])) : nbS in [1..#Snb]];
  if V eq RZ!0 then
     printf "Decoding failure.\n";
     return -1;
  end if;
  f := FactoringOnTheLeft(N,V,Sf);
  MM := ZeroMatrix(Fqn,n,n);
  for i in [1..n],j in [1..n] do
    MM[i,j] := Evaluate(f,[alpha[i],alpha[j]]);
  end for;
  return MM;
end function;

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



//====================
// DecodingAlgorithms
//====================

function DecodeTRANKbyProxy(R,mu)
  ///Given R = M + E where M in C([0..mu]²) 
  // and E an error with trank(E)<= 
  // returns M.
  
  if mthd eq "Fibre" then
     return FibreWiseDecoder(R,mu,mu);
  elif mthd eq "Rad" then
     return RadicalDecoder(R,mu);
  end if;
  printf "Wrong method chosen\n";
  return -1;
end function;
procedure TestDecodeTRANKbyProxy()
  mu := Random([1..n-2]);
  tt := Floor((n-mu-1)/2);
  printf "\n\n== TensorRankDecoderByProxy in C([0..%o]^2) in (Fq^%o)^(ox3) === \n   [Radius:%5o]\n",mu,n,tt;
  f := &+[Random(Fqn)*X^(q^(i1))*Y^(q^(i2)) : i1 in [0..mu],i2 in [0..mu]];
  M := ZeroMatrix(Fqn,n,n);
  for i in [1..n],j in [1..n] do
    M[i,j] := Evaluate(f,[alpha[i],alpha[j]]);
  end for;
  printf "  --Codeword:\n%o\n",M;
  E := RandomErrorTRANK(tt);
  printf "  --Error:\n%o\n",E;
  R := M + E;
  printf "  --Recieved:\n%o\n",R;
  MM := DecodeTRANKbyProxy(R,mu);
  printf "  --Result:\n%o   ie %o\n",MM,MM eq M; 
end procedure;

