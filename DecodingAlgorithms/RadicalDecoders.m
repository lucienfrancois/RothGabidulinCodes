//=============================================================================
//=============================================================================
//          Radical decoding algorithm for RothGabidulin codes
//=============================================================================
//=============================================================================


//===================Choice of parameters=======================
q := 3;                                            //===========
n := 6;                                            //==========   
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
procedure TestFactoringOnTheLeft()
   S := &join[ {[Random([0..n-1]),Random([0..n-1])]} : tries in [1..9]];
   V := &+[Random(Fqn)*Z^(q^i) : i in [0..2]];
   f := &+[Random(Fqn)*X^(q^(s[1]))*Y^(q^(s[2])) : s in S];
   N := Evaluate(V,f);
   fres := FactoringOnTheLeft(N,V,S);
   printf "    V(Z) = %o\n       S = %o\n  f(X,Y) = %o\n  N(X,Y) = %o\n Res(X,Y) = %o\n",V,S,f,N,fres;
   printf "%o\n",f eq fres;
end procedure;


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


function RandomErrorFs3andMinssj(t3,tj)
  //returns a random (nonuniform) error E in Fqn^(nxn) 
  //  with wfs3(E)<= t3
  //  and with minj wssj(E) <= tj
  SpanE := sub<FqnVs| {Random(FqnVs) : tries in [1..t3]} join {FqnVs!0}>;
  SpanSSj := sub<MatSpFq | {Random(MatSpFq) : tries in [1..tj]} join {MatSpFq!0}>;
  realt := Dimension(SpanE);
  BSpanE := [ Inverse(Fqn2V)(SpanE.L) : L in [1..realt]];

  E := ZeroMatrix(Fqn,n,n);
  for i in [1..n] do 
     ROW := Random(SpanSSj);
     for j in [1..n] do
        E[i,j] := &+ ([Fqn!0] cat [ ROW[j,L] * BSpanE[L] : L in [1..realt]]);
     end for;
  end for;

  if Random({true,false}) then
     return Transpose(E);
  end if;
  return E;
end function;





//====================
// DecodingAlgorithms
//====================

function RadicalDecoderFix(R,mu,t)
  //Given R = M + E where M in C([0..mu]²) 
  // and E an error with wfs3(E) <= t and minj wssj(E)<= n-t-mu-2
  // returns M.

  Sf := {[s1 ,s2 ]  : s1 in [0..mu],s2 in [0..mu]};
  SnSet := {[s1 + L,s2 + L]  : s1 in [0..mu],s2 in [0..mu], L in [0..t]};
  Sn := SetToSequence(SnSet);

  //Step1:Constructing the matrix system.
  M := ZeroMatrix(Fqn,n^2,#Sn + t + 1);
  for nbline in [1..n^2] do 
     i := Int2Tpl(nbline); // equation <i[1],i[2]>;
     for dV in [1..t+1] do
       M[nbline,dV] := R[i[1],i[2]]^(q^(dV-1));
     end for;
     for nbS in [1..#Sn] do
       M[nbline,t+1+nbS] := - alpha[i[1]]^(q^(Sn[nbS][1])) * alpha[i[2]]^(q^(Sn[nbS][2]));
     end for;
  end for;

  //Step2: Solve
  if Rank(M) eq #Sn + t + 1 then
     printf "Decoding failure.\n";
     return -1;
  end if;
  K := Kernel(Transpose(M));
  Solution := K.1;

  //Step3: Reconstruct
  V := &+ [Solution[i] * Z^(q^(i-1)) : i in [1..t+1]];
  N := &+ [Solution[t+1+nbS]* X^(q^(Sn[nbS][1])) * Y^(q^(Sn[nbS][2])) : nbS in [1..#Sn]];
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
procedure TestRadicalDecoderFix()
  t := 2;
  mu := Random([1..n-3-t]);
  printf "\n\n== Radical D.Fix in C([0..%o]^2) in (Fq^%o)^(ox3) === \n        [Parameter:%5o]    [Radius:%5o]\n",mu,n,t,n-mu-2-t;
  f := &+[Random(Fqn)*X^(q^(i1))*Y^(q^(i2)) : i1 in [0..mu],i2 in [0..mu]];
  M := ZeroMatrix(Fqn,n,n);
  for i in [1..n],j in [1..n] do
    M[i,j] := Evaluate(f,[alpha[i],alpha[j]]);
  end for;
  printf "  --Codeword:\n%o\n",M;
  E := RandomErrorFs3andMinssj(t,n-mu-t-2);
  printf "  --Error:\n%o\n",E;
  R := M + E;
  printf "  --Recieved:\n%o\n",R;
  MM := RadicalDecoderFix(R,mu,t);
  printf "  --Result:\n%o   ie %o\n",MM,MM eq M; 
end procedure;




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
procedure TestRadicalDecoder()
  mu := Random([1..n-2]);

  printf "\n\n== Radical Decoder in C([0..%o]^2) in (Fq^%o)^(ox3) === \n   [Radius:%5o]\n",mu,n,n-mu-1;
  f := &+[Random(Fqn)*X^(q^(i1))*Y^(q^(i2)) : i1 in [0..mu],i2 in [0..mu]];
  M := ZeroMatrix(Fqn,n,n);
  for i in [1..n],j in [1..n] do
    M[i,j] := Evaluate(f,[alpha[i],alpha[j]]);
  end for;
  printf "  --Codeword:\n%o\n",M;
  t := Random([0..n-2-mu]);
  E := RandomErrorFs3andMinssj(t,n-mu-t-1);
  printf "  --Error:\n%o\n",E;
  R := M + E;
  printf "  --Recieved:\n%o\n",R;
  MM := RadicalDecoder(R,mu);
  printf "  --Result:\n%o   ie %o\n",MM,MM eq M; 
end procedure;



//================
// Examples in paper 
//================
procedure TestRadicalDecoderExample1()
  if n ne 6 or q ne 3 then
    error "The parameters are not q=3 and n =6";
  end if;
  mu := 2;
  printf "\n\n== Radical Decoder in C([0..%o]^2) in (Fq^%o)^(ox3) === \n   [Radius:%5o]\n",mu,n,n-mu-1;
  f := RXY!0;
  M := ZeroMatrix(Fqn,n,n);
  printf "  --Codeword:\n%o\n",M;

  E := ZeroMatrix(Fqn,n,n);
  for j in [1..3] do
    E[1,j] := alpha[1];
    E[2,j] := alpha[1];
    E[3,j] := -alpha[1];
  end for;
  E[1,4] := -alpha[1];
  E[2,4] := alpha[1];
  E[4,4] := -alpha[1];
  E[1,5] := alpha[1];
  E[2,5] := -alpha[1];
  E[4,5] := alpha[1];
  printf "  --Error:\n%o\n",E;
  R := M + E;
  printf "  --Recieved:\n%o\n",R;
  MM := RadicalDecoder(R,mu);
  printf "  --Result:\n%o   ie %o\n",MM,MM eq M; 
end procedure;
procedure TestRadicalDecoderFixExample2()
  if n ne 6 or q ne 3 then
    error "The parameters are not q=3 and n =6";
  end if;
  mu := 2;
  printf "\n\n== Radical Decoder in C([0..%o]^2) in (Fq^%o)^(ox3) === \n   [Radius:%5o]\n",mu,n,n-mu-1;
  f := RXY!0;
  M := ZeroMatrix(Fqn,n,n);
  printf "  --Codeword:\n%o\n",M;

  E := ZeroMatrix(Fqn,n,n);
  for i in [1..3] do
    E[i,1] := alpha[1];
    E[i,2] := alpha[1];
    E[i,3] := alpha[2];
    E[i,4] := alpha[2];
  end for;
  E[4,1] := -alpha[1];
  E[4,2] := -alpha[1];
  E[4,3] := -alpha[2];
  E[4,4] := -alpha[2];
  printf "  --Error:\n%o\n",E;
  R := M + E;
  printf "  --Recieved:\n%o\n",R;
  MM := RadicalDecoderFix(R,mu,1);
  printf "  --Result:\n%o   ie %o\n",MM,MM eq M; 
end procedure;
procedure TestRadicalDecoderExample2()
  if n ne 6 or q ne 3 then
    error "The parameters are not q=3 and n=6";
  end if;
  mu := 2;
  printf "\n\n== Radical Decoder in C([0..%o]^2) in (Fq^%o)^(ox3) === \n   [Radius:%5o]\n",mu,n,n-mu-1;
  f := RXY!0;
  M := ZeroMatrix(Fqn,n,n);
  printf "  --Codeword:\n%o\n",M;

  E := ZeroMatrix(Fqn,n,n);
  for i in [1..3] do
    E[i,1] := alpha[1];
    E[i,2] := alpha[1];
    E[i,3] := alpha[2];
    E[i,4] := alpha[2];
  end for;
  E[4,1] := -alpha[1];
  E[4,2] := -alpha[1];
  E[4,3] := -alpha[2];
  E[4,4] := -alpha[2];
  printf "  --Error:\n%o\n",E;
  R := M + E;
  printf "  --Recieved:\n%o\n",R;
  MM := RadicalDecoder(R,mu);
  printf "  --Result:\n%o   ie %o\n",MM,MM eq M; 
end procedure;
