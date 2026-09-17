/******************************************************************************/
//
// Copyright (C) 2007 Istvan Montvay
// Copyright (C) 2007 Carsten Urbach
//
// This file is part of tmLQCD.
//
// tmLQCD is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// tmLQCD is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
//
//
//  Computes all roots of a Chebycheff Polynomial of a given 
//  order.
//
//  The function to be approximated is:              1/x^(1/2)
//
//  Set order, approximation inverval and the precision 
//  in chebyRoot.H
//
//  For other powers than -1/2 change the function "func" below
//  appropiately
//
//  Root pairs and square-root of roots are bit reversed ordered.
// 
//  The roots are written to the file roots.dat
//  The monomials are written to Square_root_BR_roots.dat
//  The coefficients are written to coefs.dat
//
//  The normalisation factor is written to normierungLocal.dat
//  it corresponds to the n-th root of C, where n is the order
//  of the polynomial. Note that for the monomial representation
//  you need to take the square-root of that value!
//
//  This runs with CLN, which is available from
//   http://www.ginac.de/CLN/ 
// 
//  Last changed:     Feb. 15, 2007              C. Urbach
// 
//  This is based on a code provided by Istvan Montvay for
//  least squared optimised polynomials, see
//  quadroptRoot.C
//
/******************************************************************************/

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <fstream>
#include <time.h>
#include <sys/stat.h> 
#include <unistd.h>
#include <cmath>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/shm.h>
#include <signal.h>
#include <cln/number.h>
#include <cln/float.h>
#include <cln/complex.h>
#include <cln/real.h>
#include <cln/io.h>
#include <cln/integer_io.h>
#include <cln/float_io.h>
#include <cln/complex_io.h>
#include <cln/univpoly.h>
#include <cln/univpoly_integer.h>
#include <cln/univpoly_complex.h>
#include <cln/univpoly_real.h>


using namespace cln;
using namespace std;

//  Defining input parameters in a header file

#include "chebyRoot.H"


//  Define global variables

cl_F Recb[1+MAXORD-1], Recg[1+MAXORD-2], Orth[1], Coed[1+MAXORD];
//cl_F Coef[1+MAXORD];
cl_F Coef[1+MAXORD], c[1+MAXORD];

/******************************************************************************/
//
//  Calculate the basic integral s_nu for least-square optimization.
//
//  Input parameters:
//  order of the polynomial:                               Maxpow
//
//  the power to be approximated:                          Alpha
//  lower bound of interval of approximation:              Epsilon
//  upper bound of interval of approximation:              Lambda
//
//  The result is put in                                   Sint
//  The precision:                                         Digit

void BaseIntS(int Maxpow, cl_F Alpha, cl_F Epsilon, cl_F Lambda, cl_F* Sint,
              float_format_t Digit)
 {
   int  ord;
   cl_F power, small = As(cl_F)(expt(cl_float(0.1,Digit),DIGIT));

// Loop over powers

    for(ord = 0; ord < 2*Maxpow+1; ord++)
     { power = As(cl_F)(2*Alpha+(ord+1));

        if(abs(power) < small)
        Sint[ord] = As(cl_F)(log(Lambda/Epsilon));

        else 
        Sint[ord] = As(cl_F)(expt(Lambda,power)-expt(Epsilon,power))/power; }

 }
/******************************************************************************/
//
//  Calculate the basic integral t_nu for least-square optimization.
//
//  Input parameters:
//  order of the polynomial:                               Maxpow
//
//  the power to be approximated:                          Alpha
//  lower bound of interval of approximation:              Epsilon
//  upper bound of interval of approximation:              Lambda
//
//  The result is put in                                   Tint
//  The precision:                                         Digit

void BaseIntT(int Maxpow, cl_F Alpha, cl_F Epsilon, cl_F Lambda, cl_F* Tint,
              float_format_t Digit)
 {
   int  ord;
   cl_F power, small = As(cl_F)(expt(cl_float(0.1,Digit),DIGIT));;

// Loop over powers

    for(ord = 0; ord < Maxpow+1; ord++)
     { power = As(cl_F)(Alpha+(ord+1));

        if(abs(power) < small)
        Tint[ord] = As(cl_F)(log(Lambda/Epsilon));

        else
        Tint[ord] = As(cl_F)(expt(Lambda,power)-expt(Epsilon,power))/power; }

 }
/******************************************************************************/
//
//  Evaluate the approximate polynomial up to the order        Maxpow
//  at the variable value                                      Xval
//
//  The recurrence coefficients are in                         Recb
//  and in                                                     Recg
//  Constant term of the first orthogonal polynomial:          Orth
//  Expansion coefficients in orthogonal polynomials:          Coed

cl_F Recurev(int Maxpow, cl_F Xval,
             cl_F* Recb,cl_F* Recg,cl_F* Orth,cl_F* Coed)
 {
   int  ord;
   cl_F res, orth, orthb, orthg;


// Check input

    if(Maxpow < 0 || Maxpow > MAXORD)
     { cout <<endl <<"wrong order in Recurev:" <<endl;
       cout <<endl <<Maxpow <<endl; }

// Start iteration

   orth = As(cl_F)(ONE);
   res = Coed[0]*orth;

    if(Maxpow > 0)
     { orthb = orth;
       orth = Xval+Orth[0];
       res = res+Coed[1]*orth; }

// Iteration for recurrence

    for(ord = 2; ord < Maxpow+1; ord++)
     { orthg = orthb;
       orthb = orth;
       orth = (Xval+Recb[ord-1])*orthb+Recg[ord-2]*orthg;
       res = res+Coed[ord]*orth; }

   return res;
 }
/******************************************************************************/
//
// Write elements of a list of numbers in the array       Wlist
// of length                                              Leng
// as assignment statements  
// to the elements of an array                            Arrayname
//
// The ofstream of the output file is                     Ostream
// The format of assignements is C, F(ortran) or T(ao):   Format
//  
// At the beginning write the string                      Text

void WriteAssign(ostream& Ostream, char* Format, char* Text,
                 cl_F* Wlist, int Leng, char* Arrayname)
 {
   int   ord;
   char c[] = "C", fortran[] = "Fortran", tao[] = "Tao";

   Ostream <<endl <<Text <<endl <<endl;

   Ostream.setf(ios::scientific);
   Ostream.precision(16);

    for(ord = 0; ord < Leng; ord++)
     {  if(Format[0] == c[0] || Format[0] == tao[0])
        Ostream <<"   " <<Arrayname <<"[" <<ord <<"] = "
                <<double_approx(Wlist[ord]) <<endl;

        if(Format[0] == fortran[0])
        Ostream <<"      " <<Arrayname <<"(" <<ord+1 <<") = "
                <<double_approx(Wlist[ord]) <<endl; }

 }

void WriteRoots(char * filename, cl_N *list, const int length) {
  ofstream ofs(filename, ios::out);
  ofs.setf(ios::scientific);
  ofs.precision(16);
  ofs << "Nr.           Re            Im" << endl;
  for(int i = 0; i < length; i++) {
    ofs << i << " " <<
      double_approx(realpart(list[i])) << " " <<
      double_approx(imagpart(list[i])) << endl;
  }
  ofs.close();
}

void WriteCoefs(cl_N *list, const int length) {
  ofstream ofs("coefs.dat", ios::out);
  ofs.setf(ios::scientific);
  ofs.precision(16);
  for(int i = 0; i < length; i++) {
    ofs << i << " " <<
      double_approx(realpart(list[i])) << endl;
  }
  ofs.close();
}


/******************************************************************************/
//
//  Calculating the polynomial approximation minimizing the integral
//  of quadratic deviation from x^(-Alpha) in an interval.
//
//  The weight function in least-squares optimization
//  is the inverse of the function.
//
//  Recursion relations for orthogonal polynomials are used
//  whenever possible.
//
//  Input parameters:
//  order of the optimized polynomial:                     Maxpow  
//  the power to be approximated:                          Alpha
//  lower bound of interval of approximation:              Epsilon
//  upper bound of interval of approximation:              Lambda
//  the precision:                                         Digit
//
//  name of the file for writing results:                  Filename
//  format for assignments, C, Fortran or Tao:             Format
//  Print intermediate results for                         Printing=TRUE

void Quadropt(int Maxpow, cl_F Alpha, cl_F Epsilon, cl_F Lambda,
              float_format_t Digit,
              char* Filename, char* Format, int Printing)
 {
   int  ord, or1, or2, leng;
   cl_F fmu, quad, scal, p1e, p1l;

   cl_F sint[1+2*Maxpow],
        norq[1+Maxpow],
        bint[1+Maxpow],
        tint[1+Maxpow],

        rmup[1+2*Maxpow-2],
        rmuv[1+2*Maxpow-1],
        rmum[1+2*Maxpow-2],

        bmup[1+Maxpow-2],
        bmuv[1+Maxpow-1],
        bmum[1+Maxpow-2];


// Define necessary basic integrals of the weight function

   BaseIntS(Maxpow, Alpha, Epsilon, Lambda, sint, Digit);

// Start iteration for determining the orthogonal polynomials

    if(Printing == TRUE)
    cout <<endl <<endl <<"recurrence constants for orthogonal polynomials:" 
         <<endl;

   norq[0] = sint[0];

    if(Maxpow > 0)
     { Recb[0] = -sint[1]/sint[0];
       Orth[0] = -sint[1]/sint[0];
       norq[1] = sint[2]-expt(sint[1],2)/sint[0]; }

    if(Maxpow > 1)
     {  for(ord = 0; ord < 2*Maxpow-1; ord++)
        rmum[ord] = sint[ord];

        for(ord = 0; ord < 2*Maxpow; ord++)
        rmuv[ord] = sint[ord+1]+Orth[0]*sint[ord];

       fmu = Orth[0];
       Recb[1] = -fmu-rmuv[2]/rmuv[1];
       Recg[0] = -rmuv[1]/rmum[0]; }

// Iteration for orthogonal polynomials

     for(or1 = 1; or1 < Maxpow; or1++)
      { Recb[or1] = -fmu-rmuv[or1+1]/rmuv[or1];
        Recg[or1-1] = -rmuv[or1]/rmum[or1-1];
        fmu = fmu+Recb[or1];

         for(or2 = 0; or2 < 2*Maxpow-or1; or2++)
         rmup[or2] = rmuv[or2+1]+Recb[or1]*rmuv[or2]+Recg[or1-1]*rmum[or2];

        norq[or1+1] = rmup[or1+1];

         for(or2 = 0; or2 < 2*Maxpow-or1-1; or2++)
         rmum[or2] = rmuv[or2];

         for(or2 = 0; or2 < 2*Maxpow-or1; or2++)
         rmuv[or2] = rmup[or2];

         if(Printing == TRUE)
         cout <<endl <<or1 <<endl; }

// Calculate orthogonal expansion coefficients of optimized polynomial
// Calculate basic integrals

   BaseIntT(Maxpow, Alpha, Epsilon, Lambda, tint, Digit);

// Start iteration

    if(Printing == TRUE)
    cout <<endl <<endl <<"expansion coefficients of optimized polynomials:" 
         <<endl;

   bint[0] = tint[0];
   Coed[0] = bint[0]/norq[0];

    if(Maxpow > 0)
     { bint[1] = tint[1]+Orth[0]*tint[0];
       Coed[1] = bint[1]/norq[1]; }

    if(Maxpow > 1)
     {  for(ord = 0; ord < Maxpow-1; ord++)
        bmum[ord] = tint[ord];

        for(ord = 0; ord < Maxpow; ord++)
        bmuv[ord] = tint[ord+1]+Orth[0]*tint[ord]; }

// Perform iteration

    for(or1 = 1; or1 < Maxpow; or1++)
     {
        for(or2 = 0; or2 < Maxpow-or1; or2++)
        bmup[or2] = bmuv[or2+1]+Recb[or1]*bmuv[or2]+Recg[or1-1]*bmum[or2];

       bint[or1+1] = bmup[0];
       Coed[or1+1] = bint[or1+1]/norq[or1+1];

        for(or2 = 0; or2 < Maxpow-or1-1; or2++)
        bmum[or2] = bmuv[or2];

        for(or2 = 0; or2 < Maxpow-or1; or2++)
        bmuv[or2] = bmup[or2];

         if(Printing == TRUE)
         cout <<endl <<or1 <<endl; }

// Print out numerical value at the minimum found

    if(Printing == TRUE)
    cout <<endl <<endl <<"minimum calculated with Digits [" 
         <<Digit <<"]" <<endl;

   quad = As(cl_F)(ONE);
   scal = Lambda-Epsilon;

    for(ord = 0; ord < Maxpow+1; ord++)
    quad = quad-Coed[ord]*bint[ord]/scal;

    if(Printing == TRUE)
     { cout <<endl <<"quadratic deviation: " <<quad <<endl;
       cout <<endl <<"    on interval of weight  " <<scal <<endl; }

// Check interval ends

  p1e = Recurev(Maxpow,Epsilon,Recb,Recg,Orth,Coed);
  p1e = p1e*As(cl_F)(expt(Epsilon,Alpha))-As(cl_F)(ONE);

  p1l = Recurev(Maxpow,Lambda,Recb,Recg,Orth,Coed);
  p1l = p1l*As(cl_F)(expt(Lambda,Alpha))-As(cl_F)(ONE);

    if(Printing == TRUE)
     { cout <<endl <<"relative deviation of polynomial at interval ends:"
            <<endl;
       cout <<endl <<p1e <<endl;
       cout <<endl <<p1l <<endl; }

// Write coefficients of the polynomial to file

   ofstream Ostream(Filename, ios::out);

   Ostream <<endl <<" Maxpow, Alpha, Epsilon, Lambda:" <<endl;
   Ostream <<endl <<"   " <<Maxpow <<"   " <<Alpha 
            <<"   " <<Epsilon <<"   " <<Lambda <<endl;
   Ostream <<endl <<" calculated with Digits [" <<Digit <<"]" <<endl;

   Ostream <<endl <<" Quadratic deviation: " <<quad <<endl;
   Ostream <<endl <<" on interval of weight  " <<scal <<endl;

   Ostream <<endl <<" Relative deviation of polynomial at interval ends:"
                <<endl;
   Ostream <<endl <<p1e <<endl;
   Ostream <<endl <<p1l <<endl <<endl;

   leng = 1;
   WriteAssign(Ostream,Format," Result for starting recurrence coefficient:",
               Orth,leng,"Orth1");

   leng = Maxpow+1;
   WriteAssign(Ostream,Format," Result for expansion coefficients:",
               Coed,leng,"Coed1");

   leng = Maxpow;
   WriteAssign(Ostream,Format," Result for first coefficients for recurrence:",
               Recb,leng,"Recb1");

   leng = Maxpow-1;
   WriteAssign(Ostream,Format," Result for second coefficients for recurrence:",
               Recg,leng,"Recg1");

   cout <<endl <<"polynomial is ready" <<endl;

 }
/******************************************************************************/
//
//  Restore polynomial from its expansion in orthogonal polynomials.
//
//  The input coefficients will be left unchanged.
//
//  Input parameters:
//
//  Maximal order of the polynomial:                          Maxpow
//
//  The recurrence coefficients are in                        Recb
//  and in                                                    Recg
//  Constant term of the first orthogonal polynomial:         Orth
//  Expansion coefficients in orthogonal polynomials:         Coed
//
//  the obtained coefficients of the polynomial:              Coef

void RestorePolyCoef(int Maxpow, cl_F* Recb, cl_F* Recg, cl_F* Orth, cl_F* Coed,
                     cl_F* Coef)
 {
   int ord, orv;

   cl_F poly[1+Maxpow];
   cl_F opl[1+Maxpow], oplm[1+Maxpow], oplp[1+Maxpow];


// Check input

    if(Maxpow < 0 || Maxpow > MAXORD)
     { cout <<endl <<"wrong order in RestorePolyCoef:" <<endl;
       cout <<endl <<Maxpow <<endl; }

// Start iteration

    for(ord = 0; ord < Maxpow+1; ord++)
     { poly[ord] = As(cl_F)(ZERO);
       opl[ord]  = As(cl_F)(ZERO);
       oplm[ord] = As(cl_F)(ZERO);
       oplp[ord] = As(cl_F)(ZERO); }

   oplm[0] = As(cl_F)(ONE);

   opl[0] = Orth[0];
   opl[1] = As(cl_F)(ONE);

   poly[0] = Coed[0]*oplm[0]+Coed[1]*opl[0];
   poly[1] = Coed[0]*oplm[1]+Coed[1]*opl[1];

// Loop over order of orthogonal polynomial

    for(ord = 1; ord < Maxpow; ord++)
     {
        for(orv = 0; orv < Maxpow+1; orv++)
         { oplp[orv] = Recb[ord]*opl[orv]+Recg[ord-1]*oplm[orv];
           if(orv > 0) oplp[orv] = oplp[orv]+opl[orv-1]; }

        for(orv = 0; orv < Maxpow+1; orv++)
        poly[orv] = poly[orv]+Coed[ord+1]*oplp[orv];

        for(orv = 0; orv < Maxpow+1; orv++)
        oplm[orv] = opl[orv];

        for(orv = 0; orv < Maxpow+1; orv++)
        opl[orv] = oplp[orv]; }

// Extract coefficients

    for(ord = 0; ord < Maxpow+1; ord++)
    Coef[ord] = As(cl_F)(poly[Maxpow-ord]);

 }
/******************************************************************************/
//
//  Evaluate a complex polynomial
//
//  The polynomial is                                    Poly
//  The order is                                         Maxpow
//
//  The value of the variable:                           Valu

cl_N EvalPoly(cl_N* Poly, int Maxpow, cl_N Valu)
 {
   int  pow;
   cl_N xpow, sum;


   sum = As(cl_N)(complex(ZERO,ZERO));
   xpow = As(cl_N)(complex(ONE,ZERO));

    for(pow = 0; pow < Maxpow+1; pow++)
     { sum = sum+xpow*Poly[pow];
       xpow = xpow*Valu; }

   return sum;

 }

// this routine evaluated the product representation of a
// polynomial with roots roots, order Maxpow at value x and
// with (overall) normalisation factor norma

cl_N EvalPolyProd(cl_N* roots, const int Maxpow, 
		  const cl_N x, const cl_N norma = complex(ONE, ZERO)) {

  cl_N prod = As(cl_N)(norma);
  for(int i = 0; i < Maxpow; i++) {
    prod = prod * (x-roots[i]);
  }
  return(prod);
}

/******************************************************************************/
//
//  Find a root of a complex polynomial by Laguerre iteration.
//
//  The polynomial is                                    Poly
//  The order is                                         Maxpow
//
//  The precision:                                       Digit
//
//  Print intermediate results for                       Printing=TRUE

cl_N Lasolv(cl_N* Poly, int Maxpow, float_format_t Digit, int Printing, 
	    cl_N root = complex(ZERO,ZERO), const int itemax=100)
 {
   int  pow, ite;

   cl_F angl, small = As(cl_F)(expt(cl_float(0.1,Digit),DIGIT/2));

   cl_N dif1[Maxpow], dif2[Maxpow-1];
   cl_N val0, val, val1, val2, denp, denm, las1, las2, sqrv;
   //   cl_N root;
    for(pow = 0; pow < Maxpow; pow++)
    dif1[pow] = (pow+1)*Poly[pow+1];

    for(pow = 0; pow < Maxpow-1; pow++)
    dif2[pow] = (pow+1)*dif1[pow+1];

// The maximal allowed number of iterations is set here;
// this can be chosen larger, but 100 usually suffices

//   root = As(cl_N)(complex(ZERO,ZERO));
   val0 = EvalPoly(Poly,Maxpow,root);

// Iteration

    for(ite = 0; ite < itemax; ite++)
     { 
       val = val0;
       val1 = EvalPoly(dif1,Maxpow-1,root);
       val2 = EvalPoly(dif2,Maxpow-2,root);

       sqrv = (Maxpow-1)*((Maxpow-1)*val1*val1-Maxpow*val0*val2);
       angl = HALF*cl_float(phase(sqrv),Digit);
       sqrv = sqrt(abs(sqrv))*complex(cos(angl),sin(angl));
       denp = val1+sqrv;
       denm = val1-sqrv;

        if(denp == complex(ZERO,ZERO))
        root = root-Maxpow*val0/denm;

        else
         {  if(denm == complex(ZERO,ZERO))
            root = root-Maxpow*val0/denp;

            else
             { las1 = -Maxpow*val0/denp;
               las2 = -Maxpow*val0/denm;

                if(realpart(las1*conjugate(las1)) <
                   realpart(las2*conjugate(las2)))
                root = root+las1;

                else
                root = root+las2; } }

//  Look whether the root is good enough

       val0 = EvalPoly(Poly,Maxpow,root);

        if(abs(val0) == ZERO ||
          (abs(val0) < small) && abs(val0/val) > 0.7)
         {
            if(Printing == TRUE)
             { cout << endl << "Laguerre iterations: " << ite << endl;
               cout << endl << "root = " << root << endl;
               cout << endl << "value at root: " << val0 << endl; }

           break; } }

    if(ite >= itemax) {
      cout <<endl << "Laguerre iteration did not converge" <<endl;
      exit(5);
    }

   return root;

 }
/******************************************************************************/
//
//  Find the complex roots of a complex polynomial       Poly
//  The order is                                         Maxpow
//  The result will be in the array                      Root
//
//  The precision:                                       Digit
//
//  Print intermediate results for                       Printing=TRUE

void Polyrootc(cl_N* Poly, int Maxpow, cl_N* Root,
               float_format_t Digit, int Printing) {
  
  int  ord, pow, fnd, pov, maxp;

  cl_N poly[1+Maxpow], polc[1+Maxpow], coef[1+Maxpow], coen[1+Maxpow];


  // Put coefficients in an array

  for(pow = 0; pow < Maxpow+1; pow++) { 
    poly[pow] = As(cl_N)(Poly[pow]);
    coef[pow] = As(cl_N)(Poly[pow]); 
  }

  for(pow = 0; pow < Maxpow+1; pow++) {
    polc[pow] = As(cl_N)(complex(ZERO,ZERO));
  }

  polc[0] = As(cl_N)(complex(ONE,ZERO));
  fnd = -1;

  // Loop for finding all roots

  for(ord = 0; ord < Maxpow; ord++) {
    fnd++;
    pov = Maxpow-fnd;

    if(fnd < Maxpow) {
      if(Printing == TRUE) {
	cout <<endl <<" root number: " <<fnd+1 <<endl;
      }

      if((ord%2 == 1) && (Maxpow%2 == 0) && (false)) {
	//	Root[fnd] = Lasolv(poly,pov,Digit,Printing, Root[fnd-1]);
	Root[fnd] = As(cl_N)(conjugate(Root[fnd-1]));
 	cl_N val0 = EvalPoly(poly, Maxpow, Root[fnd]);
	
 	if(Printing == TRUE) {
 	  cout << endl << "root = " << Root[fnd] << endl;
 	  cout << endl << "value at root: " << val0 << endl; 
 	}
      }
      else {
	Root[fnd] = Lasolv(poly,pov,Digit,Printing, complex(ONE,ONE), 500);
      }

      for(pow = Maxpow; pow > 0; pow--) {
	polc[pow] = polc[pow-1]-Root[fnd]*polc[pow];
      }

      polc[0] = -Root[fnd]*polc[pow];

      // Divide the polynomial by the root

      maxp = Maxpow-fnd-1;
      coen[maxp] = coef[maxp+1];

      for(pow = maxp-1; pow > -1; pow--) {
	coen[pow] = coef[pow+1]+Root[fnd]*coen[pow+1];
      }

      for(pow = 0; pow < maxp+1; pow++) {
	coef[pow] = coen[pow];
	poly[pow] = coef[pow]; 
      } 
    }

    else {
      break;
    }
  }

// Compare input with product of root factors

  if(Printing == TRUE) {
    for(pow = 0; pow < Maxpow+1; pow++) {
      polc[pow] = Poly[pow]-poly[0]*polc[pow];
    }

    cout <<endl <<"control polynomial should be close to zero:" <<endl;

    for(pow = 0; pow < Maxpow+1; pow++) {
      cout <<endl <<"  x^{" <<pow <<"}" <<endl;
          cout <<endl <<polc[pow] <<endl; 
    } 
  }
}

/******************************************************************************/
//
// Find the optimal order of the roots in             Root
// for calculating the polynomial with order          Maxpow
// by multiplication of the root factors.
//
// The overall factor in the polynomial is:           Coef0
// The function to be approximated is                 x^(-Alpha)
// Optimization is done in the interval               [Epsilon,Lambda]
//
// Check values in                                    Ncheck
// equidistant points of the interval.
//
// Print intermediate results for                     Printing=TRUE
//
// This calculation is done with precision            digit

void Optimord(int Maxpow, cl_N* Root, cl_F Coef0,
              cl_F Alpha, cl_F Epsilon, cl_F Lambda,
              int Ncheck, int Printing)
 {
   int nmax = Ncheck, digval = 40;
   float_format_t digit = float_format(digval);

   int nl, rr, rt, rv, rc, facc[Maxpow];

// Check input

    if(Ncheck < 3)
     { cout <<endl <<"number of checkpoints in Optimord should be at least 3"
            <<endl;
       nmax = 3; }

    if(Ncheck > 47)
     { cout <<endl <<"with this value of Ncheck it would take somewhat longer"
            <<endl;
       cout <<endl <<"the number of checkpoints will be set to 47" <<endl;
       nmax = 47; }

   cl_F large = As(cl_F)(expt(cl_float(10.0,digit),digval));
   cl_F small = As(cl_F)(expt(cl_float( 0.1,digit),digval));

   cl_F xx[nmax], cc[nmax], values[Maxpow][nmax], mx[Maxpow];
   cl_F max, min, mn;

   cl_N root[Maxpow], ff[Maxpow][nmax], pp[nmax];


// Define the points where comparison will be made

    for(nl = 0; nl < nmax; nl++)
    xx[nl] = cl_float(Epsilon*(nmax-nl-1)/(nmax-1)+Lambda*nl/(nmax-1),digit);

// Special case: Epsilon=0

    if(xx[0] == 0.0)
    xx[0] = xx[1]/100;

// The values to be compared

    for(nl = 0; nl < nmax; nl++)
    cc[nl] = cl_float(Coef0*exp(ln(xx[nl])*Alpha),digit);

// Loop over the root factors of the polynomial

    for(rr = 0; rr < Maxpow; rr++)
    root[rr] = complex(cl_float(realpart(Root[rr]),digit),
                       cl_float(imagpart(Root[rr]),digit));

    for(rr = 0; rr < Maxpow; rr++)
     { facc[rr] = 1;

        for(nl = 0; nl < nmax; nl++)
        ff[rr][nl] = xx[nl]-root[rr]; }

// Find optimal order by minimizing maximal ratio in partial products

    for(nl = 0; nl < nmax; nl++)
    pp[nl] = complex(cl_float(1,digit),cl_float(0,digit));

    if(Printing == TRUE)
     { cout <<endl <<endl <<"  The logarithm of the overall constant:" <<endl;
       cout <<endl <<cl_float(ln(abs(Coef0)),digit) <<endl;
       cout <<endl <<"  Maximal logarithm of ratios in partial products:"
            <<endl; }

// Loop over the roots

   rc = 0;

    for(rt = 0; rt < Maxpow; rt++)
     {
        for(rr = 0; rr < Maxpow; rr++)
         if(facc[rr] == 1)
          { max = small;
            min = large;

             for(nl = 0; nl < nmax; nl++)
              { values[rr][nl] =
                 cl_float(ln(abs(cc[nl]*pp[nl]*ff[rr][nl])),digit);

                if(max < values[rr][nl]) max = values[rr][nl];
                if(min > values[rr][nl]) min = values[rr][nl]; }

            mx[rr] = max-min; }

       mn = large;
       rv = -1;

        for(rr = 0; rr < Maxpow; rr++)
         if(facc[rr] == 1)
          if(mn > mx[rr])
          { mn = mx[rr];
            rv = rr; }

        if(Printing == TRUE)
        cout <<endl <<rt+1 <<":  " <<mn <<endl;

       Root[rc] = root[rv];
       rc++;
       facc[rv] = 0;

        for(nl = 0; nl < nmax; nl++)
        pp[nl] = pp[nl]*ff[rv][nl]; }

    if(Printing == TRUE)
     { cout <<endl <<endl <<"  Product check of ratios:" <<endl;

        for(nl = 0; nl < nmax; nl++)
         { cout <<endl <<"  x = " <<xx[nl] <<" :" <<endl;
           cout <<endl <<cc[nl]*pp[nl] <<endl; } }

 }
/******************************************************************************/
//
// Find the optimal order of the roots in             Root
// for calculating the polynomial with order          Maxpow
// by multiplication of the root factors.
//
// Complex conjugate pairs of roots are kept together:
// positive imaginary part first in pair.
//
// The overall factor in the polynomial is:           Coef0
// The function to be approximated is                 x^(-Alpha)
// Optimization is done in the interval               [Epsilon,Lambda]
//
// Check values in                                    Ncheck
// equidistant points of the interval.
//
// Print intermediate results for                     Printing=TRUE
//
// This calculation is done with precision            digit

void OptimordPair(int Maxpow, cl_N* Root, cl_F Coef0,
                  cl_F Alpha, cl_F Epsilon, cl_F Lambda,
                  int Ncheck, int Printing)
{
  int nmax = Ncheck, digval = 40;
  float_format_t digit = float_format(digval);

  int nl, rv, rn, r1, r2, rc, rt, rr, maxp;
  int facp[Maxpow], facr[Maxpow], facl[Maxpow];

  // Check input

  if(Ncheck < 3)
    { cout <<endl <<"number of checkpoints in Optimord should be at least 3"
	   <<endl;
      nmax = 3; }

  if(Ncheck > 137)
    { cout <<endl <<"with this value of Ncheck it would take somewhat longer"
	   <<endl;
      cout <<endl <<"the number of checkpoints will be set to 137" <<endl;
      nmax = 137; }

  cl_F large = As(cl_F)(expt(cl_float(10.0,digit),digval));
  cl_F small = As(cl_F)(expt(cl_float( 0.1,digit),digval));

  cl_F xx[nmax], cc[nmax], values[Maxpow][nmax], mx[Maxpow];
  cl_F max, min, mn;

  cl_N rot[Maxpow], root[Maxpow], ff[Maxpow][nmax], pp[nmax], prd;


  // Define the points where comparison will be made

  for(nl = 0; nl < nmax; nl++)
    xx[nl] = cl_float(Epsilon*(nmax-nl-1)/(nmax-1)+Lambda*nl/(nmax-1),digit);

  // Special case: Epsilon=0

  if(xx[0] == 0.0)
    xx[0] = xx[1]/100;

  // The values to be compared

  for(nl = 0; nl < nmax; nl++)
    cc[nl] = cl_float(Coef0*exp(ln(xx[nl])*Alpha),digit);

  for(rr = 0; rr < Maxpow; rr++)
    rot[rr] = complex(cl_float(realpart(Root[rr]),digit),
                      cl_float(imagpart(Root[rr]),digit));

  // Pair up indices of complex cojugate roots: 
  // positive imaginary part first in pair
  // Store indices of pairs and real roots

  for(rr = 0; rr < Maxpow; rr++)
    { facr[rr] = 1;
      facl[rr] = 0; }

  rn = 0;

  for(rr = 0; rr < Maxpow; rr++)
    if(abs(imagpart(rot[rr])) < 100*small)
      { facp[rn] = -rr;
        facr[rr] = 0;
        rn++; }

  // Check whether the complex roots are in pairs

  if((Maxpow-rn)%2 != 0)
    { cout <<endl <<"The number of complex roots is not even: " 
	   <<Maxpow-rn <<endl;
      return; }

  for(r1 = 0; r1 < Maxpow; r1++)
    if(facr[r1] == 1)
      {  for(r2 = r1+1; r2 < Maxpow; r2++)
          if(facr[r2] == 1)
	    {
              if(abs(conjugate(rot[r2])-rot[r1]) < 100*small)
		{
                  if(imagpart(rot[r1]) > 0)
		    { facp[rn] = r1;
		      facr[r1] = 0;
		      rn++;
		      facp[rn] = r2;
		      facr[r2] = 0;
		      rn++; }
                  else
		    { facp[rn] = r2;
		      facr[r2] = 0;
		      rn++;
		      facp[rn] = r1;
		      facr[r1] = 0;
		      rn++; }

		  break; } } }

  // Check whether the complex roots are in complex conjugate pairs

  if((Maxpow-rn) != 0)
    { cout <<endl <<"The roots are not in complex conjugate pairs: " 
	   <<Maxpow-rn <<endl;
      return; }

  // Reorder roots: first single real then complex pairs

  rc = 0;

  for(rn = 0; rn < Maxpow; rn++)
    {  if(facp[rn] < 0)
	{ root[rc] = rot[-facp[rn]];
	  facl[rc] = -1;
	  rc++; }
      else
	{ root[rc] = rot[facp[rn]];
	  facl[rc] = (imagpart(root[rc]) > 0);
	  rc++; } }

  // Calculate root factors

  for(rr = 0; rr < Maxpow; rr++)
    for(nl = 0; nl < nmax; nl++)
      ff[rr][nl] = xx[nl]-root[rr];

  for(nl = 0; nl < nmax; nl++)
    pp[nl] = complex(cl_float(1,digit),cl_float(0,digit));

  if(Printing == TRUE)
    { cout <<endl <<endl <<"  The logarithm of the overall constant:" <<endl;
      cout <<endl <<cl_float(ln(abs(Coef0)),digit) <<endl;
      cout <<endl <<"  Maximal logarithm of ratios in partial products:"
	   <<endl; }

  // Find optimal order by minimizing maximal ratio in partial products
  // Loop over complex conjugate root pairs and real roots

  maxp = 0;

  for(rr = 0; rr < Maxpow; rr++)
    { facr[rr] = 1;
      if(facl[rr] != 0) maxp++; }

  rc = 0;

  for(rt = 0; rt < maxp; rt++)
    {
      for(rr = 0; rr < Maxpow; rr++)
	if(facl[rr] != 0 && facr[rr] == 1)
          { max = small;
            min = large;

	    for(nl = 0; nl < nmax; nl++)
              { prd = ff[rr][nl];
                if(facl[rr] == 1) prd = prd*ff[rr+1][nl];

                values[rr][nl] =
		  cl_float(ln(abs(cc[nl]*pp[nl]*prd)),digit);

                if(max < values[rr][nl]) max = values[rr][nl];
                if(min > values[rr][nl]) min = values[rr][nl]; }

            mx[rr] = max-min; }

      mn = large;
      rv = -1;

      for(rr = 0; rr < Maxpow; rr++)
	if(facl[rr] != 0 && facr[rr] == 1)
          if(mn > mx[rr])
	    { mn = mx[rr];
	      rv = rr; }

      if(Printing == TRUE)
        cout <<endl <<rt+1 <<":  " <<mn <<endl;

      Root[rc] = root[rv];
      rc++;
      facr[rv] = 0;

      for(nl = 0; nl < nmax; nl++)
        pp[nl] = pp[nl]*ff[rv][nl]; 

      if(facl[rv] == 1)
	{ rv++;
	  Root[rc] = root[rv];
	  rc++;
	  facr[rv] = 0;

	  for(nl = 0; nl < nmax; nl++)
            pp[nl] = pp[nl]*ff[rv][nl]; } }

  if(Printing == TRUE)
    { cout <<endl <<endl <<"  Product check of ratios:" <<endl;

      for(nl = 0; nl < nmax; nl++)
	{ cout <<endl <<"  x = " <<xx[nl] <<" :" <<endl;
	  cout <<endl <<cc[nl]*pp[nl] <<endl; } }

}

// Comparison function for naive ordering, See hep-lat/9805026

bool CompSortNaiv(const cl_N & x, const cl_N & y) {

  if((imagpart(x) <= ZERO) && (imagpart(y) <= ZERO)) {
    if(realpart(x) < realpart(y)) {
      return(true);
    }
  }
  if((imagpart(x) <= ZERO) && (plusp(imagpart(y)))) {
    return(true);
  }
  if((plusp(imagpart(x))) && (plusp(imagpart(y)))) {
    if(realpart(x) > realpart(y)) {
      return(true);
    }
  }
  return(false);
}

// This routine returns for idx in [0,2^digits-1] the corresponding
// bit reversal index 

int bitReversalRepresentation(const int idx, const int digits) {
  int res[digits];
  int num = idx;
  for(int i = 0; i < digits; i++) {
    res[i] = 0;
  }
  for(int i = 0; i < digits; i++) {
    int k = digits - i;
    int p = int(pow(2., (digits-(i+1))) );
    res[k] = num/p;
    num = num - res[k]*p;
  }
  num = 0;
  for(int i = 0; i < digits; i++) {
    num = num + res[i]*int(pow(2.,digits-(i+1)));
  }
  return(num);
}

// quicksort template for comparison function (*comp)(const T&, const T&)

template<class T> void quicksort(const int n, T arr[], int idx[], bool (*comp)(const T&, const T&)){
  T v, td;
  int i, j, l, r, ti, tos, stack[32];
  
  l = 0; r = n-1; tos = -1;
  for (;;) {
    while (r > l) {
      v = arr[r]; i = l; j = r-1;
      for (;;){
	while (comp(arr[i], v)) i ++;
	/* j > l prevents underflow */
	while (!comp(arr[j], v) && j > l) j --;
	if (i >= j) break;
	td = arr[i]; arr[i] = arr[j]; arr[j] = td;
	ti = idx[i]; idx[i] = idx[j]; idx[j] = ti;
      }
      td = arr[i]; arr[i] = arr[r]; arr[r] = td;
      ti = idx[i]; idx[i] = idx[r]; idx[r] = ti;
      if (i-l > r-i){
	stack[++tos] = l; stack[++tos] = i-1; l = i+1;
      }
      else{
	stack[++tos] = i+1; stack[++tos] = r; r = i-1;
      }
      if(tos > 31) {
	cerr << "Error in quicksort! Aborting...!" << endl;
	exit(31);
      }
    }
    if (tos == -1) break;
    r = stack[tos--]; l = stack[tos--];
  }
}

// This bring the complex list Roots of length Maxpow to
// naive order.
// See hep-lat/9805026

void NaiveOrder(cl_N * Roots, const int Maxpow) {
  int idx[Maxpow];

  for(int i =0; i < Maxpow; i++) {
    idx[i] = i;
  }
  quicksort(Maxpow, Roots, idx, &CompSortNaiv);
  cout << "Naivly ordered roots" << endl;
  for(int i = 0; i < Maxpow; i++) {
    cout << i << " " << double_approx(realpart(Roots[i])) << " " <<
      double_approx(imagpart(Roots[i])) << endl;
  }
  cout << endl;
}

// This bring the complex list Roots of length Maxpow to
// bit reversal order.
// See hep-lat/9805026

void BitReversalOrder(cl_N *Roots, const int Maxpow, bool Printing=false) {

  int digits = 2;
  int power = 2;
  while(power < Maxpow) {
    power=power*2;
    digits++;
  }

  //  cout << "digits = " << digits << " power " << power << " " << pow(2.,digits-1) << endl;
  
  cl_N paddedRoots[power];
  cl_N reversedRoots[power];

  for(int i = 0; i < Maxpow; i++) {
    paddedRoots[i] = Roots[i];
  }
  for(int i = Maxpow; i < power; i++) {
    paddedRoots[i] = complex(-HUND, ZERO);
  }
  NaiveOrder(paddedRoots, power);
  for(int i = 0; i < power; i++) {
    reversedRoots[i] = paddedRoots[bitReversalRepresentation(i, digits)];
    // cout << i << " " << bitReversalRepresentation(i, digits) << endl;
  }
  
  if(Printing) {
    cout << "Bit reversed ordered roots" << endl;
    for(int i = 0, j=0; i < power; i++) {
      if((realpart(reversedRoots[i]) != -HUND)) {
	Roots[j] = reversedRoots[i];
	cout << j << " " << double_approx(realpart(Roots[j])) << " " <<
	double_approx(imagpart(Roots[j])) <<  endl; 
	j++;
      }
    }
    cout << endl;
  }
  return;
}


/******************************************************************************/
//
//  Calculating the roots of a polynomial approximation minimizing
//  the integral of relative quadratic deviation from x^(-Alpha)
//  in an interval.
//
//  The coefficients of the polynomials are assumed to be known.

//  Input parameters:
//  order of the polynomial:                                Maxpow
//  the (negative) power to be approximated:                Alpha
//  lower bound of interval of approximation:               Epsilon
//  upper bound of interval of approximation:               Lambda
//
//  The name of a array containing the coefficients:        Coef
//
//  The precision:                                          Digit
//
//  name of the file for writing results:                   Filename
//  start file for Start=yes, otherwise append              Start
//  format for assignments, fortran or tao:                 Format
//  print intermediate results for Printing=yes             Printing


void ApproxiRootr(int Maxpow, cl_F Epsilon, cl_F Lambda, cl_F* Coef,
                  float_format_t Digit,
                  char* Filename, char* Format, int Printing)
{
  int  ord, leng;
  
  cl_N Poly[1+Maxpow], Root[Maxpow], Rho[Maxpow];
  cl_F wlst[1+Maxpow];
  cl_N roots[2*Maxpow];
  cl_F pi2 = cl_float(realpart(acos(ZERO)))/HALF/HALF;
  cl_F rr, ang, coef;
  
  
  // Check input
  
  if(Maxpow < 0 || Maxpow > MAXORD) { 
    cout <<endl <<"wrong order in ApproxiRootr:" <<endl;
    cout <<endl <<Maxpow <<endl; 
  }
  
  if(Printing == TRUE)
    cout <<endl <<endl <<" Roots of the polynomial:" <<endl;
  
  // Define polynomial with coefficient array
  
  for(ord = 0; ord < Maxpow+1; ord++) {
    Poly[ord] = complex(Coef[ord],ZERO);
  }
  // Find roots of the polynomial
  
  Polyrootc(Poly,Maxpow,Root,Digit,Printing);

  for(int i = 0; i < Maxpow; i++) {
    Root[i] = HALF*(Root[i]*(Lambda-Epsilon)+Epsilon+Lambda);
  }
  
  // Put them in pairwise optimal order for numerical evaluation
  // OptimordPair(2*Maxpow,roots,Coef[0],Alpha,Epsilon,Lambda,Maxpow/2+1,Printing);  

  // Put them in bit reversal order!

  BitReversalOrder(Root, Maxpow, Printing);

  // Compute monomials
  for(int i = 0; i < Maxpow; i++) {
    roots[i] = sqrt(Root[i]);
    if(!plusp(imagpart(roots[i]))) {
      roots[i] = -roots[i];
    }
    roots[2*Maxpow-i-1] = conjugate(roots[i]);
  }

  cl_N a = EvalPoly(Poly, Maxpow, (TWO*Epsilon-Epsilon-Lambda)/(Lambda-Epsilon));
  cl_N b = EvalPolyProd(Root, Maxpow, Epsilon);
  cl_N norma = a/b;
  cout << endl << a << endl << b << endl;
  if(Printing) {
    cout << endl << "Normierungsfaktor = "<< double_approx(realpart(norma)) << endl
	 << " local = " << double_approx(abs(realpart(exp(log(norma)/cl_R(Maxpow))))) << endl;
  }

  leng = Maxpow;
  WriteRoots("roots.dat", Root, Maxpow);
  WriteRoots("Square_root_BR_roots.dat", roots, 2*Maxpow);
  WriteCoefs(Coef, Maxpow+1);

  // Write the local normalisation factor to a file.
  
  ofstream ofs("normierungLocal.dat", ios::out);
  //  ofs.setf(ios::scientific);
  ofs.precision(20);
  ofs << double_approx(abs(realpart(exp(log(norma)/cl_R(Maxpow))))) << endl;
  ofs.close();

 }
/******************************************************************************/

template <class T> cl_F func(T &x) {
  return(ONE/sqrt(x));
}

// This routine produces coefficients for Chebycheff pol.
// of a given order and interval [epsilon, lambda]

void ChebyCoeff(const int order, 
		const cl_F &epsilon, const cl_F &lambda,
		cl_F * Coeff, cl_F * c) {

  cl_F bma = HALF*(lambda-epsilon);
  cl_F bpa = HALF*(lambda+epsilon);
  cl_F y;
  cl_F ftable[5000];

  for(int i = 0; i < order+1; i++) {
    y = cos(pi(bma)*As(cl_F)((cl_R(i)+HALF)/cl_R(order+1)));
    ftable[i] = func(y*bma+bpa);
  }

  cl_F fac = As(cl_F)(TWO/cl_R(order+1));
  for(int i = 0; i < order+1; i++) {
    cl_F sumit = ZERO;
    for(int j = 0; j < order+1; j++) {
      sumit = sumit + ftable[j]*cos(pi(bma)*(cl_R(i)*(cl_R(j)+HALF)/cl_R(order+1)));
    }
    Coeff[i] = fac*sumit;
  }

  cl_univpoly_real_ring PR = find_univpoly_ring(cl_R_ring);
  cl_UP_R b = PR->create(order);
  for(int i = 0; i < order+1; i++) {
    b.set_coeff(i, ZERO);
  }
  b.set_coeff(0, -HALF*Coeff[0]);
  b.finalize();

  for(int i=0; i < order+1; i++) {
    cl_UP_I C = tschebychev(i);
    for(int j=0; j < i+1; j++) {
      cl_R c = coeff(b, j);
      b.set_coeff(j, c + Coeff[i]*cl_R(coeff(C, j)));
    }
    b.finalize();
  }
  for(int i = 0; i < order+1; i++) {
    c[i] = Coeff[i];
    Coeff[i] = As(cl_F)(coeff(b, i));
    cout << i << " : " << c[i] << endl;
  }
}

// here we use clenshaw to evaluate the polynomial

cl_N EvalCheby(const int order, cl_N * Coeff, cl_N &x, 
	       const cl_N & epsilon, const cl_N lambda) {
  cl_N d=complex(ZERO,ZERO), dd=complex(ZERO,ZERO), sv, z, res;
  int j;

  z = (TWO*x - epsilon - lambda)/(lambda-epsilon);

  for(j=order; j>=1; j--) {
    sv = d;
    d = TWO*z*d - dd + Coeff[j];
    dd = sv;
  }

  res = z*d - dd + HALF*Coeff[0];

  return(res);
}

int main(int argc, char *argv[])
 {
   int  Maxpow = MAXPOW;
   int  Printing = TRUE;

   float_format_t Digit = float_format(DIGIT);

   cl_F Epsilon = EPSILON, Lambda =  LAMBDA;

   double sec(-(double(clock()))/double(CLOCKS_PER_SEC));


// Check order

    if(Maxpow > MAXORD)
     { cout <<"Polynomial order is too large: " 
            <<Maxpow <<" > " <<MAXORD <<endl;
       exit(EXIT_FAILURE); }

// Print out basic parameters

   cout <<"Maxpow  ="  <<Maxpow  <<endl;
   cout <<"Epsilon =" <<Epsilon <<endl;
   cout <<"Lambda  ="  <<Lambda  <<endl;

   ChebyCoeff(Maxpow, Epsilon, Lambda, Coef, c);

//    cl_F x = "0.1e+0_1500";
//    cl_F y = "0.1e+0_1500";
//    for(int i = 0; i < 5; i++) {
//      cout << EvalPoly(Coef, Maxpow, x) << endl;
//      x = x + y;
//    }
//    exit(5);

//    Quadropt(Maxpow,Alpha,Epsilon,Lambda,Digit,Filename,Format,Printing);

//    sec = sec+double(clock())/double(CLOCKS_PER_SEC);
//    cout << endl <<"cpu time for Quadropt: " <<sec <<"s" <<endl;

//    sec = -(double(clock()))/double(CLOCKS_PER_SEC);

//    RestorePolyCoef(Maxpow,Recb,Recg,Orth,Coed,Coef);

   
   ApproxiRootr(Maxpow,Epsilon,Lambda,Coef,Digit,
		Filenamr,Format,Printing);
   
   cout << EvalPoly(Coef, Maxpow, complex((ONE-Epsilon-Lambda)/(Lambda-Epsilon), ZERO)) << endl <<
     EvalCheby(Maxpow, c, HALF, Epsilon, Lambda) << endl;

   sec = sec+double(clock())/double(CLOCKS_PER_SEC);
   cout << endl <<"cpu time for roots: " <<sec <<"s" <<endl;

   exit(EXIT_SUCCESS);
 }

/******************************************************************************/
