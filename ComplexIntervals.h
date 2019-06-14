/*

iRRAM/c_int.h -- header file for class of complex intervals for the iRRAM library
 
Copyright (C) 2015 Norbert Mueller

The underlying formulae are taken from the dissertation of Götz Alefeld
 
This file is part of the iRRAM Library.
 
The iRRAM Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Library General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at your
option) any later version.
 
The iRRAM Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
License for more details.
 
You should have received a copy of the GNU Library General Public License
along with the iRRAM Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. 
*/

/****************************************************************************/
// template class for closed complex intervals 
// 
// C should at least be a ring with unity
// there should at least be a cast of integers into C,
// furthermore x+y, x-y, x*y and norm(x) should exist 
// we do not require completeness of C, so this is slighty more general that a Banach space?
/****************************************************************************/

#ifndef iRRAM_ComplexIntervals_H
#define iRRAM_ComplexIntervals_H

namespace iRRAM {

class c_int
{
public:

// Constructors: -------------------------------

// value: point c_int  [0,0]
c_int():_Ireal(0),_Iimag(0){}

// value: c_int AxB from REAL intervals A,B
c_int(const INTERVAL& A, const INTERVAL& B):_Ireal(A),_Iimag(B){}

// value: c_int Ax{0} from REAL interval A
c_int(const INTERVAL& A):_Ireal(A),_Iimag(INTERVAL(0)){}

// value: c_int from COMPLEX point
c_int(const COMPLEX& z):c_int(INTERVAL(real(z)),INTERVAL(imag(z))){}


// value: c_int from REAL point
c_int(const REAL& z):c_int(INTERVAL(z),INTERVAL(REAL(0))){}

// Standard Arithmetic: ------------------------
//  x°y must be a set containing all points a°b where point a is from x and
// point b is from y

friend c_int  operator +  (const c_int& x, const c_int& y);
friend c_int  operator +  (const COMPLEX& x, const c_int& y){return c_int(x) - y;}
friend c_int  operator +  (const c_int& x, const COMPLEX& y){return x - c_int(y);}
friend c_int  operator +  (const REAL& x, const c_int& y){return c_int(x) + y;}
friend c_int  operator +  (const c_int& x, const REAL& y){return x + c_int(y);}
friend c_int& operator += (c_int& x, const c_int& y){x = x + y; return x;}
friend c_int& operator += (c_int& x, const COMPLEX& y){x= x + c_int(y); return x;}

friend c_int  operator -  (const c_int& x,  const c_int& y);
friend c_int  operator -  (const COMPLEX& x, const c_int& y){return c_int(x) - y;}
friend c_int  operator -  (const c_int& x, const COMPLEX& y){return x - c_int(y);}
friend c_int  operator -  (const REAL& x, const c_int& y){return c_int(x) - y;}
friend c_int  operator -  (const c_int& x, const REAL& y){return x - c_int(y);}
friend c_int& operator -= (c_int& x, const c_int& y){x = x - y; return x;}

friend c_int  operator -  (const c_int& x);

friend c_int  operator *  (const c_int& x, const c_int& y);
friend c_int  operator *  (const COMPLEX& x, const c_int& y){return c_int(x) * y;}
friend c_int  operator *  (const c_int& x, const COMPLEX& y){return x * c_int(y);}
friend c_int  operator *  (const REAL& x, const c_int& y){return c_int(x) * y;}
friend c_int  operator *  (const c_int& x, const REAL& y){return x * c_int(y);}
friend c_int& operator *= (c_int& x, const c_int& y){x = x * y; return x;}
friend c_int& operator *= (c_int& x, const COMPLEX& y){x = x * c_int(y); return x;}

friend c_int  operator /  (const c_int& x, const c_int& y);
friend c_int  operator /  (const COMPLEX& x, const c_int& y){return c_int(x) / y;}
friend c_int  operator /  (const c_int& x, const COMPLEX& y){return x / c_int(y);}
friend c_int  operator /  (const REAL& x, const c_int& y){return c_int(x) / y;}
friend c_int  operator /  (const c_int& x, const REAL& y){return x / c_int(y);}
friend c_int& operator /= (c_int& x, const c_int& y){x = x / y; return x;}
friend c_int& operator /= (c_int& x, const COMPLEX& y){x = x / c_int(y); return x;}

//friend REAL wid(const c_int& x);
//friend REAL inf(const c_int& x);
//friend REAL sup(const c_int& x);
//friend REAL mid(const c_int& x);
friend REAL mag(const c_int& x);
//friend REAL mig(const c_int& x);

// friend LAZY_BOOLEAN superset (const c_int& x,
//                              const c_int& y);
// friend LAZY_BOOLEAN proper_superset (const c_int& x,
//                              const c_int& y);
// friend LAZY_BOOLEAN subset (const c_int& x,
//                              const c_int& y);
// friend LAZY_BOOLEAN proper_subset (const c_int& x,
//                              const c_int& y);
// friend LAZY_BOOLEAN in_interior (const c_int& x,
//                              const c_int& y);
// friend LAZY_BOOLEAN disjoint (const c_int& x,
//                              const c_int& y);
// friend LAZY_BOOLEAN in (const C& x,
//                              const c_int& y);
// 
// friend c_int interval_hull (const c_int& x,
//                              const c_int& y);
friend c_int intersect (const c_int& x, const c_int& y);
friend c_int join (const c_int& x, const c_int& y);
friend std::vector<c_int> split (const c_int& x);
 
// friend c_int fabs(const c_int& x);

friend c_int power(const c_int& x,int n);
friend c_int square(const c_int& x);

//friend c_int exp(const c_int& x);
//friend c_int log(const c_int& x);
//friend c_int sin(const c_int& x);
//friend c_int cos(const c_int& x);
//friend c_int tan(const c_int& x);
//friend c_int asin(const c_int& x);
//friend c_int acos(const c_int& x);
//friend c_int atan(const c_int& x);

//private :

INTERVAL _Ireal,_Iimag;
};




c_int operator + (const c_int & A, const c_int & B){
	return c_int( A._Ireal + B._Ireal, A._Iimag + B._Iimag);
}


c_int operator - (const c_int & A, const c_int & B){
	return c_int( A._Ireal - B._Ireal, A._Iimag - B._Iimag);
}

c_int operator - (const c_int & B){
	return c_int( - B._Ireal, - B._Iimag);
}


c_int operator * (const c_int & A, const c_int & B){
	return c_int( A._Ireal * B._Ireal - A._Iimag * B._Iimag,
		      A._Ireal * B._Iimag + A._Iimag * B._Ireal);
}

c_int operator / (const c_int & A, const c_int & B){
        INTERVAL denom(square(B._Ireal)+ square(B._Iimag));
	return c_int( (A._Ireal * B._Ireal + A._Iimag * B._Iimag)/denom,
		      (A._Iimag * B._Ireal - A._Ireal * B._Iimag)/denom);
}

c_int square2(const c_int & x){
  return x*x;
}

c_int square3(const c_int & x){
  std::vector<c_int> v=split(x);
  return join(join(v[0]*v[0],v[1]*v[1]),join(v[2]*v[2],v[3]*v[3]));
}

c_int square(const c_int & x){
  return c_int(square(x._Ireal)-square(x._Iimag),REAL(2)*(x._Ireal*x._Iimag));
}

c_int inv2(const c_int & x){
  std::vector<c_int> v=split(x);
  c_int one(COMPLEX(1));
  return join(join(one/v[0],one/v[1]),join(one/v[2],one/v[3]));
}




c_int power(const c_int & x,int n){
  if (n==0)  return c_int(REAL(1));
  if (n==1)  return x;
  if (n==2)  return square(x);
  if (n==3)  return square(x)*x;
  if (n<0)   return power(COMPLEX(1)/x,-n);
  if (n%2==1) return square(power(x,n/2))*x;
  return square(power(x,n/2));
}

REAL mag(const c_int & x){
  return maximum(
	 maximum(abs(COMPLEX(x._Ireal.low,x._Iimag.low)),
		 abs(COMPLEX(x._Ireal.low,x._Iimag.upp))),
	 maximum(abs(COMPLEX(x._Ireal.upp,x._Iimag.low)),
		 abs(COMPLEX(x._Ireal.upp,x._Iimag.upp))));
}

c_int intersect (const c_int& x, const c_int& y){
   return c_int(INTERVAL(maximum(x._Ireal.low,y._Ireal.low),minimum(x._Ireal.upp,y._Ireal.upp)),
	INTERVAL(maximum(x._Iimag.low,y._Iimag.low),minimum(x._Iimag.upp,y._Iimag.upp)));}

c_int join (const c_int& x, const c_int& y){
   return c_int(INTERVAL(minimum(x._Ireal.low,y._Ireal.low),maximum(x._Ireal.upp,y._Ireal.upp)),
	INTERVAL(minimum(x._Iimag.low,y._Iimag.low),maximum(x._Iimag.upp,y._Iimag.upp)));}

std::vector<c_int> split (const c_int& x){
   std::vector<c_int> v(4);
   v[0]=c_int(INTERVAL(x._Ireal.low,(x._Ireal.low+x._Ireal.upp)/2),INTERVAL(x._Iimag.low,(x._Iimag.low+x._Iimag.upp)/2));
   v[1]=c_int(INTERVAL(x._Ireal.low,(x._Ireal.low+x._Ireal.upp)/2),INTERVAL((x._Iimag.low+x._Iimag.upp)/2,x._Iimag.upp));
   v[2]=c_int(INTERVAL((x._Ireal.low+x._Ireal.upp)/2,x._Ireal.upp),INTERVAL((x._Iimag.low+x._Iimag.upp)/2,x._Iimag.upp));
   v[3]=c_int(INTERVAL((x._Ireal.low+x._Ireal.upp)/2,x._Ireal.upp),INTERVAL(x._Iimag.low,(x._Iimag.low+x._Iimag.upp)/2));
   return v;
}

} /* ! namespace iRRAM */

#endif /* !  iRRAM_ComplexIntervals_H */
