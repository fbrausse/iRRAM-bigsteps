/*

TaylorModel.h - Prototype for an implementation of Taylor Models in iRRAM

Copyright (C) 2015 Norbert Mueller, Franz Brausse
 
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

/* iRRAM version of simplified Taylor Models.
 * Each error symbols is implemented as a REAL variable, (i.e. an interval!)
 * 
 * TODO: 
 * (a) This version is restricted to only a subset of interesting functions
 * and to order 1. 
 * (b) Currently, there is no division of models, and no versions for
 * trigonometric functions.
 * (c) We should assign error symbols to basic constants like pi, ln(2) and e
 */


#ifndef iRRAM_TaylorModel_H
#define iRRAM_TaylorModel_H


#include <iRRAM/lib.h>
#include <vector>


namespace iRRAM {




class TM {
public: 
	struct I {
		unsigned long long id;
		REAL ci;
	};

private:
/* The max_id has to be unique to distinguish different error symbols.
 * TODO: For multithreaded versions of the iRRAM, this should be atomic!
 */
	static unsigned long long max_id;
	static unsigned default_sweep;
	static int prec_diff;
	static unsigned sweepto;  /* For testing purposes: where should sweeping be directed to? */

/* helper function to destruct the internal representation of REAL x into
 * center and radius, where both are REAL themselves.
 * 
 * TODO: write a version that does not have to transform a double interval into
 * a multiple precision interval first.
 * 
 * TODO: There is a transformation of the error mantissa into a signed version.
 * Maybe we should a conversion of unsigned to REAL
 */
	static void my_to_formal_ball(const REAL&x, REAL&center, REAL&radius) {
		DYADIC d;
		sizetype r;
		x.to_formal_ball(d,r);
		center=REAL(d);
		radius=scale(REAL(signed(r.mantissa)),r.exponent);
	}


public:              /* TODO: This should be private in the new future! */
	REAL c0;           /* center of the Taylor Model */
	std::vector<I> c;  /* error symbols with their IDs */

//   TM(TM&& t): c(std::move(t.c)),sweepto(t.sweepto){c0=t.c0;}
//   TM(const TM& t): c(t.c),sweepto(t.sweepto){c0=t.c0;}
public:
	/* default constructor */
	TM(){}

	/* Constructor for conversion from REAL to TM
	 * This does NOT introduce an error symbol!
	 * The purpose is just to allow an automatic conversion and 
	 * to simplify the formulation of algorithms.
	 */
	TM(const REAL& x) : c0(x) {}

	// conversion from a TM to a REAL
	REAL to_real() const {
		REAL zero_one=0;
		sizetype zo;
		sizetype_set(zo,1,0);
		zero_one.seterror(zo);
		REAL s=0;
		for (const I &i : c) {
			s += i.ci * zero_one;
		}
		return s + c0;
	}

	// explicit cast from a TM to a REAL
	explicit inline operator REAL(void) const {
		return to_real();
	}

	static void set_default_sweep(unsigned swp) { default_sweep=swp; }
	static void set_prec_diff    (unsigned prd) { prec_diff=prd; }


	void geterror_info(sizetype& error) {
		c0.geterror(error);
		sizetype error2;
		for (I &i : c) {
			i.ci.getsize(error2);
			sizetype_inc(error,error2);
		}
	}


	/* The following routines manipulate the Taylor Model in different ways 
	 * trying to find some useful "rounding" procedures
	 * None of them is perfect...
	 */

	//merge the rounding error and all error symbols into a new error symbol
	void round();

	//merge the rounding error and all but the largest error symbol and into a new error symbol
	void round0();

	//merge the rounding error and the smallest error symbol into a new error symbol
	void round1();

	//merge the rounding error and the two smallest error symbols into a new error symbol
	void round2();

	//merge the rounding error and all but the largest error symbol and into a new error symbol
	// in case the own error in the largest symbol is quite large, merge this symbol as well.
	void round5();

	//merge the rounding error and all but the largest error symbol and into a new error symbol
	// in case the own error in the largest symbol is quite large, merge this symbol as well and 
	// try to keep another large symbol.
	void round6();

	//merge the rounding error and  all but the two largest error symbols into a new error symbol
	void round7();


	//merge the rounding error and the smallest error symbol into a new error symbol
	//but only if the rounding error is larger than the largest symbol
	void check();


	//test whether the rounding error is larger than the largest symbol
	/* private */
	bool test2();

	//merge the rounding error and all error symbols into a new error symbol
	//but only if the rounding error is larger than the largest symbol
	static void check2(TM& x) {
		if (x.test2())
			x.round();
	}

	//vector version of check2()
	static void check2(std::vector<TM>& x) {
		bool test = false;
		for (unsigned i=0; i<x.size(); i++) {
			test = test || x[i].test2();
		}
		if (test) { 
	//    cerr << "# check\n";
			for (unsigned i=0; i<x.size(); i++) {
				x[i].round();
			}
		}
	}


	//test whether the rounding error is larger than the accumulated error symbols
	// or whether an error symbol is close to exhaustion (i.e. it's own error is too large) 
	/* private */
	bool test5();


	//merge the rounding error and all but the largest error symbols into a new error symbol
	//but only if the rounding error is larger than the accumulated error symobls
	// also merge the largest symbol if its own error is too large
	static void check5(TM & x){
		if (x.test5())
			x.round5();
	}

	//vector version of check5()
	static void check5(std::vector<TM>& x){
		bool test=false;
		for (unsigned i=0; i< x.size(); i++){
			test=test || x[i].test5();
		}
		if ( test ){ 
			for (unsigned i=0; i< x.size(); i++){
				x[i].round5();
			}
		}
	}


	//test whether the rounding error is larger than the accumulated error symbols
	// or whether an error symbol is close to exhaustion (i.e. it's own error is too large) 
	/* private */
	bool test6();

	static void check6(TM & x){
		if ( x.test6() ){ 
			x.round6();
		}
	}
	//vector version of check6()
	static void check6(std::vector<TM>& x){
		bool test=false;
		for (unsigned i=0; i< x.size(); i++){
			test=test || x[i].test6();
		}
		if ( test ){ 
			for (unsigned i=0; i< x.size(); i++){
				x[i].round6();
			}
		}
	}


	//test whether the rounding error is larger than the error of one of the error symbols
	// or whether an error symbol is close to exhaustion (i.e. it's own error is too large) 
	/* private */
	bool test7();

	static void check7(TM & x){
		if ( x.test7() ){ 
			x.round7();
		}
	}
	//vector version of check6()
	static void check7(std::vector<TM>& x){
		bool test=false;
		for (unsigned i=0; i< x.size(); i++){
			test=test || x[i].test7();
			if (test) break;
		}
		if ( test ){ 
			for (unsigned i=0; i< x.size(); i++){
				x[i].round7();
			}
		}
	}





#ifndef TM_CHECK
	inline static void polish(TM& x){ check7(x); }
	inline static void polish(std::vector<TM>& x){ check7(x); }
#else
	inline static void polish(TM& x){ TM_CHECK(x); }
	inline static void polish(std::vector<TM>& x){ TM_CHECK(x); }
#endif



	//test whether the rounding error is larger than the smallest symbol
	/* private */
	bool test3();

	//merge the rounding error and all error symbols into a new error symbol
	//but only if the rounding error is larger than the smallest symbol
	static void check3(TM& x) {
		if (x.test3())
			x.round();
	}

	//vector version of check3()
	static void check3(std::vector<TM>& x) {
		bool test=false;
		for (unsigned i=0; i < x.size(); i++) {
			test = test || x[i].test3();
		}
		if (test) { 
			cerr << "# check\n";
			for (unsigned i=0; i < x.size(); i++) {
				x[i].round();
			}
		}
	}

	//test whether the rounding error is larger than the accumulated error symbols
	/* private */
	bool test4();

	//merge the rounding error and all error symbols into a new error symbol
	//but only if the rounding error is larger than the accumulated error symbols
	static void check4(TM& x) {
		if (x.test4())
			x.round();
	}


	TM & operator+=(const TM &tm)
	{
		{ stiff code(prec_diff); c0+=tm.c0; }
		for (const I &i : tm.c) {
			for (I &j : c) {
				if (i.id == j.id) {
					j.ci += i.ci;
					goto ok;
				}
			}
			c.push_back(i);
ok:;
		}
		return *this;
	}

	TM & operator-=(const TM &tm)
	{
		{ stiff code(prec_diff); c0 -= tm.c0; }
		for (const I &i : tm.c) {
			for (I &j : c) {
				if (i.id == j.id) {
					j.ci -= i.ci;
					goto ok;
				}
			}
			c.push_back({i.id,-i.ci});
ok:;
		}
		return *this;
	}

	TM & operator*=(const REAL &r) {
		{ stiff code(prec_diff); c0 *=r; }
		for (I &l : c)
				l.ci *= r;
		return *this;
	}


	friend TM operator*(const TM &q, const TM &r) {
		TM f = q*r.c0;
		REAL q_real = q.to_real();
		if (q.sweepto == 0) {
			for (const I &i : r.c)
				f.c0 += q_real*i.ci;
			return f;
		} else {  
			TM g = TM(REAL(0));
			for (const I &i : r.c)
				g.c.push_back(I({i.id,q_real*i.ci}));
			return f+g;
		}
	}


	friend TM operator+(const TM &q, const TM &r)   { return TM(q) += r; }
	friend TM operator-(const TM &q, const TM &r)   { return TM(q) -= r; }
	friend TM operator*(const TM &q, const REAL &r) { return TM(q) *= r; }


	friend orstream & operator<<(orstream &o, const TM &p)
	{
		o << swrite(p.c0,20,iRRAM_float_show)
		  << " swpt: " << p.sweepto
		  << " i.e. "
		  << p.c0.vsize.mantissa << "*2^" << p.c0.vsize.exponent
		  << " +- "
		  << p.c0.error.mantissa << "*2^" << p.c0.error.exponent;
		for (const I &i : p.c) {
			o << "\n   + (" << swrite(i.ci,25,iRRAM_float_show) << ")"
			  << "*[" << i.id << "]"
			  <<" i.e. "
			  << i.ci.vsize.mantissa << "*2^" << i.ci.vsize.exponent
			  << " +- "
			  << i.ci.error.mantissa << "*2^" << i.ci.error.exponent;
		}
		return o;
	}
}; // !class TM 

} // !namespace iRRAM

#endif /* !  iRRAM_TaylorModel_H */
