/*

PIVP.cc - Solver for Polynomial Initial Value Problems

Copyright (C) 2010-2015 Norbert Mueller, Franz Brausse
 
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



#include <iRRAM.h>
#include <vector>
#include <list>
#include <cassert>
#include <iRRAM/limit_templates.h>

#include <sys/time.h>
#include <cstdarg>

#include "Polynomials.h"
#include "TaylorModel.h"
#include "TaylorSeries.h"
#include "VectorExtensions.h"

#define dbg(fmt, ...) iRRAM_DEBUG0(2,{fprintf(stderr, (fmt), ##__VA_ARGS__);})


using namespace iRRAM;


class POLYNOMIAL_FLOW {
public:
	typedef VI I_type;

	class Iterator {
		const POLYNOMIAL_FLOW &F;
		const unsigned nu;
		unsigned state;
	public:
		Iterator(const POLYNOMIAL_FLOW &F, unsigned nu)
		: F(F), nu(nu), state(0) {}
		void operator++() { state++; }
		explicit operator bool() const { return state < F.c[nu].size(); }
		unsigned size() const { return F.dimension() + 1; }
		const I_type * operator->() const { return &**this; }
		const I_type & operator*() const { return F.c[nu][state]; }
		int operator[](unsigned j) const { return (*this)->ik[j]; }
	};

	typedef Iterator iterator_type;

	POLYNOMIAL_FLOW(unsigned dimension)
	: _d(dimension), _mu(0), _autonomous(true), c(dimension)
	{	for (unsigned nu=0; nu<dimension; nu++)
			poly.push_back(POLYNOMIAL2(dimension+1)); //no content yet, just to get right size!
	}

	POLYNOMIAL_FLOW(const POLYNOMIAL_FLOW &old, const REAL &told)
	: _d(old._d), _mu(0), _autonomous(old.is_autonomous()), c(old._d)
	{
		unsigned d = _d;
		for (unsigned nu=0; nu<d; nu++)
			poly.push_back(POLYNOMIAL2(d));
		for (unsigned nu=0; nu<_d; nu++) {
			POLYNOMIAL2 p(d+1);
			for (Iterator it = old.iterator(nu, old.mu()); it; ++it) {
				unsigned k = it[d];
				REAL c = old(nu,it);
				std::vector<unsigned> i = it->ik;
				POLYNOMIAL2 q(d+1, k+1);
				INTEGER bin_coeff = 1;
				REAL told_power = 1;
				for (unsigned l=0; l<=k; l++) {
					i[d] = k-l;
					// q.c.push_back({ i, times_bin_coeff(c, k, l) * power(told, k-l) });
					q.add_coeff({ i, c * REAL(bin_coeff) * told_power });
					bin_coeff *= k - l;
					bin_coeff /= l + 1;
					told_power *= told;
				}
				p += q;
			}
			poly[nu]=p;
			for (const POLYNOMIAL2::I &ik : p.coefficients())
				c[nu].push_back(ik);
			_mu = std::max(_mu, p.degree());
		}
	}

	void add_coeff(unsigned nu, const I_type &coeff)
	{
		assert(nu < _d);
		assert(coeff.ik.size() == _d + 1);
		c[nu].push_back(coeff);
		poly[nu].add_coeff(coeff);
		for (unsigned j=0; j<=_d; j++)
			_mu = std::max(_mu, coeff.ik[j]);
		if (coeff.ik[_d] != 0)
			_autonomous = false;
	}

	inline unsigned dimension()     const { return _d; }
	inline int      mu()            const { return _mu; }
	inline bool     is_autonomous() const { return _autonomous; }

	iterator_type iterator(unsigned nu, unsigned l) const
	{
		return Iterator(*this, nu);
	}

	REAL operator()(unsigned nu, const iterator_type &idx) const
	{
		return idx->c;
	}

	REAL UF_bruteforce(
		std::vector<REAL> w, REAL t0, const REAL &delta, const REAL &eps
	) const {
		for (REAL &wj : w)
			wj = abs(wj) + eps;
		t0 = abs(t0) + delta;

		REAL max_result(0);
		for (unsigned nu=0; nu<dimension(); nu++) {
			REAL result(0);
			for (Iterator it = iterator(nu, _mu); it; ++it) {
				REAL c = (*this)(nu, it);
				REAL m(1);
				if (it->ik[_d] > 0)
					m *= power(t0, it->ik[_d]);
				for (unsigned j : it->idx_i_gt0)
					m *= power(w[j], it->ik[j]);
				result += m * abs(c);
			}
			max_result = nu == 0 ? result : maximum(max_result, result);
		}
		return max_result;
	}

// determine upperbound of the right hand side for complex parameters w (with eps) and t (with delta)
	REAL UF_interval(
		std::vector<REAL> w, const REAL &eps, REAL t, const REAL &delta
	) const {
		std::vector<c_int> w_int(w.size()+1);
		for (unsigned i=0; i<w.size(); i++)
			w_int[i] = c_int(INTERVAL(w[i]-eps,w[i]+eps,true),
			                 INTERVAL(    -eps,     eps,true));
		w_int[w.size()] = c_int(INTERVAL(t-delta,t+delta,true),
		                        INTERVAL( -delta,  delta,true));
		REAL max_result(0);
		for (unsigned nu=0; nu<dimension(); nu++) {
			c_int res = poly[nu](w_int);
			REAL  mx  = mag(res);
			max_result = maximum(mx, max_result);
		}
		return max_result;
	}

	void get_RM(
		const std::vector<REAL> &w,
		const REAL &t0,
		const REAL &delta,
		const REAL &eps,
		REAL &R,
		REAL &M
	) const {
		REAL _UF = UF_bruteforce(w,t0,delta,eps);
		// cout << "U: " << _UF << "\n";
		R = minimum(delta, eps/_UF);
		REAL abs_w(0);
		for (unsigned j = 0; j < w.size(); j++) {
			abs_w += square(w[j]);
		}
		abs_w = sqrt(abs_w);
		//M = abs_w + R*_UF;
		// besser: _UF neu berechnen, aber das macht nur bei nicht-autonomen
		// Systemen einen Unterschied
		// NM, 2013-08-31
		M = abs_w + R*UF_bruteforce(w,t0,R,eps);
	}

	REAL f_abs(
		unsigned nu, const std::vector<REAL> &w_abs, const REAL &eps,
		const REAL &t_abs, const REAL &delta
	) const {
		REAL r(0);
		for (Iterator it = iterator(nu, _mu); it; ++it) {
			REAL m = abs((*this)(nu, it));
			if (it[dimension()] != 0)
				m *= power(t_abs+delta, it[dimension()]);
			for (unsigned j=0; j<it->ni_gt0; j++)
				m *= power(w_abs[it->idx_i_gt0[j]] + eps,
				           it[it->idx_i_gt0[j]]);
			r += m;
		}
		return r;
	}

	REAL f_max(const std::vector<REAL> &w_abs, const REAL &eps, const REAL &t_abs, const REAL &delta=0) const
	{
		REAL m = f_abs(0, w_abs, eps, t_abs, delta);
		for (unsigned nu=1; nu<dimension(); nu++)
			m = maximum(m, f_abs(nu, w_abs, eps, t_abs, delta));
		return m;
	}

	void get_RM2(
		std::vector<REAL> w,
		const REAL &t0,
		const REAL &delta,
		REAL &eps,
		REAL &R,
		REAL &M, 
	        const int step_control_alg
	) const {
//		stiff area(-ACTUAL_STACK.prec_step/2);
		REAL abs_w(0);
		for (unsigned j = 0; j < w.size(); j++) {
			abs_w += square(w[j]);
			w[j] = abs(w[j]);
		}
		REAL t = abs(t0) + delta;
		abs_w = sqrt(abs_w);

		REAL lower_sum = 0;
		REAL rect_w = eps;
		REAL R_simple;
		REAL eps_simple=eps;
		REAL eps_opt;
		int kl=0,kr=0,ks=0;

		unsigned k =100;
		for (unsigned j=0; j<k; j++) {
			rect_w = 9*rect_w/8;
			REAL tmp=rect_w / f_max(w, rect_w, t);
			if (tmp > R_simple || tmp > 999*R_simple/1000) {R_simple=tmp; eps_simple=rect_w;ks=-j;}
//			cout << "# stepcontrol " <<j<< " # "<<rect_w<< " : "<< tmp/9<< " sum: " << lower_sum << " vs. " <<  tmp<<"\n";;
			lower_sum += tmp/9;
			if ( tmp < lower_sum/30 || tmp< lower_sum/40) { kl=-j; break;}
		}
		eps_opt=rect_w;
		rect_w = eps;
		for (unsigned j=0; j<k; j++) {
			REAL tmp=rect_w / f_max(w, rect_w, t);
			if (tmp > R_simple || tmp > 999*R_simple/1000) {R_simple=tmp; eps_simple=rect_w;ks=j;}
//			cout << "# stepcontrol " <<j<< " # " <<rect_w<< " : "<< tmp/16<< " sum: " << lower_sum << " vs. " <<  tmp<<"\n";;
			lower_sum += tmp/16;
			rect_w =15*rect_w/16;
			if ( tmp < lower_sum/60 || tmp< lower_sum/70) { kr=j; break;}
		}
		lower_sum += rect_w / f_max(w,rect_w, t);

		REAL R_opt = minimum(delta, lower_sum);
		REAL M_opt = abs_w + eps_opt;
		eps=eps_simple;
		REAL M_simple=abs_w + eps_simple;
		cout << "# Radius  " << R_opt << " vs. " <<  R_simple<<"\n";
		cout << "# Maximum " << M_opt << " vs. " <<  M_simple<<"\n";
		cout << "# Epsilon " << eps_opt << " vs. " <<  eps_simple<<"\n";
		cout << "# kl/ks/kr "<< kl << " " << ks << " " <<  kr<<"\n";
		eps=eps_simple;

		if (step_control_alg) {
			R=R_simple; M=M_simple;
		}  else {
			R=R_opt; M=M_opt;
		}
	}

	void get_RM3(
		const std::vector<TM>& w_tm,
		const REAL &t,
		const REAL &delta,
		const REAL &R_scale,
		REAL &eps,
		REAL &R,
		REAL &M, 
		REAL &Lo,
		REAL &Rs, 
		REAL &Ro,
		const FUNCTION<std::vector<TM>,unsigned int> a
	) const {
		std::vector<REAL> w_copy; w_copy.reserve(w_tm.size());
		std::vector<REAL> w_abs; w_abs.reserve(w_tm.size());
		REAL abs_w(0);
		for (unsigned j = 0; j < w_tm.size(); j++) {
			w_copy.push_back(w_tm[j].to_real());
			w_abs.push_back(abs(w_copy[j]));
///**/			abs_w=maximum(abs_w, abs(w_copy[j]));
/**/			abs_w += square(w_copy[j]);
		}

		REAL t_abs = abs(t);
/**/		abs_w = sqrt(abs_w);


		REAL lower_sum = 0;
		REAL R_simple;
		REAL eps_opt;
		REAL eps_simple=eps;
		REAL eps_current = eps;
		int test_size=100;
		int simple=0;
		int low=0,high=0; // the sampling fo the integral will be saved in [low..high] (including)

		REAL R_delta[2*test_size]; 
		REAL epstest[2*test_size]; // R_delta[j] contains \int_{epstest[j-1]}^{epstest[j]} 1/U(s) ds

		REAL Rtest[2*test_size];
		int factor=20;

		for (int j=test_size;; j++) { // increasing epsilon
			eps_current = factor*eps_current/(factor-1);
			REAL _UF2 = UF_interval(w_copy,eps_current,t,delta);
			REAL _UF3 = f_max(w_abs, eps_current, t_abs, delta);
			REAL tmp = eps_current /minimum(_UF2,_UF3);
			if (tmp > R_simple || tmp > 999*R_simple/1000) {
				R_simple=tmp;
				eps_simple=eps_current;
				simple=j;
			}
			// R_simple: here we want to find the best combinaton radius via ordinary Picard-Lindeloef
			// at the same time, we determine the corresponding height of the rectangle that leads to R_simple
			lower_sum += tmp/factor;// lower_sum is just for termination!
			R_delta[j]=tmp/factor;
			epstest[j]=eps_current;

			if ( bool(tmp < lower_sum/25 || tmp< lower_sum/30) || j>=2*test_size-1) { high=j; break;} // lower_sum is just for termination!
		}
//		R_delta epstest have been initialized for indices [test_size..high]

		eps_opt=eps_current; // the maximal epsilon taken into account = the height of the final rectangle

		eps_current = eps; // for decreasing epsilon,we again start at eps
		for (int j=test_size-1;; j--) { // decreasing epsilon
			REAL _UF2 = UF_interval(w_copy,eps_current,t,delta);
			REAL _UF3 = f_max(w_abs, eps_current, t_abs, delta);
			REAL tmp = eps_current / minimum(_UF2,_UF3);
			if (tmp > R_simple || tmp > 999*R_simple/1000) {
				R_simple=tmp;
				eps_simple = eps_current;
				simple = j;
			}
			// R_simple: here we want to find the best combinaton radius via ordinary Picard-Lindeloef
			// at the same time, we determine the corresponding height of the rectangle that leads to R_simple
			lower_sum += tmp/factor;// lower_sum is just for termination!
			R_delta[j] = tmp/factor;
			epstest[j] = eps_current;

			eps_current = (factor-1)*eps_current/factor;
			if (bool(tmp < lower_sum/90 || tmp< lower_sum/100) || j==1) { // lower_sum is just for termination!
				low=j-1;
				break;
			}
		}
//		R_delta epstest have been initialized for indices [low+1..test_size-1]and [test_size..high]


		REAL _UF2 = UF_interval(w_copy,eps_current,t,delta);
		REAL _UF3 = f_max(w_abs, eps_current, t_abs, delta);
		REAL tmp=eps_current / minimum(_UF2,_UF3);

		R_delta[low]=tmp; //initial part, \int_{}^{epstest[j]} 1/U(s) ds
		epstest[low]=eps_current;
//		R_delta epstest have been initialized for indices [low..high]

		REAL Rsum= 0;
		for (int j= low;j<=high;j++){
			Rsum+=R_delta[j];Rtest[j]=Rsum; 
		}

		REAL R_opt = minimum(delta, Rsum);
		REAL M_opt = abs_w + eps_opt;

		// Now compute the Lipschitz bound L_opt for radius R_opt
		REAL L_opt=(abs_w+epstest[high])/epstest[low];
		{
			int k=high,j=k-1;
			while (j>low) {
				if (Rtest[k]-Rtest[j] < R_scale*R_opt)
					j--;
				else {
					L_opt= minimum(L_opt,(abs_w+epstest[k])/epstest[j]);
					k--;
				}
			}
		}
		// R_simple, M_simple and eps_simple correspond to the values for ordinary Picard-Lindeloef
		REAL M_simple=abs_w + eps_simple;
//debug		cout << "# Radius  " << R_opt << " vs. " <<  R_simple<< " at w0: "<< w_copy[0]<< " w1: "<< w_copy[1]<<"\n";
//debug		cout << "# Maximum " << M_opt << " vs. " <<  M_simple<< " at w0: "<< w_copy[0]<< " w1: "<< w_copy[1]<<"\n";
//debug		cout << "# Epsilon " << eps_opt << " vs. " <<  eps_simple<< " at  w0: "<< w_copy[0]<< " w1: "<< w_copy[1]<<"\n";
//debug		cout << "# l/s/r "<< low-test_size << " " << simple-test_size << " " <<  high-test_size<< " at  w0: "<< w_copy[0]<< " w1: "<< w_copy[1]<<"\n";

		eps=eps_simple; // return eps as the initial start value for the next step


		// Now we try to improve the radius of convergence using information about the 
		// real trajectory (that we can approximate reliably using and R_opt,M_opt 
		REAL R0 = R_opt,R2=0;
		REAL y_delta, R_test;
		for (int k=100; k<130; k++){
			R_test=k*R0/150;
			REAL yval=poly_bound(a,30,R_test,0,false);
			REAL yerr=M_opt*R0/(R0-R_test)*power(R_test/R0,30+1);
			y_delta= yval+yerr;
			int i=low;
			while ( (i<high) && bool(y_delta>epstest[i] && y_delta>epstest[i]*0.999) ) i++;
			R2=  minimum(delta,R_test+Rtest[high]-Rtest[i]);
//			cerr << "(R',M',R,M)= "<<R_test<<" "<< y_delta<<" "<<R2<<" "<<M_opt<<"\n";
			R0=maximum(R2,R0);
		}


		R=R0;
		M=M_opt;

		Rs = R_simple;
		Ro = R_opt;
		Lo = L_opt;
	}

private:
	unsigned _d;
	unsigned _mu;
	bool _autonomous;
	std::vector<std::vector<I_type>> c;
public:
	std::vector<POLYNOMIAL2> poly;
};

static REAL parse_REAL(const char *s)
{
	if (!strcmp(s, "pi"))
		return pi();
	else if (!strcmp(s, "e"))
		return euler();
	else if (strchr(s, '/'))
		return RATIONAL(s);
	else
		return REAL(s);
}

POLYNOMIAL_FLOW read_poly_flow(const char *fname)
{
	REAL r;
	unsigned dimension, nu;
	const char *msg = NULL;

	irstream in(fname);

	in >> dimension;

	POLYNOMIAL_FLOW F(dimension);

	std::vector<unsigned> ik(dimension+1);
	std::string s;
	unsigned line = 1, column = 0;

	if (!in)
		{ msg = "d"; goto err; }

	while (++line, (in >> nu)) {				/* nu */
		if (!(in >> ik[dimension]))			/* k */
			{ column = 2; msg = "k"; goto err; }
		for (unsigned j=0; j<dimension; j++)
			if (!(in >> ik[j]))			/* i1,...,id */
				{ column = 3+j; msg = "i_j"; goto err; }
		if (!(in >> s) || !s.length())			/* c_{nu,k,i1,...,id} */
			{ column = 3+dimension; msg = "real number c_{nu,k,i_1,...,i_d}"; goto err; }

		F.add_coeff(nu, VI(parse_REAL(s.c_str()), ik));
	}

	return F;
err:
	fprintf(stderr, "invalid coefficient input file, line %u, value %u: expected <%s>\n", line, column, msg);
	exit(1);
}


struct Timer {
	unsigned long t = 0;

	void start() { struct timeval tv; gettimeofday(&tv, NULL); t -= tv.tv_sec * 1000000UL + tv.tv_usec; }
	void stop() { struct timeval tv; gettimeofday(&tv, NULL); t += tv.tv_sec * 1000000UL + tv.tv_usec; }
} t_auto, t_rest0, t_rest1;



struct F_REAL {
	TM x;
	bool valid;

	inline F_REAL() : valid(false) {}
	inline explicit F_REAL(const TM &v) : x(v), valid(true) {}

/***************************************************************
 *  Usage of the below stuff indeed slows down the computation *
 ***************************************************************

	inline explicit F_REAL(const REAL &v, bool is_zero = false)
	: x(v), valid(true), is_zero(is_zero) {}

	inline F_REAL operator+(const F_REAL &b) const
	{
		return is_zero ? b : b.is_zero ? *this : F_REAL(x + b.x);
	}

	inline F_REAL & operator+=(const F_REAL &b)
	{
		if (b.is_zero)
			return *this;
		if (is_zero)
			*this = b;
		else
			x += b.x;
		return *this;
	}

	inline F_REAL operator*(const F_REAL &b) const
	{
		return is_zero || b.is_zero ? F_REAL(0) : F_REAL(x * b.x);
	}

	inline F_REAL & operator*=(const F_REAL &b)
	{
		if (is_zero)
			return *this;
		if (b.is_zero)
			*this = b;
		else
			x *= b.x;
		return *this;
	}

	inline F_REAL operator/(const F_REAL &b) const
	{
		return is_zero ? F_REAL(0) : F_REAL(x / b.x);
	}

	inline F_REAL & operator/=(const F_REAL &b)
	{
		if (!is_zero)
			x /= b.x;
		return *this;
	}

	inline F_REAL & operator=(const F_REAL &b)
	{
		x = b.x;
		valid = b.valid;
		is_zero = b.is_zero;
		return *this;
	}
*/
};

template <typename F>
class FUNCTIONAL_IVP_SOLVER_RECURSIVE : public FUNCTIONAL_object<std::vector<TM>,unsigned int > {
public:
	const F _flow;
	TM initial_value;
	std::vector<std::vector<std::vector<F_REAL>>> _a;

	const unsigned int _dimension;
	const bool iv_is_zero;

	/*Computes a_{\nu,n}^{(i+1)}*/
	TM a37part2(unsigned nu, unsigned n, unsigned i)
	{
		TM sum=REAL(0);
		for (unsigned j=0; j <= n; j++)
			sum += a(nu,1,j) * a(nu,i,n-j);
		return sum;
	}

	template <typename V, typename W>
	inline TM combination_recursive4(
		unsigned j, unsigned d_m_1, const V &i, const W &idx, unsigned l
	) {
		if (j == d_m_1) {
			unsigned k = l;
			return a(i[j],idx[i[j]],(iv_is_zero ? idx[i[j]] : 0) + k);
		} else {
			TM r = REAL(0);
			for (unsigned k=0; k<=l; k++) {
				r += a(i[j],idx[i[j]],(iv_is_zero ? idx[i[j]] : 0) + k)
				   * combination_recursive4(j+1, d_m_1, i, idx, l-k);
			}
			return r;
		}
	}


	TM auto_ivp_38(unsigned int nu, int l)
	{
		TM sum = REAL(0);

		dbg("auto_ivp(nu=%u, l=%u)\n", nu, l);
		for (typename F::iterator_type idx = _flow.iterator(nu, l); idx; ++idx) {
			unsigned i;

			if (_flow.mu() >= 0) {
				for (i=0; i<=_dimension; i++)
					if (_flow.mu() < idx[i]) /* XXX: problem for FLOW_PROXY when we get mu dependant on nu */
						break;
				if (i <= _dimension) {
					for (unsigned j=0; j<_dimension; j++)
						dbg("i_%u = %d, ", j+1, idx[j]);
					dbg("skipped: idx[%u] = %d > %d = mu\n",
						i, idx[i], _flow.mu());
					continue;
				}
			}

			int l_ = l - idx[_dimension];
			if (iv_is_zero)
				for (unsigned j=0; j<_dimension; j++)
					l_ -= idx[j];

			for (unsigned j=0; j<_dimension; j++)
				dbg("i_%u = %d, ", j+1, idx[j]);
			dbg("k = %d, l_ = %d", idx[_dimension],  l_);

			if (l_ < 0) {
				dbg("skipped: l_ < 0\n");
				continue;
			}
			dbg("\n");

			REAL c = _flow(nu, idx); /* TODO: TM? */

			const typename F::I_type &I = *idx;
			if (I.ni_gt0 > 0) {
				sum += combination_recursive4(0, I.ni_gt0-1, I.idx_i_gt0, I.ik, l_) * c;
			} else if (l_ == 0) {
				sum += c;
			} else {
				/* some a_{\nu,n}^{(i)} is 0 whereas n >= 1:
				 * nothing to add */
			}
		}

		return sum *(REAL(1)/ (l+1));
	}

	/*a(nu,i,n) returns a_{\nu,n}^{(i)}*/
	const TM & a(const unsigned nu, const unsigned i, const unsigned n)
	{
		F_REAL &r = _a[nu][i][n];
		if (r.valid != true) {
			if (i == 0) {
				r = F_REAL(REAL((n == 0) ? 1 : 0));
			} else if (i == 1) {
				/* Assumption: n > 0 since for n=0, a_{\nu,n} is
				 * valid since it is the \nu-th component of the
				 * initial value */
				assert(n > 0);
				r.x = auto_ivp_38(nu,n-1);
				r.valid = true;
			} else {
				r.x = a37part2(nu,n,i-1);
				r.valid = true;
			}
		}
		return r.x;
	}

	FUNCTIONAL_IVP_SOLVER_RECURSIVE(
		const F &flow,
		const std::vector<TM> &w, bool iv_is_zero
	) : _flow(flow), _dimension(flow.dimension()), iv_is_zero(iv_is_zero) 
	{
		initial_value=w[0]; //just to get the Taylor model!
		_a.resize(_dimension);
		for (unsigned nu = 0; nu < _dimension; nu++) {
			_a[nu].resize(2);
			_a[nu][0].resize(1);
			_a[nu][0][0] = F_REAL(REAL(1));
			_a[nu][1].resize(1);
			_a[nu][1][0].x = w[nu];
			_a[nu][1][0].valid = true;
		}
	}

	~FUNCTIONAL_IVP_SOLVER_RECURSIVE() {}

	std::vector<TM> eval(const unsigned int & n)
	{
		unsigned max_power = (_flow.mu() < 0) ? n : _flow.mu();

		for (unsigned int nu = 0; nu < _dimension; nu++) {
			if (max_power >= _a[nu].size())
				_a[nu].resize(max_power + 1);
			for(unsigned int j = 0; j < _a[nu].size(); j++)
				if(_a[nu][j].size() < n+1)
					_a[nu][j].resize(n+1);
		}

		//Timer t;
		//t.start();

		std::vector<TM> result(_dimension);
		for (unsigned int nu = 0; nu < _dimension; nu++) {
			result[nu] = a(nu,1,n);
//			iRRAM_DEBUG0(1,{cerr << "a(" << nu << ",1," << n << ") = " << result[nu]() << "\n";});
		}

		//t.stop();
		//fprintf(stderr, "n: %u,\tt: %luµs   \r", n, t.t);
		//fflush(stderr);
		return result;
	}
};


template <typename F>
inline FUNCTION<std::vector<TM>,unsigned int > ivp_solver_recursive(
	const F &flow,
	const std::vector<TM> &w,
	bool iv_is_zero
) {
	return new FUNCTIONAL_IVP_SOLVER_RECURSIVE<F>(flow, w, iv_is_zero);
}



template <typename Flow>
static void print_iterator(const Flow &F)
{
	cout << "## FLow under consideration:\n";
	cout << "## dimension: " << F.dimension() << "\n";
	for (unsigned nu=0; nu<F.dimension(); nu++) {
		auto it = F.iterator(nu,F.mu());
		cout << "## it ? " << (it ? "true" : "false") << "\n";
		for (; it; ++it) {
			unsigned i_ne0 = 0;
			cout << "## c_{" << nu << "," << it[it.size()-1];
			for (unsigned j=0; j<it.size()-1; j++) {
				cout << "," << it[j];
				if (it[j])
					i_ne0++;
			}
			REAL out = F(nu,it);
			cout << "} = " << out << ", i != 0: " << i_ne0 << "\n";
		}
	}
	for (unsigned nu=0; nu<F.dimension(); nu++)
		cout << "#RHS     : " << F.poly[nu] << "\n";
	cout << "## autonomous: " << (F.is_autonomous() ? "true" : "false") << "\n";
	assert(F.is_autonomous());//Jacobian not correct for non-autonomous!
}

struct Smallstep_Control {
	RATIONAL small_delta_t;
	unsigned n_smallsteps;
	bool n_not_delta;

	explicit Smallstep_Control(const char *init)
	{
		n_not_delta = !strchr(init, '/');
		if (n_not_delta) {
			n_smallsteps = atoi(init);
		} else {
			small_delta_t = RATIONAL(init);
		}
	}
};


static FUNCTION<std::vector<TM>,REAL> bigstep(
	const std::vector<TM> &w, const POLYNOMIAL_FLOW &F,
	const REAL &R2, const REAL &M2
) {
	FUNCTION<std::vector<TM>,unsigned int> a;
	a=ivp_solver_recursive(F, w, false);
	return taylor_sum(a, R2, M2,0);
}



struct Input {
	int p;
	Smallstep_Control ssteps;
	REAL delta, eps, final_t, R_scale;
	POLYNOMIAL_FLOW F;
	std::vector<REAL> w;
	int step_control_alg;
	int R_steps;
	int sweepto;
	int prcdiff;
};

template <bool autonomous>
void plot_output(const Input &in)
{
	FUNCTION<std::vector<TM>,REAL> taylor;

	POLYNOMIAL_FLOW F = in.F;
	std::vector<TM> w(in.w.size());
	TM::set_default_sweep(in.sweepto);
	TM::set_prec_diff(in.prcdiff);
	for (unsigned i=0;i<in.w.size(); i++){
		w[i]=TM(in.w[i]);
		w[i].round();
	}


	DYADIC current_t = 0;
	DYADIC small_t = 0;

	Timer t;
	t.start();
	REAL eps_local=in.eps;

	int it=0;

	unsigned bigsteps;

	REAL R2, M2;
	REAL Rs,Ro,Lo; //informational parameters
	DYADIC delta_t;
	sizetype max_err;

	for (bigsteps = 0; positive(in.final_t-current_t,-20); bigsteps++) {

		if (in.step_control_alg == 0 || bigsteps % in.step_control_alg == 0) {
			TM::polish(w);
		//	for (unsigned i=0;i<in.w.size(); i++) w[i].check4();
		//	for (unsigned i=0;i<in.w.size(); i++) w[i]=TM::round0(w[i]);
		//	for (unsigned i=0;i<in.w.size(); i++) w[i]=TM(w[i].to_real(),in.sweepto);
		//	for (unsigned i=0;i<in.w.size(); i++) w[i]=TM::round2(w[i]);
		} else {
			for (unsigned i=0;i<in.w.size(); i++) w[i].round();
		//		for (unsigned i=0;i<in.w.size(); i++) w[i]=TM::round2(w[i]);
		//		for (unsigned i=0;i<in.w.size(); i++) w[i]=TM::round2(w[i]);
		}

		t.stop();
		cout << "# --------------------------------------------\n";
		cout << "# big step " << bigsteps
		     << ", last max coeff: " << iRRAM_DEBUG_last_taylor_coeff
		     << ", " << (t.t / 1e3) << "ms\n";
		t.t = 0;
		t.start();


		if (it++ % in.R_steps == 0)
		{
			stiff code(-5);
			std::vector<REAL> wr(w.size());
			for (unsigned i=0;i<w.size();i++)
				wr[i] = w[i].to_real();
			FUNCTION<std::vector<TM>,unsigned int> a
				= ivp_solver_recursive(F, w, false);
			F.get_RM3(w, 0, in.delta, in.R_scale, eps_local,
			          R2, M2, Lo, Rs, Ro, a);
		}


		delta_t = approx(R2 * in.R_scale, size(R2)-20)
		        - scale(DYADIC(1), size(R2)-20);

		{
			sizetype err;
			w[0].to_real().geterror(err);
			max_err = err;
			for (unsigned k=1; k<w.size(); k++) {
				w[k].to_real().geterror(err);
				if (sizetype_less(max_err, err))
					max_err = err;
			}
			while (max_err.mantissa > 1024) {
				max_err.mantissa >>= 1;
				max_err.exponent++;
			}
		}

		if (in.step_control_alg > 0) {
			FUNCTION<std::vector<TM>,unsigned int> a
				= ivp_solver_recursive(F, w, false);
			taylor = taylor_sum(a, R2, M2,0);
		} else {

			auto bs = std::bind(bigstep, std::placeholders::_1,
			                    std::cref(F), std::cref(Ro), std::cref(M2));

			std::function<FUNCTION<std::vector<TM>,REAL>(const std::vector<TM> &)> bsf = bs;
			auto on_domain = from_value<LAZY_BOOLEAN,std::vector<TM>>(true);


			auto lip_alg = [Lo](const std::vector<TM>& w,const REAL& t) {
				return Lo;
			};
			std::function<REAL(const std::vector<TM> &,const REAL &)> lipf = lip_alg;

			taylor = lipschitzify(from_algorithm(bsf),
			                      from_algorithm(lipf),
			                      on_domain,
			                      w);
		} // ! in.step_control_alg

		if (in.ssteps.n_not_delta) {
			for (unsigned j=0; j<in.ssteps.n_smallsteps; j++) {
				REAL t = (delta_t / (int)in.ssteps.n_smallsteps) * (int)j;
				std::vector<TM> y = taylor(t);
				cout << "taylor( " << (t+current_t) << " ) = ( ";
				cout << y[0].to_real();
				for (unsigned k=1; k<y.size(); k++)
					cout << " , " << y[k].to_real();
				cout << "\n" << std::flush;
			}
		} else {
			for (; small_t < current_t + delta_t; small_t = small_t + REAL(in.ssteps.small_delta_t).as_DYADIC()) {
				std::vector<TM> y = taylor(REAL(small_t) - current_t);
				cout << "taylor( " << small_t << " ) = ( ";
				cout << y[0].to_real();
				for (unsigned k=1; k<y.size(); k++)
					cout << " , " << y[k].to_real();
				cout << "\n" << std::flush;
			}
		}

		REAL old_t = current_t;
		{
			DYADIC_precision dyadic_prec(size(R2)-20);
			current_t = current_t + delta_t;
		}

		std::vector<TM> w_old = w;

		w = taylor(current_t-old_t);

		if (it % in.R_steps == 0) {
			cerr << "# t = " <<swrite(old_t,15,iRRAM_float_show)
			     << " (R2,M2,L,Ro,Rs,delta_t) = ( "
			     << swrite(R2,20,iRRAM_float_show) << ", "
			     << swrite(M2,20,iRRAM_float_show) << ", "
			     << swrite(Lo,20,iRRAM_float_show) << ", "
			     << swrite(Ro,20,iRRAM_float_show) << ", "
			     << swrite(Rs,20,iRRAM_float_show) << ", "
			     << swrite(delta_t,20,iRRAM_float_show)
			     << ")\n";
			for (int i=0; i<w.size(); i++)
				cerr << " w_old[" << i << "]= " << w_old[i] << "\n";
			for (int i=0; i<w.size(); i++)
				cerr << " w_new[" << i << "]= " << w[i] << "\n";
// 		    for (int i=0; i< w.size();i++) 
// 		    for (int j=0; j< w.size();j++) 
// 		      for (int k=0; k< w.size();k++)
// 		      for (int l=0; l< w.size();l++)
// 
// 		         cerr << " di(k)/dj(l)[" << i<<","<<j<<","<<k<<","<<l<<"]= "<<w_old[i].c[k].ci/w[j].c[l].ci;
// 		    cerr << "\n";
		}

		if (!autonomous && !F.is_autonomous()) {
			F = POLYNOMIAL_FLOW(F, delta_t);
			print_iterator(F);
		}

		fprintf(stderr, "bs %d, prec %d*2^(%d), max_coeff %d   \n",
			bigsteps, max_err.mantissa, max_err.exponent,
			iRRAM_DEBUG_last_taylor_coeff);
		fflush(stderr);
	}
	R2 = R2 - delta_t;

	t.stop();
	cout << "# no. bigsteps: " << bigsteps << ", " << (t.t / 1e3) << "ms   \n";
}

__attribute__((format(printf,2,3)))
static void die(int exitval, const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	exit(exitval);
}

static const char usage[] =
"usage: echo <p> { <small-steps> | <delta_t_num>/<delta_t_den> } <delta> <eps>"
	" <R-scale> <x> <path/to/coeffs> <iv1> ... <ivd> <step-size-alg> <R_steps> <sweepto>  <prcdiff>| PIVP [-iRRAM-OPTS]\n"
"\n"
"  p               precision of output values in decimal places\n"
"  small-steps     # of steps to divide each big-step into\n"
"  delta_t_*       fraction defining the evenly spaced step-size for printing results\n"
"  delta           t-radius for estimation of (R,M)\n"
"  eps             inital solution-space-radius for estimation of (R,M)\n"
"  R_scale         proceed with each big-step only (<R-scale>)*R in t\n"
"  x               final t to stop the iteration at (or around)\n"
"  path/to/coeffs  coefficients of the flow function in the format '<d> (<nu> <k> <i1> ... <id>)*'\n"
"  iv*             initial values at t=0\n"
"  step-size-alg   0: Lipschitz, >0: Taylor model with recentering\n"
"  R_steps         recompute (R,M) only every R_steps (i.e decrease R to R*(1-R_scale) )\n"
"  sweepto         0: traditional, 1: to left monomial\n"
"  prcdiff         #steps to increase the precision of the TM centers\n"
;

__attribute__((format(printf,1,2)))
static void input_error(const char *which, ...)
{
	va_list ap;
	va_start(ap, which);
	fprintf(stderr, "input error: ");
	vfprintf(stderr, which, ap);
	va_end(ap);
	die(1, " invalid\n\n%s", usage);
}

static Input read_input()
{
	int p;
	std::string ssteps;
	REAL delta, eps, final_x, alpha, R_scale;
	std::string fname;
	int step_control_alg;
	int R_steps;
	int sweepto;
	int prcdiff;
	cout << "# input: p, { small-steps | delta_t_num/delta_t_den }, delta, eps, R-scale, x, fname, w1, ..., wd, sc, R_steps, sw, prc\n";
	cout << "# fname file format: <d> (<nu> <k> <i1> ... <id>)*\n";

	if (!(cin >> p)) input_error("<p>");
	if (!(cin >> ssteps)) input_error("{ <smalls-steps> | <delta_t_num>/<delta_t_den> }");
	if (!(cin >> delta)) input_error("<delta>");
	if (!(cin >> eps)) input_error("<eps>");
	if (!(cin >> R_scale)) input_error("<R-scale>");
	if (!(cin >> final_x)) input_error("<x>");

	if (!(cin >> fname)) input_error("<path/to/coeffs>");

	POLYNOMIAL_FLOW F = read_poly_flow(fname.c_str());
	std::vector<REAL> w(F.dimension());
	unsigned i;
	for (i=0; i<F.dimension(); i++) {
		std::string s;
		if (!(cin >> s) || !s.length()) input_error("<iv%d>", i);
		w[i] = parse_REAL(s.c_str());
	}
	if (!(cin >>  step_control_alg)) input_error("<step control parameter>");
	if (!(cin >>  R_steps)) input_error("<R_steps>");
	if (!(cin >>  sweepto)) input_error("<sweepto>");
	if (!(cin >>  prcdiff)) input_error("<prcdiff>");

	return {
		p, Smallstep_Control(ssteps.c_str()), delta, eps, final_x,
		R_scale, F, w, step_control_alg, R_steps, sweepto, prcdiff
	};
}

void compute()
{
  
	Input in = read_input();

	cout << setRwidth(in.p);

	print_iterator(in.F);

	if (in.F.is_autonomous())
		plot_output<true>(in);
	else
		plot_output<false>(in);
}

