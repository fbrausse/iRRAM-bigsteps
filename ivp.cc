#include "iRRAM.h"
#include <vector>
#include <list>
#include <cassert>
#include "iRRAM/limit_templates.h"
#include "extension-vector-2.h"

#include <sys/time.h>
#include <cstdarg>

#ifndef METHOD_PICARD
# define METHOD_PICARD	0
#endif

#define dbg(...)	iRRAM_DEBUG2(2,__VA_ARGS__)

unsigned int max_taylor_coeff=0; // only statistical purpose!
unsigned int last_taylor_coeff=0; // only statistical purpose!


using namespace iRRAM;
using std::vector;
using std::pair;

//********************************************************************************
// Small helper functions
// not yet part of the published version of the iRRAM
//********************************************************************************

//********************************************************************************
// Absolute value of vector in Euclidean space
//********************************************************************************

REAL abs(const std::vector<REAL>& x)
{
	unsigned int n=x.size();
	REAL sqrsum=0;
	for (unsigned i=0;i<n;i++) {
		sqrsum += square(x[i]);
	}
	return sqrt(sqrsum);
}
//********************************************************************************
//********************************************************************************
//********************************************************************************

/****************************************************************************************/
/* summation of Taylor sequences							*/
/* The constructor needs a sequences of numbers with a radius of convergence 		*/
/* larger(!) than the given value 'radius'.						*/
/* 'bound' has to be an upper bound for the absolute value of the sum function 		*/
/* on a circle with the given 'radius' value						*/
/****************************************************************************************/

class FUNCTIONAL_taylor_sum :public FUNCTIONAL_object<REAL,REAL> {

/***********************************************************/
/* local data: coefficients (as a function),  radius and corresponding bound */

	FUNCTION<REAL,unsigned int> _coeff;
	REAL _radius;
	REAL _bound;

public:
	FUNCTIONAL_taylor_sum(
			FUNCTION<REAL,unsigned int> coeff,
			const REAL& radius,
			const REAL& bound
	) {
		_coeff=coeff;
		_radius=radius;
		_bound=bound;
	}

	~FUNCTIONAL_taylor_sum() {}

/****************************************************************************************/
/* Evaluation:										*/
/* Compute coefficients until the truncation error term is smaller than actual_prec	*/
/* or until the accumulated error is already of the order of the truncation error	*/
/* (whichever comes first)								*/
/****************************************************************************************/

	REAL eval(const REAL &x)
	{
		ITERATION_STACK SAVED_STACK;
		ACTUAL_STACK.inlimit+=1;
		REAL sum=0;
		REAL best=0;
		REAL factor=1;
		REAL error=_bound*_radius/(_radius-abs(x));
		REAL errorfactor=abs(x)/_radius;
		iRRAM_DEBUG0(2,{cerr << "FUNCTIONAL_taylor_sum starting with precision "<<ACTUAL_STACK.actual_prec
			<< " at ratio "<< errorfactor.vsize.mantissa*pow(2,errorfactor.vsize.exponent)<<"\n";});

		sizetype sum_error,trunc_error,best_error;

		sizetype_add(best_error,error.vsize,error.error);
		int best_index=-1;

		for (unsigned int i=0;;i++){
			sum= sum + _coeff(i)*factor;
			error= error*errorfactor;

			sum.geterror(sum_error);
			sizetype_add(trunc_error,error.vsize,error.error);

			sizetype local_error;
			sizetype_add(local_error,trunc_error,sum_error);
			if (sizetype_less(local_error,best_error)) {
				best=sum;
				best_error=local_error;
				best.seterror(best_error);
				best_index=i;
			}
			if (trunc_error.exponent<ACTUAL_STACK.actual_prec ||
			     sizetype_less(trunc_error,sum_error)) {
				iRRAM_DEBUG0(2,{cerr << "FUNCTIONAL_taylor_sum: stop at "
						<<i<< " with best at "<<best_index<<"\n";});
				if (i>max_taylor_coeff) max_taylor_coeff=i; // only statistical purpose!
				last_taylor_coeff = i;
				break;
			}
			factor=factor*x;
		}
		return best;
	}
};

/********************************************************/
/* corresponding constructor for FUNCTION		*/
/********************************************************/
inline FUNCTION<REAL,REAL> taylor_sum (
		FUNCTION<REAL,unsigned int> coeff,
		const REAL& radius,
		const REAL& bound
) {
	return new FUNCTIONAL_taylor_sum(coeff,radius,bound);
}

//********************************************************************************
// Class for summation of a vector of Taylor series
// 'bound' must be valid for 'radius' in each dimension
//********************************************************************************

class FUNCTIONAL_vector_taylor_sum :public FUNCTIONAL_object<std::vector<REAL>,REAL > {

	REAL _radius;
	REAL _bound;
	FUNCTION<std::vector<REAL>,unsigned int > _coeff;

public:
	FUNCTIONAL_vector_taylor_sum(
			FUNCTION<std::vector<REAL>,unsigned int > coeff,
			const REAL& radius,
			const REAL& bound
	) {
		_coeff=coeff;
		_radius=radius;
		_bound=bound;
	}

	~FUNCTIONAL_vector_taylor_sum() {}

	std::vector<REAL> eval(const REAL& t)
	{
		std::vector<REAL> result=_coeff(0); // just to get the dimension
		for (unsigned int nu=0;nu<result.size();nu++){
			FUNCTION<REAL,unsigned int>  coeff_nu=projection(_coeff,nu);
			FUNCTION<REAL,REAL> f_nu=taylor_sum(coeff_nu,_radius,_bound);
			result[nu]=f_nu(t);
		}
		return result;
	}
};

/********************************************************/
/* corresponding constructor for FUNCTION		*/
/********************************************************/
inline FUNCTION<std::vector<REAL>,REAL > taylor_sum (
		FUNCTION<std::vector<REAL>,unsigned int > coeff,
		const REAL& radius,
		const REAL& maximum
) {
	return new FUNCTIONAL_vector_taylor_sum(coeff,radius,maximum);
}

//********************************************************************************
//********************************************************************************


template <typename K>
class POLYNOMIAL {
	std::vector<K> c;
public:
	POLYNOMIAL() : c{0} {assert(c.size());}
	inline POLYNOMIAL(const K &c0) : c{c0} {assert(c.size());}
	inline explicit POLYNOMIAL(const std::vector<K> &c) : c(c)
	{
		if (c.size() == 0)
			this->c.push_back(K(0));
		assert(this->c.size());
	}
	template <typename K2>
	POLYNOMIAL(const POLYNOMIAL<K2> &p) : c(p.degree() + 1)
	{
		for (unsigned i=0; i<c.size(); i++)
			c[i] = p[i];
		assert(c.size());
	}

	void set_coeff(unsigned long i, const K &v)
	{
		c.resize(std::max(i+1, c.size()));
		c[i] = v;
	}

	template <typename K2>
	auto operator()(const K2 &x) const -> decltype(c[0]*x)
	{
		using KR = decltype(c[0]*x);
		KR r(0);
		K2 xi(1);
		for (const K &ci : c) {
			r  += ci * xi;
			xi *= x;
		}
		return r;
	}

	template <typename K2>
	auto operator()(const POLYNOMIAL<K2> &b) const
	-> decltype((*this)*(b))
	{
		using PR = decltype((*this) * b);
		PR r;
		POLYNOMIAL<K2> bi = K2(1);
		for (const K &ci : c) {
			r  += POLYNOMIAL(ci) * bi;
			bi *= b;
		}
		return r;
	}

	/* beware: returns just upper bound, unknown whether c[degree()] == 0 */
	inline unsigned degree() const { return c.size()-1; }

	inline const K & operator[](unsigned i) const { return c[i]; }

	template <typename K2>
	POLYNOMIAL & operator+=(const POLYNOMIAL<K2> &b)
	{
		c.resize(1+std::max(degree(), b.degree()));
		for (unsigned i=0; i<=b.degree(); i++)
			c[i] += b[i];
		return *this;
	}

	template <typename K2>
	inline auto operator+(const POLYNOMIAL<K2> &b) const
	-> POLYNOMIAL<decltype(c[0]+b[0])>
	{
		return POLYNOMIAL<decltype(c[0]+b[0])>(*this) += b;
	}

	template <typename K2>
	POLYNOMIAL & operator-=(const POLYNOMIAL<K2> &b)
	{
		c.resize(1+std::max(degree(), b.degree()));
		for (unsigned i=0; i<=b.degree(); i++)
			c[i] = c[i] - b[i];
		return *this;
	}

	template <typename K2>
	inline auto operator-(const POLYNOMIAL<K2> &b) const
	-> POLYNOMIAL<decltype(c[0]-b[0])>
	{
		return POLYNOMIAL<decltype(c[0]-b[0])>(*this) -= b;
	}

	POLYNOMIAL operator-() const
	{
		POLYNOMIAL r = *this;
		for (K &ci : r.c)
			ci = -ci;
		return r;
	}

	template <typename K2>
	auto operator*(const POLYNOMIAL<K2> &b) const
	-> POLYNOMIAL<decltype(c[0]*b[0])>
	{
		using KR = decltype(c[0]*b[0]);
		std::vector<KR> r(degree() + b.degree() + 1);
		for (unsigned i=0; i<=degree(); i++)
			for (unsigned j=0; j<=b.degree(); j++)
				r[i+j] += c[i] * b[j];
		return POLYNOMIAL<KR>(r);
	}

	template <typename K2>
	inline POLYNOMIAL & operator*=(const POLYNOMIAL<K2> &b)
	{
		return *this = *this * b;
	}

	POLYNOMIAL primitive(const K &C = 0) const
	{
		POLYNOMIAL b(std::vector<K>(c.size() + 1));
		b.c[0] = C;
		b.c[1] = c[0];
		for (unsigned i=1; i<c.size(); i++)
			b.c[i+1] = c[i] / (int)(i+1);
		return b;
	}

	POLYNOMIAL derivative() const
	{
		POLYNOMIAL b(std::vector<K>(c.size() - 1));
		for (unsigned i=1; i<c.size(); i++)
			b.c[i-1] = c[i] * (int)i;
		return b;
	}

	friend POLYNOMIAL mul_imod(const POLYNOMIAL &a, const POLYNOMIAL &b, unsigned mod_degree)
	{
		std::vector<K> r(std::min(mod_degree, degree() + b.degree() + 1));
		unsigned amin = std::min(mod_degree, a.degree()+1);
		for (unsigned i=0; i<amin; i++) {
			unsigned bmin = std::min(mod_degree-i, b.degree()+1);
			for (unsigned j=0; j<bmin; j++)
				r[i+j] += a[i] * b[j];
		}
		return POLYNOMIAL(r);
	}

	friend POLYNOMIAL imod(const POLYNOMIAL &p, unsigned degree)
	{
		if (degree > p.degree())
			return p;
		decltype(p.c) c(degree);
		for (unsigned i=0; i<degree; i++)
			c[i] = p.c[i];
		return POLYNOMIAL(c);
	}

	friend orstream & operator<<(orstream &o, const POLYNOMIAL &p)
	{
		o << "(" << p[0] << ")";
		for (unsigned i=1; i<=p.degree(); i++)
			o << "+(" << p[i] << ")*x^" << i;
		return o;
	}

	friend POLYNOMIAL abs(const POLYNOMIAL &p)
	{
		POLYNOMIAL r = p;
		for (K &ci : r.c)
			ci = abs(ci);
		return r;
	}
};


class POLYNOMIAL2 {
public:
	struct I {
		std::vector<unsigned> i;
		REAL ci;
	};

	explicit POLYNOMIAL2(unsigned d, unsigned n_coeffs = 0) : mu(0), d(d)
	{
		c.reserve(n_coeffs);
	}

	unsigned degree() const { return mu; }

	POLYNOMIAL2 & operator+=(const POLYNOMIAL2 &p)
	{
		for (const I &i : p.c) {
			for (I &j : c) {
				if (i.i == j.i) {
					j.ci += i.ci;
					goto ok;
				}
			}
			c.push_back(i);
ok:;
		}
		mu = std::max(mu, p.mu);
		return *this;
	}

	const std::vector<I> & coefficients() const { return c; }

	void add_coeff(const I &i)
	{
		assert(i.i.size() == d);
		c.push_back(i);
		for (const unsigned &ij : i.i)
			mu = std::max(mu, ij);
	}

	POLYNOMIAL2 & operator*=(const REAL &r)
	{
		for (I &l : c)
			l.ci *= r;
		return *this;
	}
	POLYNOMIAL2 operator*(const REAL &r) const { return POLYNOMIAL2(*this) *= r; }

	POLYNOMIAL2 & operator*=(const I &k)
	{
		for (I &l : c) {
			assert(l.i.size() == k.i.size());
			for (unsigned j = 0; j < l.i.size(); j++) {
				l.i[j] += k.i[j];
				if (l.i[j] > mu) mu = l.i[j];
			}
			l.ci *= k.ci;
		}
		return *this;
	}
	POLYNOMIAL2 operator*(const I &k) const { return POLYNOMIAL2(*this) *= k; }

	POLYNOMIAL2 & operator*=(const POLYNOMIAL2 &p)
	{
		for (const I &i : p.c) {
			*this *= i;
		}
		return *this;
	}
	POLYNOMIAL2 operator*(const POLYNOMIAL2 &p) const { return POLYNOMIAL2(*this) *= p; }

	/* partial derivative */
	POLYNOMIAL2 derivative(unsigned component) const
	{
		POLYNOMIAL2 r(d);
		for (const I &i : c) {
			if (!i.i[component])
				continue;
			I j = i;
			j.ci *= j.i[component]--;
			r.add_coeff(j);
		}
		return r;
	}

	template <typename K2>
	auto operator()(const std::vector<K2> &x) -> decltype(REAL(0) * x[0])
	{
		using KR = decltype(REAL(0) * x[0]);
		assert(x.size() == d);
		KR r = 0;
		for (const I &i : c) {
			K2 s = 1;
			for (unsigned j=0; j<d; j++)
				s *= power(x[j], i.i[j]);
			r += i.ci * s;
		}
		return r;
	}

	friend orstream & operator<<(orstream &o, const POLYNOMIAL2 &p)
	{
		bool first = true;
		for (const I &i : p.c) {
			if (!first)
				o << "+";
			o << "(" << i.ci << ")";
			for (unsigned k=0; k<p.d; k++)
				o << "*x" << k << "^" << i.i[k];
			first = false;
		}
		return o;
	}

private:
	unsigned mu;
	unsigned d;
	std::vector<I> c;
};

struct VI {
	std::vector<unsigned> ik;
	std::vector<unsigned> idx_i_gt0;
	unsigned ni_gt0;
	REAL c;

	VI(const REAL &c, const std::vector<unsigned> &ik)
	: ik(ik), c(c)
	{
		for (unsigned j=0; j<ik.size()-1; j++)
			if (ik[j] > 0)
				idx_i_gt0.push_back(j);
		ni_gt0 = idx_i_gt0.size();
	}

	VI(const POLYNOMIAL2::I &pi) : VI(pi.ci, pi.i) {}

	operator POLYNOMIAL2::I() const
	{
		return { ik, c };
	}
};

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
	{}

	POLYNOMIAL_FLOW(const POLYNOMIAL_FLOW &old, const REAL &told)
	: _d(old._d), _mu(0), _autonomous(old.is_autonomous()), c(old._d)
	{
		unsigned d = _d;
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
		for (unsigned j=0; j<=_d; j++)
			_mu = std::max(_mu, coeff.ik[j]);
		if (coeff.ik[_d] != 0)
			_autonomous = false;
	}

	inline unsigned dimension() const { return _d; }
	inline int      mu()        const { return _mu; }
	inline bool     is_autonomous() const { return _autonomous; }

	iterator_type iterator(unsigned nu, unsigned l) const
	{
		return Iterator(*this, nu);
	}

	REAL operator()(unsigned nu, const iterator_type &idx) const
	{
		return idx->c;
	}

	REAL UF(std::vector<REAL> w, REAL t0, const REAL &delta, const REAL &eps) const
	{
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

	void get_RM(
		const std::vector<REAL> &w,
		const REAL &t0,
		const REAL &delta,
		const REAL &eps,
		REAL &R,
		REAL &M
	) const {
		REAL _UF = UF(w,t0,delta,eps);
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
		M = abs_w + R*UF(w,t0,R,eps);
	}
#if 0
	REAL f_abs(unsigned nu, REAL x) const
	{
		REAL r(0);
		for (Iterator it = iterator(nu, _mu); it; ++it) {
			unsigned sum_i = 0;
			for (unsigned j=0; j<dimension(); j++)
				sum_i += it[j];
			r += abs((*this)(nu, it)) * power(x, sum_i);
		}
		return r;
	}

	REAL f_max(const REAL &x) const
	{
		REAL m = f_abs(0, x);
		for (unsigned nu=1; nu<dimension(); nu++)
			m = maximum(m, f_abs(nu, x));
		return m;
	}
#endif
	REAL f_abs(unsigned nu, const std::vector<REAL> &w_abs, const REAL &s, const REAL &t) const
	{
		REAL r(0);
		for (Iterator it = iterator(nu, _mu); it; ++it) {
			REAL m = abs((*this)(nu, it));
			if (it[dimension()] != 0)
				m *= power(t, it[dimension()]);
			for (unsigned j=0; j<it->ni_gt0; j++)
				m *= power(w_abs[it->idx_i_gt0[j]] + s, it[it->idx_i_gt0[j]]);
			r += m;
		}
		return r;
	}

	REAL f_max(const std::vector<REAL> &w, const REAL &s, const REAL &t) const
	{
		REAL m = f_abs(0, w, s, t);
		for (unsigned nu=1; nu<dimension(); nu++)
			m = maximum(m, f_abs(nu, w, s, t));
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

private:
	unsigned _d;
	unsigned _mu;
	bool _autonomous;
	std::vector<std::vector<I_type>> c;
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
		return s;
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



template <typename F>
class FUNCTIONAL_ivp_solver_auto :public FUNCTIONAL_object<std::vector<REAL>,unsigned int >
{
	F _flow;

/* The vector a[nu][i] contains an initial segment of the Taylor coefficients
   for (y_nu)^i.
   As soon as taylor[nu][i][1] is computed, we will also compute 
   a[nu][j][i] for as many values i (<=j) as needed by the flow.
*/
	std::vector< std::vector< std::vector<REAL> > >a;

	unsigned int _dimension;
	const bool iv_is_zero;


	REAL a_vector_power(const std::vector<int> &n, const typename F::iterator_type &i) const {
		REAL ergebnis = REAL(1);
		for (unsigned int nu = 0; nu < _dimension; nu++) {
			ergebnis *= a[nu][i[nu]][n[nu]];
		}
		return ergebnis;
	}

	inline void combination_recursive(std::vector<int> &n, unsigned j, REAL &sum, const typename F::iterator_type &idx, unsigned l) const
	{
		if (j == n.size() -1) {
			n[j] = (iv_is_zero ? idx[j] : 0) + l;
			sum += a_vector_power(n, idx);
			return;
		}

		if (idx[j] == 0) {
			n[j] = 0;
			combination_recursive(n, j+1, sum, idx, l);
		} else {
			for (unsigned k=0; k <= l; k++) {
				n[j] = (iv_is_zero ? idx[j] : 0) + k;
				combination_recursive(n, j+1, sum, idx, l-k);
			}
		}
	}

	inline void combination_recursive2(const REAL &mul, unsigned j, unsigned d, REAL &sum, const typename F::iterator_type &idx, unsigned l) const
	{
		if (j == d) {
			if (l == 0)
				sum += mul;
			return;
		}

		if (idx[j] == 0) {
			combination_recursive2(mul, j+1, d, sum, idx, l);
		} else {
			for (unsigned k = 0; k <= l; k++)
				combination_recursive2(mul * a[j][idx[j]][(iv_is_zero ? idx[j] : 0)+k], j+1, d, sum, idx, l-k);
		}
	}

	inline void combination_recursive3(
		const REAL &mul,
		unsigned j,
		unsigned d,
		REAL &sum,
		const std::vector<int> &i,
		const typename F::iterator_type &idx,
		unsigned l
	) const {
		if (j == d-1) {
			unsigned k = l;
			sum += mul * a[i[j]][idx[i[j]]][(iv_is_zero ? idx[i[j]] : 0)+k];
		} else {
			for (unsigned k = 0; k <= l; k++) {
				combination_recursive3(mul * a[i[j]][idx[i[j]]]
				[(iv_is_zero ? idx[i[j]] : 0)+k], j+1, d, sum, i, idx, l-k);
			}
		}
	}

	REAL auto_ivp_38(unsigned int nu, int l) const
	{
		REAL sum(0);

		std::vector<int> n;
		std::vector<int> n_;
		n.resize(_dimension);
		n_.resize(_dimension);

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

			REAL c = _flow(nu, idx);
			if (l_ == 0) {
				for (unsigned j=0; j<_dimension; j++)
					n[j] = iv_is_zero ? idx[j] : 0;
				sum += c * a_vector_power(n, idx);
				continue;
			}

			REAL sum_a_ni = 0;
#if 0
			combination_recursive(n_, 0, sum_a_ni, idx, l_);
#elif 0
			combination_recursive2(REAL(1), 0, _dimension, sum_a_ni, idx, l_);
#else
			unsigned j=0;
			for (unsigned i=0; i<_dimension; i++)
				if (idx[i] > 0)
					n[j++] = i;
			if (j == 0)
				sum_a_ni = 0;
			else
				combination_recursive3(REAL(1), 0, j, sum_a_ni, n, idx, l_);
#endif
			sum += c * sum_a_ni;
		}

		return sum/(l+1);
	}

	/* berechnet a_{\nu,n}^{(i+1)} */
	REAL a37part2(unsigned nu, unsigned n, unsigned i)
	{
		REAL sum(0);
		for (unsigned j=0; j<=n; j++)
			sum += a[nu][1][j] * a[nu][i][n-j];
		return sum;
	}

public:
	FUNCTIONAL_ivp_solver_auto(const F &flow, const std::vector<REAL> &w, bool iv_is_zero)
	: _flow(flow), iv_is_zero(iv_is_zero)
	{
		_dimension = _flow.dimension();

		a.resize(_dimension);
		for (unsigned int nu = 0; nu < _dimension; nu++) {
			a[nu].resize(2);
			a[nu][0].resize(1);
			a[nu][0][0] = REAL(1);
			a[nu][1].resize(1);
			a[nu][1][0] = w[nu];
		}
	}

	~FUNCTIONAL_ivp_solver_auto() {}

	const vector<vector<vector<REAL>>> & get_a() const { return a; }

	/* |{(n_1,...,n_d): n_k>=i_k, \sum_k n_k=l}| = (d+l'-1 choose d-1)
	 * mit l' = l - \sum_k i_k
	 *
	 * gibt die Anzahl Summanden für a_{\nu,l+1} für ein (i_1,...,i_d) an.
	 */

	std::vector<REAL> eval(const unsigned int & n)
	{
		unsigned l_old = a[0][0].size()-1;
		unsigned max_power = (_flow.mu() < 0) ? n : _flow.mu();

		for (unsigned int nu = 0; nu < _dimension; nu++) {
			if (max_power >= a[nu].size())
				a[nu].resize(max_power + 1);
			for(unsigned int j = 0; j < a[nu].size(); j++)
				if(a[nu][j].size() < n+1)
					a[nu][j].resize(n+1);
		}

		for (unsigned int l = l_old; l < n; l++) {/*
			fprintf(stderr,
				"l: %u, auto: %luµs, rest0: %luµs, rest1: %luµs   \r",
				l, t_auto.t, t_rest0.t, t_rest1.t);
			fflush(stderr);
			t_auto.t = t_rest0.t = t_rest1.t = 0;*/
			for (unsigned int nu = 0; nu < _dimension; nu++) {
				/* Um in der nächsten äußeren Iteration
				 * a_{\nu,l+1} berechnen zu können, benötigen
				 * wir a_{j,k}^{(i)} mit 1<=j<=d, 0<=k<=l und
				 *             _
				 *            /  l, falls µ < 0
				 * 0 <= i <= <
				 *            \_ µ, sonst.
				 *
				 * Berechne also in aktueller Iteration (l,nu)
				 * a_{nu,k}^{(i)} wie unten mit 'x' dargestellt.
				 *
				 * µ < 0:             bzw.  µ >= 0:
				 * ~~~~~~                   ~~~~~~~
				 *    |      l                 |      l
				 * ---+------+------> k     ---+-------------> k
				 *    |      x                 |      x
				 *    |      x                 |      x
				 *    |      x                 |      x
				 *  l +xxxxxxx                 |      x
				 *    |                      µ +      x
				 *    |                        |
				 *    v                        v
				 * 
				 *    i                        i
				 */

				unsigned int i, k;                  /* (3.7) Teil 2 */
				
				/* TODO: Vertauschen der Indizierung a[nu][i][k]
				 *       zu a[nu][k][i] wg. Zugriffsgeschwindigkeit
				 *       in Schleife for (i=2; i<=i_max; i++)
				 */

				//t_rest0.start();
				k = l;
				unsigned i_max = _flow.mu() >= 0 ? _flow.mu()
				                                 : std::min(max_power, l);
				for (i=2; i <= i_max; i++)
					a[nu][i][k] = a37part2(nu, k, i-1);
				//t_rest0.stop();

				//t_rest1.start();
				if (_flow.mu() < 0) {
					i = l;
					if (i <= max_power) {
						for (k=0; k< l; k++) /* "<" because we already have a[nu][i][l], see above */
							a[nu][i][k] = a37part2(nu, k, i-1);
					}
				}
				//t_rest1.stop();

				a[nu][0][l+1] = 0;                  /* (3.7) Teil 1 */
				//t_auto.start();
				a[nu][1][l+1] = auto_ivp_38(nu, l); /* (3.8) */
				//t_auto.stop();
			}
		}

		std::vector<REAL> result(_dimension);
		for (unsigned int nu = 0; nu < _dimension; nu++) {
			result[nu] = a[nu][1][n];
			iRRAM_DEBUG0(3,{cerr << "a(" << nu << ",1," << n << ") = " << result[nu] << "\n";});
		}
		return result;
	}
};

struct F_REAL {
	REAL x;
	bool valid;
	bool is_zero;

	inline F_REAL() : valid(false), is_zero(false) {}
	inline explicit F_REAL(int v) : x(v), valid(true), is_zero(v == 0) {}

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
class FUNCTIONAL_IVP_SOLVER_RECURSIVE : public FUNCTIONAL_object<std::vector<REAL>,unsigned int > {
public:
	const F _flow;
	std::vector< std::vector< std::vector<F_REAL> > > _a;

	const unsigned int _dimension;
	const bool iv_is_zero;

	/*Computes a_{\nu,n}^{(i+1)}*/
	REAL a37part2(unsigned nu, unsigned n, unsigned i)
	{
		REAL sum(0);
		for (unsigned j=0; j <= n; j++)
			sum += a(nu,1,j) * a(nu,i,n-j);
		return sum;
	}

	template <typename V, typename W>
	inline REAL combination_recursive4(
		unsigned j, unsigned d_m_1, const V &i, const W &idx, unsigned l
	) {
		if (j == d_m_1) {
			unsigned k = l;
			return a(i[j],idx[i[j]],(iv_is_zero ? idx[i[j]] : 0) + k);
		} else {
			REAL r = 0;
			for (unsigned k=0; k<=l; k++) {
				r += a(i[j],idx[i[j]],(iv_is_zero ? idx[i[j]] : 0) + k)
				   * combination_recursive4(j+1, d_m_1, i, idx, l-k);
			}
			return r;
		}
	}


	REAL auto_ivp_38(unsigned int nu, int l)
	{
		REAL sum(0);

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

			REAL c = _flow(nu, idx);

			const typename F::I_type &I = *idx;
			if (I.ni_gt0 > 0) {
				sum += c * combination_recursive4(0, I.ni_gt0-1, I.idx_i_gt0, I.ik, l_);
			} else if (l_ == 0) {
				sum += c;
			} else {
				/* some a_{\nu,n}^{(i)} is 0 whereas n >= 1:
				 * nothing to add */
			}
		}

		return sum / (l+1);
	}

	/*a(nu,i,n) returns a_{\nu,n}^{(i)}*/
	const REAL & a(const unsigned nu, const unsigned i, const unsigned n)
	{
		F_REAL &r = _a[nu][i][n];
		if (r.valid != true) {
			if (i == 0) {
				r = F_REAL((n == 0) ? 1 : 0);
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
		const std::vector<REAL> &w, bool iv_is_zero
	) : _flow(flow), _dimension(flow.dimension()), iv_is_zero(iv_is_zero) 
	{
		_a.resize(_dimension);
		for (unsigned nu = 0; nu < _dimension; nu++) {
			_a[nu].resize(2);
			_a[nu][0].resize(1);
			_a[nu][0][0] = F_REAL(1);
			_a[nu][1].resize(1);
			_a[nu][1][0].x = w[nu];
			_a[nu][1][0].valid = true;
			_a[nu][1][0].is_zero = iv_is_zero;
		}
	}

	~FUNCTIONAL_IVP_SOLVER_RECURSIVE() {}

	std::vector<REAL> eval(const unsigned int & n)
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

		std::vector<REAL> result(_dimension);
		for (unsigned int nu = 0; nu < _dimension; nu++) {
			result[nu] = a(nu,1,n);
			iRRAM_DEBUG0(3,{cerr << "a(" << nu << ",1," << n << ") = " << result[nu] << "\n";});
		}

		//t.stop();
		//fprintf(stderr, "n: %u,\tt: %luµs   \r", n, t.t);
		//fflush(stderr);

		return result;
	}
};

class FUNCTIONAL_IVP_SOLVER_PICARD : public FUNCTIONAL_object<std::vector<REAL>,unsigned int> {
public:
	const POLYNOMIAL_FLOW F;
	/* p[n][nu] */
	std::vector<std::vector<POLYNOMIAL<REAL>>> p;
	unsigned n = 0;

	const std::vector<REAL> w;
	const bool iv_is_zero;

	FUNCTIONAL_IVP_SOLVER_PICARD(
		const POLYNOMIAL_FLOW &F,
		const std::vector<REAL> &w,
		bool iv_is_zero
	) : F(F), p(1), w(w), iv_is_zero(iv_is_zero)
	{
		p[0].resize(F.dimension());
		for (unsigned nu=0; nu<F.dimension(); nu++)
			p[0][nu].set_coeff(0, w[nu]);
	}

	~FUNCTIONAL_IVP_SOLVER_PICARD() {}

	void step()
	{
		/* p(t) <- w + \int_0^t F(s,p(s)) mod s^n ds */
		REAL t0 = 0;
		unsigned d = F.dimension();
		p.resize(n+2);
		std::vector<POLYNOMIAL<REAL>> &r = p[n+1];
		r.resize(d);
		for (unsigned nu=0; nu<d; nu++) {
			for (auto it = F.iterator(nu, F.mu()); it; ++it) {
				POLYNOMIAL<REAL> p/*(F(nu,it))*/;
				p.set_coeff(it->ik[d], F(nu, it)); /* constant poly: c_{\nu,k,i} */
				for (unsigned xi : it->idx_i_gt0) {
					POLYNOMIAL<REAL> q;
					q.set_coeff(it->ik[xi], 1);
					p *= q(this->p[n][xi]);
				}
				r[nu] += imod(p, n+1);
			}
			POLYNOMIAL<REAL> q = r[nu].primitive();
			q += POLYNOMIAL<REAL>(w[nu] - q(t0));
			r[nu] = q;
		}
		n++;
	}

	std::vector<REAL> eval(const unsigned int &n)
	{
		std::vector<REAL> result(F.dimension());

		while (this->n <= n)
			step();

		for (unsigned nu=0; nu<F.dimension(); nu++)
			result[nu] = p[n+1][nu].degree() < n ? REAL(0) : p[n+1][nu][n];

		return result;
	}
};

template <typename F>
inline FUNCTION<std::vector<REAL>,unsigned int > ivp_solver_recursive(
	const F &flow,
	const std::vector<REAL> &w,
	bool iv_is_zero
) {
	return new FUNCTIONAL_IVP_SOLVER_RECURSIVE<F>(flow, w, iv_is_zero);
}

inline FUNCTION<std::vector<REAL>,unsigned int > ivp_solver_picard(
	const POLYNOMIAL_FLOW &flow,
	const std::vector<REAL> &w,
	bool iv_is_zero
) {
	return new FUNCTIONAL_IVP_SOLVER_PICARD(flow, w, iv_is_zero);
}

template <typename F>
inline FUNCTION<std::vector<REAL>,unsigned int > ivp_solver_auto(
	const F &flow,
	const std::vector<REAL> &w,
	bool iv_is_zero
) {
	return new FUNCTIONAL_ivp_solver_auto<F>(flow, w, iv_is_zero);
}


template <typename Flow>
static void print_iterator(const Flow &F)
{
	cout << "# dimension: " << F.dimension() << "\n";
	for (unsigned nu=0; nu<F.dimension(); nu++) {
		auto it = F.iterator(nu,F.mu());
		cout << "# it ? " << (it ? "true" : "false") << "\n";
		for (; it; ++it) {
			unsigned i_ne0 = 0;
			cout << "# c_{" << nu << "," << it[it.size()-1];
			for (unsigned j=0; j<it.size()-1; j++) {
				cout << "," << it[j];
				if (it[j])
					i_ne0++;
			}
			REAL out = F(nu,it);
			cout << "} = " << out << ", i != 0: " << i_ne0 << "\n";
		}
	}
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

#if 1
template <bool picard>
static FUNCTION<std::vector<REAL>,REAL> bigstep(
	const std::vector<REAL> &w, const POLYNOMIAL_FLOW &F,
	const REAL &R2, const REAL &M2
) {
	FUNCTION<std::vector<REAL>,unsigned int> a;

	a = picard ? ivp_solver_picard(F, w, false)
	           : ivp_solver_recursive(F, w, false);
	return taylor_sum(a, R2, M2);
}
#else
template <bool picard>
static std::vector<REAL> bigstep(
	const std::vector<REAL> &w, const POLYNOMIAL_FLOW &F,
	const REAL &R2, const REAL &M2,
	const DYADIC &current_t, DYADIC &delta_t,
	FUNCTION<std::vector<REAL>,REAL> &taylor, sizetype &max_err
) {
	FUNCTION<std::vector<REAL>,unsigned int> a;

	a = picard ? ivp_solver_picard(F, w, false)
	           : ivp_solver_recursive(F, w, false);
	taylor = taylor_sum(a, R2, M2);
	return taylor((current_t + delta_t) - REAL(current_t));
}
#endif

struct Input {
	int p;
	Smallstep_Control ssteps;
	REAL delta, eps, final_t, R_scale;
	POLYNOMIAL_FLOW F;
	std::vector<REAL> w;
	int step_control_alg;
};

template <bool autonomous,bool picard>
void plot_output(const Input &in)
{
	FUNCTION<std::vector<REAL>,REAL> taylor;

	const int cmp_p = -10;
	const int delta_t_p = -53;
	const DYADIC end_t = approx(in.final_t, cmp_p) + scale(INTEGER(2), cmp_p);

	POLYNOMIAL_FLOW F = in.F;
	std::vector<REAL> w = in.w;

	DYADIC_precision dyadic_prec(delta_t_p);

	DYADIC current_t = 0;
	DYADIC small_t = 0;

	Timer t;
	t.start();
	REAL eps_local=in.eps;
	sizetype err;

	unsigned bigsteps;
	// for (bigsteps = 0; !positive(REAL(current_t) - final_t, cmp_p); bigsteps++) {
	for (bigsteps = 0; current_t <= end_t /* final_t || current_t <=final_t*1.00001*/; bigsteps++) {
		REAL R, R_test, M, M_test;
		REAL R2, M2;
		DYADIC delta_t;
		sizetype max_err;

		t.stop();
		cout << "# --------------------------------------------\n";
		cout << "# big step " << bigsteps
		     << ", last max coeff: " << last_taylor_coeff
		     << ", " << (t.t / 1e3) << "ms\n";
		t.t = 0;
		t.start();

		F.get_RM(w, 0, in.delta, in.eps, R, M);
		cout << "#  "<<current_t << " (R ,M ) = ( " << R << ", " << M << ")\n";

		F.get_RM2(w, 0, in.delta, eps_local, R2, M2, in.step_control_alg);
		cout << "#  " << current_t << " (R2,M2) = ( " << R2 << ", " << M2 << ")\n";

		/* TODO: while (INTEGER(R2 * R_scale * 2**n)-1 <= 0) n++;
		 *  und  DYADIC_precision(n) */
		delta_t = approx(R2 * in.R_scale, delta_t_p) - scale(DYADIC(1), delta_t_p);
		cout << "# t = " << current_t << ", delta_t = " << delta_t << "\n";
		cout << "# w = (" << w[0];
		for (unsigned k=1; k<w.size(); k++)
			cout << ", " << w[k];
		w[0].geterror(err);
		cout << ") +/- (" << err.mantissa << "*2^(" << err.exponent << ")";
		max_err = err;
		for (unsigned k=1; k<w.size(); k++) {
			w[k].geterror(err);
			cout << ", " << err.mantissa << "*2^(" << err.exponent << ")";
			if (sizetype_less(max_err, err))
				max_err = err;
		}
		cout << ")\n";

#if 1
		auto bs = std::bind(bigstep<picard>, std::placeholders::_1,
		                    std::cref(F), std::cref(R2), std::cref(M2));
		std::function<FUNCTION<std::vector<REAL>,REAL>(const std::vector<REAL> &)> bsf = bs;
		auto on_domain = from_value<LAZY_BOOLEAN,std::vector<REAL>>(true);
		taylor = lipschitzify(from_algorithm(bsf), from_value<REAL,std::vector<REAL>,REAL>(REAL(1.01)), on_domain, w);
#else
		auto bs = std::bind(bigstep<picard>, std::placeholders::_1,
		                    std::cref(F), std::cref(R2), std::cref(M2),
		                    std::cref(current_t), std::ref(delta_t),
		                    std::ref(taylor), std::ref(max_err));
		std::function<std::vector<REAL>(const std::vector<REAL> &)> bsf = bs;
/*
		w = bigstep<picard>(w, F, in.delta, eps_local, in.R_scale,
		                    in.step_control_alg, current_t, delta_t,
		                    taylor, max_err);
*/
		auto on_domain = from_value<LAZY_BOOLEAN,std::vector<REAL>>(true);
		w = lipschitz_maxnorm(from_algorithm(bsf), REAL(1.01), on_domain, w);
#endif

		while (max_err.mantissa) {
			max_err.mantissa >>= 1;
			max_err.exponent++;
		}

		if (in.ssteps.n_not_delta) {
			for (unsigned j=0; j<in.ssteps.n_smallsteps; j++) {
				DYADIC t = (delta_t / (int)in.ssteps.n_smallsteps) * (int)j;
				std::vector<REAL> y = taylor(REAL(t));
				cout << "taylor( " << (t+current_t) << " ) = ( ";
				cout << y[0];
				for (unsigned k=1; k<y.size(); k++)
					cout << " , " << y[k];
				cout << " ), max_coeff = " << last_taylor_coeff << "\n" << std::flush;
			}
		} else {
			for (; small_t < current_t + delta_t; small_t = small_t + REAL(in.ssteps.small_delta_t).as_DYADIC()) {
				std::vector<REAL> y = taylor(REAL(small_t) - current_t);
				cout << "taylor( " << small_t << " ) = ( ";
				cout << y[0];
				for (unsigned k=1; k<y.size(); k++)
					cout << " , " << y[k];
				cout << " ), max_coeff = " << last_taylor_coeff << "\n" << std::flush;
			}
		}

#if 0
		// w = taylor(delta_t);
		current_t = current_t + delta_t;
#else
		REAL old_t=current_t;
		current_t= current_t + delta_t;
		w = taylor(current_t-old_t);
#endif
		/*
		for (REAL &wj : w)
			wj = approx(wj, -24);*/

		if (!autonomous && !F.is_autonomous()) {
			F = POLYNOMIAL_FLOW(F, delta_t);
			print_iterator(F);
		}

		fprintf(stderr, "bs %d, prec 2^(%d), max_coeff %d   \r",
			bigsteps, max_err.exponent, last_taylor_coeff);
		fflush(stderr);
	}
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
	" <R-scale> <x> <path/to/coeffs> <iv1> ... <ivd> <step-size-alg> | ivp [-iRRAM-OPTS]\n"
"\n"
"  p               precision of output values in decimal places\n"
"  small-steps     # of steps to divide each big-step into\n"
"  delta_t_*       fraction defining the evenly spaced step-size for printing results\n"
"  delta           t-radius for estimation of (R,M)\n"
"  eps             solution-space-radius for estimation of (R,M)\n"
"  R-scale         proceed with each big-step only (<R-scale>)*R in t\n"
"  x               final t to stop the iteration at (or around)\n"
"  path/to/coeffs  coefficients of the flow function in the format '<d> (<nu> <k> <i1> ... <id>)*'\n"
"  iv*             initial values at t=0\n"
"  step-size-alg   0: improved, 1: simple version\n"
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
	cout << "# input: p, { small-steps | delta_t_num/delta_t_den }, delta, eps, R-scale, x, fname, w1, ..., wd, sc\n";
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

	return { p, Smallstep_Control(ssteps.c_str()), delta, eps, final_x, R_scale, F, w, step_control_alg };
}

void compute()
{
	Input in = read_input();

	cout << setRwidth(in.p);

	print_iterator(in.F);

	if (in.F.is_autonomous())
		plot_output<true,!!(METHOD_PICARD-0)>(in);
	else
		plot_output<false,!!(METHOD_PICARD-0)>(in);
}

