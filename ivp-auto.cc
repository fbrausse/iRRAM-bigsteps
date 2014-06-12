#include "iRRAM.h"
#include <vector>
#include <list>
#include <cassert>
#include "iRRAM/limit_templates.h"

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
	virtual void clear()
	{
		if (this->release_check())
			return;
		delete this; 
	}

/***********************************************************/
/* local data: coefficients (as a function),  radius and corresponding bound */

	FUNCTION<unsigned int,REAL> _coeff;
	REAL _radius;
	REAL _bound;

public:
	FUNCTIONAL_taylor_sum(
			FUNCTION<unsigned int,REAL> coeff,
			const REAL& radius,
			const REAL& bound
	) {
		_coeff=coeff;
		_radius=radius;
		_bound=bound;
	}

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
		DEBUG0(2,{cerr << "FUNCTIONAL_taylor_sum starting with precision "<<ACTUAL_STACK.actual_prec
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
				DEBUG0(2,{cerr << "FUNCTIONAL_taylor_sum: stop at "
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
		FUNCTION<unsigned int,REAL> coeff,
		const REAL& radius,
		const REAL& bound
) {
	return new FUNCTIONAL_taylor_sum (coeff,radius,bound);
}

//********************************************************************************
// Class for summation of a vector of Taylor series
// 'bound' must be valid for 'radius' in each dimension
//********************************************************************************

class FUNCTIONAL_vector_taylor_sum :public FUNCTIONAL_object<REAL,std::vector<REAL> > {
	virtual void clear()
	{
		if (this->release_check()) return;
		delete this; 
	}

	REAL _radius;
	REAL _bound;
	FUNCTION<unsigned int,std::vector<REAL> > _coeff;

public:
	FUNCTIONAL_vector_taylor_sum(
			FUNCTION<unsigned int,std::vector<REAL> > coeff,
			const REAL& radius,
			const REAL& bound
	) {
		_coeff=coeff;
		_radius=radius;
		_bound=bound;
	}

	std::vector<REAL> eval(const REAL& t)
	{
		std::vector<REAL> result=_coeff(0); // just to get the dimension
		for (unsigned int nu=0;nu<result.size();nu++){
			FUNCTION<unsigned int,REAL>  coeff_nu=projection(_coeff,nu);
			FUNCTION<REAL,REAL> f_nu=taylor_sum(coeff_nu,_radius,_bound);
			result[nu]=f_nu(t);
		}
		return result;
	}
};

/********************************************************/
/* corresponding constructor for FUNCTION		*/
/********************************************************/
inline FUNCTION<REAL,std::vector<REAL> > taylor_sum (
		FUNCTION<unsigned int,std::vector<REAL> > coeff,
		const REAL& radius,
		const REAL& maximum
) {
	return new FUNCTIONAL_vector_taylor_sum(coeff,radius,maximum);
}

//********************************************************************************
//********************************************************************************

class FLOW {
	public:
//local data:
	unsigned int _dimension;
	vector< vector<REAL> > _coeff;
	REAL _coeff_linear_bound;
	REAL _coeff_const_bound;


//information on the state space of the flow
	unsigned int dimension(){return _dimension;};

// default constructor, empty system
	FLOW()
	{
		_dimension=0;
		_coeff_linear_bound=0;
		_coeff_const_bound=0;
	}

// trivial constructor, here only for linear systems
	FLOW(const vector< vector<REAL> > coeff)
	{
		cout << "# This is a prototypical implementation of flows!\n";
		_dimension=coeff.size();
		_coeff=coeff;

		for (unsigned int nu=0;nu<_dimension;nu++)
			for (unsigned int j=1;j<_dimension+1;j++)
				_coeff_linear_bound = maximum(_coeff_linear_bound,abs(coeff[nu][j]));

		for (unsigned int nu=0;nu<_dimension;nu++)
			_coeff_const_bound = maximum(_coeff_const_bound,abs(coeff[nu][0]));
	}

// bound function, simplified algorithm!
	REAL bound_of_flow(const std::vector<REAL> x,const REAL& radius){
		return abs(x)+_coeff_const_bound+radius*_coeff_linear_bound;
	}
	REAL bound_of_solution(const std::vector<REAL> x,const REAL& radius){
		return radius * (abs(x) + _coeff_const_bound)
				+ exp( abs(x)+radius * _coeff_linear_bound);
	}

	REAL coeff(unsigned int i, unsigned int j){
		return _coeff[i][j];
	}
};

class FUNCTIONAL_ivp_solver_simple :public FUNCTIONAL_object<unsigned int,std::vector<REAL> > 
{
public:
/***********************************************************/
/* local data: flow "function" */

	FLOW _flow;

/* The vector taylorpower[nu][i] contains an initial segment of the Taylor coefficients
   for (y_nu)^i.
   As soon as taylor[nu][i][1] is computed, we will also compute 
   taylorpower[nu][j][i] for as many values i (<=j) as needed by the flow.
*/
	std::vector< std::vector< std::vector<REAL> > >taylorpower;

	unsigned int _dimension;

/***********************************************************/
/* trivialer Objekt-Konstruktor: 
   es werden nur die Daten kopiert
*/
	FUNCTIONAL_ivp_solver_simple(
			const std::vector<REAL> & x0,
			const FLOW & flow
	) {
		_flow=flow;
		_dimension=_flow.dimension();

/* kopiere den Startwert x0 in die Tayloreihen als jeweils 0-te Koeffizienten bei Potenz 1*/ 
		taylorpower.resize(_dimension);
		for (unsigned int nu=0;nu<_dimension;nu++){
			taylorpower[nu].resize(2);
			taylorpower[nu][0].resize(1);
			taylorpower[nu][0][0]=REAL(1);
			taylorpower[nu][1].resize(1);
			taylorpower[nu][1][0]=x0[nu];
		}

	}

	/* (3.13) */
	REAL simple_IVP(unsigned int nu,int l)
	{
		REAL sum=REAL(0);
		for (unsigned int j=0;j<_dimension;j++){
			sum += _flow.coeff(nu,j+1)*taylorpower[j][1][l];
		}
		if (l==0) return sum+_flow.coeff(nu,0);
		return sum /(l+1);
	}


/*****************************************************************/
/* Auswertungsfunktion: Taylorreihe in jeder Dimension bestimmen */
/* Werte werden dabei zwischengespeichert                        */
/*****************************************************************/

	std::vector<REAL>  eval(const unsigned int& n)
	{
/* determine Taylor coefficients up to index n (inclusive) */

/* for the simple linear systems we do not need higher powers... */
/* in the general case we would have max_power=n*/
		unsigned int max_power=1; 

		unsigned int l_old=taylorpower[0][0].size();

/* increase the space used to store the  coefficients and their powers*/
		for (unsigned int nu=0;nu<_dimension;nu++){
			taylorpower[nu].resize(max_power+1);
			for (unsigned int l=0;l<=max_power;l++){
				taylorpower[nu][l].resize(n+1);
			}
		}

/* first determine the coefficients for exponent 1, then later the higher powers */
		for (unsigned int l=l_old;l<=n;l++){
			for (unsigned int nu=0;nu<_dimension;nu++){
/* for power 0: we only need to append with sufficiently many zeroes */
				taylorpower[nu][0][l]=REAL(0);
/* for power 1: we have to compute the values explicitly using a recursive scheme */
/*        the case l=0 has already been treated in the constructor */
 				if (l>0) taylorpower[nu][1][l]=simple_IVP(nu,l-1);

/* for powers >1: the values are determined using folding*/
				for (unsigned int m=2;m<=max_power;m++){
					unsigned int n0=0;
					if (m<l)
						n0=l;
					for (unsigned int n1=n0;n1<=l;n1++){
						REAL sum;
						for (unsigned int j=0;j<=n1;j++){
							sum = sum + (taylorpower[nu][1][j] *
							             taylorpower[nu][m-1][n1-j]);
						}
						taylorpower[nu][m][n1]=sum;
					}
				}
			}
		}
		std::vector<REAL> result(_dimension);
		for (unsigned int nu=0;nu<_dimension;nu++){
			result[nu]=taylorpower[nu][1][n];
		}
		return result;
 	}

	virtual void clear()
	{
		if (this->release_check()) return;
		delete this; 
	}
};

/*****************************************************/
/* FUNCTION - constructor for "ivp_solver_simple"    */
/*****************************************************/

inline FUNCTION<unsigned int,std::vector<REAL> > ivp_solver_simple (
	const std::vector<REAL> & x0,
	const FLOW & flow 
) {
	return new FUNCTIONAL_ivp_solver_simple (x0,flow);
}





static inline INTEGER INTEGER_xor(const INTEGER &a, const INTEGER &b)
{
	INTEGER r;
	mpz_xor(r.value, a.value, b.value);
	return r;
}

static inline INTEGER INTEGER_and(const INTEGER &a, const INTEGER &b)
{
	INTEGER r;
	mpz_and(r.value, a.value, b.value);
	return r;
}

/* See http://en.wikipedia.org/wiki/Combinatorial_number_system#Applications */
// find next k-combination
void next_combination(INTEGER& x) // assume x has form x'01^a10^b in binary
{
	// extract rightmost bit 1; u =  0'00^a10^b
	INTEGER u = INTEGER_and(x, -x);
	// set last non-trailing bit 0, and clear to the right; v=x'10^a00^b
	INTEGER v = u +  x;
	// v^x = 0'11^a10^b, (v^x)/u = 0'0^b1^{a+2}, and x ← x'100^b1^a
	x = v + ((INTEGER_xor(v,x) / u) >> 2);
}

/***********************************************************
 * Autonomous flows
 *********************************************************** */

/*
 * To enable doing this:
 * for (Integer_Tuple v({1,n_1,...,n_d}); v; ++v)
 * 	REAL c = AUTONOMOUS_FLOW.c(*v);
*/
template <typename T>
class Integer_Tuple {
protected:
	std::vector<T> v;

public:
	explicit Integer_Tuple(const std::vector<T> &v)
	: v(v) {}

	virtual void set(const std::vector<T> &value)
	{
		v = value;
	}

	virtual void set(unsigned i, const T &value)
	{
		v[i] = value;
	}

	unsigned size() const { return v.size(); }

	const T &              operator[](unsigned i) const { return v[i]; }
	const std::vector<T> * operator->() const { return &v; }
	const std::vector<T> & operator*() const { return v; }
};

class Coeff_Address : public Integer_Tuple<int> {

	bool  invalid;
	const std::vector<int> n;

	void check_validity()
	{
		invalid = (n.size() != v.size());
		for (unsigned i=0; !invalid && i<n.size(); i++)
			if (v[i] > n[i])
				invalid = true;
	}
public:
	Coeff_Address(const vector<int> &value, const vector<int> &ranges)
	: Integer_Tuple(value), n(ranges)
	{
		check_validity();
	}

	explicit Coeff_Address(const vector<int> &ranges)
	: Coeff_Address(std::vector<int>(ranges.size()), ranges)
	{}

	void set(const std::vector<int> &value)
	{
		Integer_Tuple::set(value);
		check_validity();
	}

	void set(unsigned i, const int &value)
	{
		Integer_Tuple::set(i, value);
		check_validity();
	}

	void reset()
	{
		set(std::vector<int>(n.size()));
	}

	void operator++()
	{
		for (unsigned i=0; i<n.size(); i++) {
			if ((v[i] += 1) <= n[i])
				return;
			v[i] = 0;
		}
		invalid = true;
	}

	explicit operator bool() const
	{
		return !invalid;
	}
};

struct auto_f_param {
	unsigned int nu, k;
	Integer_Tuple<int> i;
};

class AUTONOMOUS_FLOW {
	unsigned _dimension;
	FUNCTION<auto_f_param,REAL> c;

public:
	typedef Coeff_Address iterator_type;

	AUTONOMOUS_FLOW(unsigned dimension, const FUNCTION<auto_f_param,REAL> &c)
	: _dimension(dimension), c(c) {}
	
	Coeff_Address iterator(const std::vector<int> &n) const
	{
		assert(n.size() == _dimension);
		return Coeff_Address(n);
	}

	unsigned dimension() const { return _dimension; }

	REAL operator()(unsigned nu, unsigned k, const Coeff_Address &p)
	{
		return c({nu, k, p});
	}
};

template <typename ITuple>
static REAL tan_c(unsigned nu, unsigned k, const ITuple &i)
{
	if (nu != 0) return 0;
	if (k  != 0) return 0;
	if (i.size() != 1) return 0;
	switch (i[0]) {			/* x' = x^0 + x^2 */
	case 0:
	case 2:
		return 1;
	default:
		return 0;
	}
}

REAL tan_c(const auto_f_param &p) { return tan_c(p.nu, p.k, p.i); }

/* Spezielle AUTONOMOUS_FLOW-Implementierung für tan(x):
 * Dimension 1, eigener Iterator. */
class TAN_AUTONOMOUS_FLOW {
public:
	/* Iteriert maximal bis (inkl.) 2, beginnend von 0, in 2er-Schritten.
	 * Feste Dimension 1. */
	class Iterator {
		const int max;
		int i1;
	public:
		Iterator(int max) : max(max), i1(0) {}

		explicit operator bool() const { return i1 <= max; }
		void     operator++()          { i1 += 2; }

		unsigned size()                 const { return 1; }
		int      operator[](unsigned j) const { return i1; }
	};

	/* FUNCTIONAL_ivp_solver_auto<F> nutzt F::iterator_type */
	typedef Iterator iterator_type;

	int mu() const { return 2; }

	/* n[0..(dimension()-1)] = n_1,...,n_{dimension()}, n[dimension()] = k */
	inline iterator_type iterator(const std::vector<int> &n) const
	{
		assert(n.size() == dimension()+1);
		int n1 = n[0];
		int k = n[1];
		return Iterator(k > 0 ? -1 : std::min(n1, mu()));
	}
	inline unsigned dimension() const { return 1; }
	inline REAL operator()(unsigned nu, const iterator_type &idx) const
	{
		return tan_c(nu, 0, idx);
	}

	/* x = tan(0) = 0 */
	REAL bound_of_solution(const std::vector<REAL> &x, const REAL &radius)
	{
		return abs(tan(radius));
	}
};

INTEGER fac(const INTEGER &n)
{
  if (n < 0) return 0;
  INTEGER result(1);
  for (INTEGER j = 1; j <= n; j = j + 1) result = result*j;
  return result;
}

static REAL cos_c(unsigned nu, unsigned k, int i)
{
	if (nu != 0) return 0;
	if (k  != 0) return 0;
	
	if ( (i % 2) == 0) {
	  bool even = ((i/2) % 2) == 0;
	  if (even) return (REAL(1)/(REAL)fac(i));
	  else return (REAL(-1)/(REAL)fac(i));
	}
	
	return 0;
}

// REAL cos_c(const auto_f_param &p) { return cos_c(p.nu, p.k, p.i); }


/* Spezielle AUTONOMOUS_FLOW-Implementierung für die Lösungsfunktion zur Flussfunktion F(x) = cos(x);
 * Lösung zu y'(x) = cos(y(x)) ist y(x) = 2*arctan(tanh(x/2))
 * Dimension 1, eigener Iterator. */
class SOLUTION_AUTONOMOUS_FLOW_COS {
public:
	constexpr int mu() const { return -1; }

#if 0
	class Iterator {
		const int max;
		int i1;
	public:
		Iterator(int max) : max(max), i1(0) {}

		explicit operator bool() const { return i1 <= max; }
		void     operator++()          { i1 += 2; }

		unsigned size()                 const { return 1; }
		int      operator[](unsigned j) const { return i1; }
	};

	typedef Iterator iterator_type;

	/* n[0..(dimension()-1)] = n_1,...,n_{dimension()}, n[dimension()] = k */
	inline Iterator iterator(const std::vector<int> &n) const
	{
		assert(n.size() == dimension()+1);
		int n1 = n[0];
		int k = n[1];
		return Iterator( k != 0 ? -1 : n1);
	}
#else
	class Iterator_reverse_sums {
		const int max;
		int i1;
	public:
		Iterator_reverse_sums(int max) : max(max), i1(0) {}

		explicit operator bool() const { return i1 <= max; }
		void     operator++()          { i1 += 2; }

		unsigned size()                 const { return 1+1; }
		int      operator[](unsigned j) const { return j == 0 ? i1 : 0; }
	};

	typedef Iterator_reverse_sums iterator_type;

	inline Iterator_reverse_sums iterator(unsigned nu, int l) const
	{
		return Iterator_reverse_sums(l);
	}
#endif

	inline unsigned dimension() const { return 1; }
	inline REAL operator()(unsigned nu, const iterator_type &idx) const
	{
		return cos_c(nu, 0, idx[0]);
	}

	/* x = tan(0) = 0 */
	REAL bound_of_solution(const std::vector<REAL> &x, const REAL &radius) const
	{
#if 1
		/* Polstelle bei +/-i*pi/2 */
		COMPLEX pol(0, pi()/2);
		COMPLEX x0(x[0], 0);
		COMPLEX v = pol - x0;
		COMPLEX x_boundary = x0 + (radius / abs(v)) * v;
		return abs(COMPLEX(2)*atan(tanh(x_boundary/COMPLEX(2))));
#else
		return radius * cosh(radius);
#endif
	}
};







/* injectively maps {0,...,dl} to
 * {(i_1,...,i_d,k) where k = 0 and 0 <= i_1,...,i_d <= l}
 *
 * x = t*l+s, s in {1,...,l}, t in {0,...,d-1} */
class Slot_Iterator {

	const unsigned d, l;
	unsigned x;
public:
	Slot_Iterator(unsigned d, unsigned l)
	: d(d), l(l), x(0)
	{}

	unsigned size()          const { return d+1; }
	explicit operator bool() const { return x <= d * l; }
	void     operator++()          { x++; }

	int      operator[](unsigned j) const
	{
		return x == 0 ? 0 : ((x-1)/l == j) ? 1+((x-1)%l) : 0;
	}
};

class TAN_AUTONOMOUS_FLOW_nD {

	const unsigned d;

public:
	TAN_AUTONOMOUS_FLOW_nD(unsigned d = 1) : d(d) {}
#if 0
	typedef Coeff_Address iterator_type;

	iterator_type iterator(std::vector<int> nk) const
	{
		assert(nk.size() == dimension() + 1);
		int k = nk[d];
		for (int &n : nk)
			n = k > 0 ? -1 : min(n, mu());
		nk[d] = 0;
		Coeff_Address ret(nk);
		ret.set(d, k);
		return ret;
	}
#elif 0
	typedef Slot_Iterator iterator_type;
	
	iterator_type iterator(unsigned nu, int l) const
	{
		return Slot_Iterator(d, std::min(l, mu()));
	}
#else
	class Iterator {
		const unsigned d;
		const unsigned nu;
		const unsigned l;
		unsigned v;
	public:
		Iterator(unsigned d, unsigned nu, unsigned l) : d(d), nu(nu), l(l), v(0) {}
		unsigned size() const { return d+1; }
		explicit operator bool() const { return v <= l; }
		void operator++() { v += 2; }
		int operator[](unsigned j) const { return nu == j ? v : 0; }
	};

	typedef Iterator iterator_type;

	iterator_type iterator(unsigned nu, int l) const
	{
		return Iterator(d, nu, std::min(l, mu()));
	}
#endif

	int mu() const { return 2; }

	unsigned dimension() const { return d; }

	REAL operator()(unsigned nu, const iterator_type &it) const
	{
		/*
		fprintf(stdout, "c_{%u,%d", nu, it[d]);
		for (unsigned i=0; i<d; i++)
			fprintf(stdout, ",%d", it[i]);
		fprintf(stdout, "}");
		for (unsigned i=0; i<d; i++)
			fprintf(stdout, ", n_%u = %d", i+1, it.n[i]-1);
		fprintf(stdout, "\n");*/

		int k = it[d];
		int i_nu = it[nu];

		/* this check is not needed anymore in case of Slot_Iterator */
		for (unsigned i=0; i<d; i++)
			if (i != nu && it[i] != 0)
				return 0;

		return tan_c(0, k, std::vector<int>{ i_nu });
	}

	REAL bound_of_solution(const std::vector<REAL> &x, const REAL &radius) const
	{
		return abs(tan(radius));
	}
};

template <typename Flow>
class PARALLEL_FLOW {

	Flow F;
	unsigned d;
public:
	PARALLEL_FLOW(const Flow &F, unsigned d) : F(F), d(d)
	{
		assert(F.dimension() == 1);
	}

	int mu() const { return F.mu(); }

	unsigned dimension() const { return d; }

	class Iterator {
	public:
		typename Flow::iterator_type fit;
		const unsigned d;
		const unsigned nu;
		Iterator(
			const typename Flow::iterator_type &fit,
			unsigned d, unsigned nu
		) : fit(fit), d(d), nu(nu)
		{}

		unsigned size() { return d+1; }
		explicit operator bool() const { return fit ? true : false; }
		void operator++() { ++fit; }
		int operator[](unsigned j) const { return j == nu ? fit[0] : j == d ? fit[1] : 0; }
	};

	typedef Iterator iterator_type;

	iterator_type iterator(unsigned nu, int l) const
	{
		return Iterator(F.iterator(0, l), d, nu);
	}

	REAL operator()(unsigned nu, const iterator_type &it) const
	{
		return it.nu == nu ? F(0, it.fit) : REAL(0);
	}

	REAL bound_of_solution(const std::vector<REAL> &x, const REAL &radius) const
	{
		return F.bound_of_solution(x, radius);
	}
};


/*
 * Riccati-Beispiel
 * y'(x) = 2+2xy(x)+(y(x))^2, y(0) = 0
 * y(x) = (-2)*(e^{-x^2}*E(x))/(1-x*e^{-x^2}*E(x)), E(x) = \int_0^x e^{t^2} dt
 *
 * c_{1,0,0} = 2
 * c_{1,1,1} = 2
 * c_{1,0,2} = 1
 */
class RICCATI_FLOW {
  
public:
  
	/*mu aus Abschnitt 3.5; Falls mu < 0, so sind die Werte der Indizes unbeschränkt*/
	static const int _mu = 2;
  
	int mu()
	{
	  return _mu;
	}
  
	class Iterator {
	public:
		const int k;
		const int n1;
		int i1;
	
		Iterator(int k, int n1) : k(k), n1(n1), i1(0) {}

		explicit operator bool() const  { return i1 <= RICCATI_FLOW::_mu; }
		void     operator++()           { i1++; }
		unsigned dimension() const      { return 1; }
		int      operator[](unsigned j) const { return i1; }
	};

	typedef Iterator iterator_type;

	inline iterator_type iterator(const std::vector<int> &nk) const
	{
		return Iterator(nk[dimension()], nk[0]);
	}

	inline unsigned dimension() const { return 1; }

	inline REAL operator()(unsigned nu, const iterator_type &idx) const
	{
		switch (idx.k) {
		case 0:
			return (idx.i1 == 0) ? 2 : (idx.i1 == 2) ? 1 : 0;
		case 1:
			return (idx.i1 == 1) ? 2 : 0;
		default:
			return 0;
		}
	}

	/* w=0, radius <= 0.9 */
	REAL bound_of_solution(const std::vector<REAL> &w, const REAL &radius) const
	{
		return 41;
	}
};

class SIN_COS {
public:
	static const unsigned _mu = 1;

	inline int mu() const { return _mu; }

	class Iterator {
	public:
		const unsigned nu;
		bool valid;

		Iterator(unsigned nu) : nu(nu), valid(true) {}

		explicit operator bool() const { return valid; }
		void     operator++()          { valid = false; }
		unsigned size()          const { return 4+1; }
		int      operator[](unsigned j) const
		{
			return j == (nu+1) % 4 ? 1 : 0;
		}
	};

	typedef Iterator iterator_type;

	inline iterator_type iterator(unsigned nu, int l) const
	{
		return Iterator(nu);
	}

	inline unsigned dimension() const { return 4; }

	inline REAL operator()(unsigned nu, const iterator_type &idx) const
	{
		return nu == idx.nu ? 1 : 0;
	}

	REAL bound_of_solution(const std::vector<REAL> &w, const REAL &radius) const
	{
		return radius * abs(w) + exp(abs(w)+radius);
	}
};

/*Flow function F(t,x) = t^2 + x^2 -> yields Ricatti-ODE y'(x) = x^2 + (y(x))^2
 * Using the notation of section 4 U_{F}(t_0, x_0, \delta, \varepsilon) := (\delta + abs(t_0))^2 + (\varepsilon + abs(x_0))^2
 * is an upper bound for F on C*/
/* Maxima solution formula: y(x):=(gamma(3/4)*(x^2*bessel_j(-5/4,x^2/2)-x^2*bessel_j(3/4,x^2/2)+bessel_j(-1/4,x^2/2))-x^2*gamma(1/4)*bessel_j(-3/4,x^2/2))/(x*(gamma(1/4)*bessel_j(1/4,x^2/2)-2*gamma(3/4)*bessel_j(-1/4,x^2/2))); */
class RICCATI_FLOW_2
{
public:
	static const unsigned _mu = 2;
	inline int mu() const { return _mu;}

	class Iterator
	{
	public:
		int state;

		Iterator(int l) : state(l >= 2 ? 0 : 2) {}

		explicit operator bool() const { return state < 2; }
		void operator++()
		{
			state++;
		}
		unsigned size()          const { return 2; }
		int operator[](unsigned j) const
		{
			switch(state) {
			case 0: return j == 0 ? 2 : 0;
			case 1: return j == 0 ? 0 : 2;
			}
			return 0;
		} 
	};

	typedef Iterator iterator_type;

	inline iterator_type iterator(unsigned nu, int l) const
	{
		return Iterator(l);
	}

	inline unsigned dimension() const { return 1; }

	inline REAL operator() (unsigned nu, const iterator_type &idx) const
	{
		return 1;
	}

	REAL min(const REAL &a, const REAL &b) { return (a+b-abs(a-b))/2; }

	void get_RM(const std::vector<REAL> &w, const REAL &x0, const REAL &t0, const REAL &delta, const REAL &eps, REAL &R, REAL &M)
	{
		REAL UF = (delta + abs(t0))*(delta + abs(t0)) + (eps + abs(x0))*(eps + abs(x0));
		cout << "U: " << UF << "\n";
		R = min(delta, eps/UF);
		M = abs(w[0]) + R*UF;
	}
};


class EXP {
public:
	static const unsigned _mu = 1;

	inline int mu() const { return _mu; }

	class Iterator {
		bool valid;
	public:
		Iterator() : valid(true) {}
		explicit operator bool() const { return valid; }
		void     operator++()          { valid = false; }
		unsigned size()          const { return 2; }
		int      operator[](unsigned j) const
		{
			return j == 0 ? 1 : 0;
		}
	};

	typedef Iterator iterator_type;

	inline iterator_type iterator(unsigned nu, int l) const
	{
		return Iterator();
	}

	inline unsigned dimension() const { return 1; }

	inline REAL operator()(unsigned nu, const iterator_type &idx) const
	{
		return 1;
	}

	REAL bound_of_solution(const std::vector<REAL> &w, const REAL &radius) const
	{
		return w[0] * exp(radius);
	}
};

template <unsigned _n>
struct I {
	static const unsigned n = _n;
	static_assert(n >= 1, "n = d+1 must be >= 1");
	unsigned ik[n];
	unsigned idx_i_gt0[n];
	unsigned ni_gt0;
};

class VANDERPOL_FLOW
{
public:
	typedef I<3> I_type;

private:
	static const constexpr I_type c[] = {
		{ {0,1,0}, {1  }, 1 }, /* 0+0 */
		{ {0,1,0}, {1  }, 1 }, /* 1+0 */
		{ {2,1,0}, {0,1}, 2 }, /* 1+1 */
		{ {1,0,0}, {0  }, 1 }, /* 1+2 */
	};

	const REAL alpha;

public:
	static const unsigned _mu = 2;

	VANDERPOL_FLOW(const REAL &alpha) : alpha(alpha) {}
	
	inline int mu() const { return _mu; }
	
	class Iterator
	{
	public:
		const unsigned nu;
		int state;

		Iterator(unsigned nu) : nu(nu), state(0) {}

		explicit operator bool() const {
			if (nu == 0) return state == 0 ? true : false;
			if (nu == 1) return state <= 2 ? true : false;
			return false;
		}
		
		void operator++()
		{
			state++;
		}
		
		unsigned size()          const { return 3; }

		const I_type * operator->() const { return &**this; }
		const I_type & operator*() const { return c[nu+state]; }
		
		int operator[](unsigned j) const
		{
#if 1
			return (*this)->ik[j];
#elif 1
			int ret[2][3][3] = { [0]={ {0,0,1} }, [1] = { {0,0,1},{0,2,1},{0,1,0} } };
			return ret[nu][state][j==2 ? 0 : (j+1)];
#endif
		}
	};
	 
	typedef Iterator iterator_type;
	
	inline iterator_type iterator(unsigned nu, int l) const { return Iterator(nu); }
	
	inline unsigned dimension() const { return 2; }
	
	inline REAL operator()(unsigned nu, const iterator_type &idx) const {
		assert(nu == idx.nu);
		
		if (nu == 0) {
			if (idx.state == 0) return 1;
		}
		
		if (nu == 1) {
			switch(idx.state) {
			case 0: return alpha;
			case 1: return 0-alpha;
			case 2: return -1;
			}
		}
		
		return REAL(0);
	}
	
	REAL UF_1(const std::vector<REAL> &w, const REAL &t0, const REAL &delta, const REAL &eps) const
	{
		return abs(w[1]) + eps;
	}
	
	REAL UF_2(const std::vector<REAL> &w, const REAL &t0, const REAL &delta, const REAL &eps) const
	{
		REAL d0 = eps + abs(w[0]);
		REAL d1 = eps + abs(w[1]);
		return alpha * (d1 + square(d0) * d1) + d0;
	}
	
	REAL UF(const std::vector<REAL> &w, const REAL &t0, const REAL &delta, const REAL &eps) const
	{
		return maximum(
			UF_1(w,t0,delta,eps),
			UF_2(w,t0,delta,eps));
	}
	
	void get_RM(const std::vector<REAL> &w, const REAL &t0, const REAL &delta, const REAL &eps, REAL &R, REAL &M) const
	{
		REAL _UF = UF(w,t0,delta,eps);
		// cout << "U: " << _UF << "\n";
		R = minimum(delta, eps/_UF);
		REAL abs_w(0);
		for (unsigned j = 0; j < w.size(); j++) {
			abs_w += square(w[j]);
		}
		abs_w = sqrt(abs_w);
		//M = abs_w + R*_UF;
		// besser: _UF neu berechnen, aber das macht nur bei nicht-autonomen Systemen einen Unterschied
		// NM, 2013-08-31
		M = abs_w + R*UF(w,t0,R,eps);
	}
};

constexpr typename VANDERPOL_FLOW::I_type VANDERPOL_FLOW::c[];


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
		const REAL &eps,
		REAL &R,
		REAL &M
	) const {
		REAL abs_w(0);
		for (unsigned j = 0; j < w.size(); j++) {
			abs_w += square(w[j]);
			w[j] = abs(w[j]);
		}
		REAL t = abs(t0) + delta;
		abs_w = sqrt(abs_w);

		unsigned k = 50;
		REAL rect_w = scale(eps, -(int)k);
		REAL total_w = 0;
		REAL lower_sum = 0;
		for (unsigned j=0; j<k; j++) {
			lower_sum += rect_w / f_max(w, total_w += rect_w, t);// f_max(abs_w + (total_w += rect_w));
			rect_w *= 2;
		}

		R = minimum(delta, lower_sum);
		M = abs_w + eps;
	}

private:
	unsigned _d;
	unsigned _mu;
	bool _autonomous;
	std::vector<std::vector<I_type>> c;
};

template <typename Flow>
REAL flow_lipschitz_C(
	const Flow &F, unsigned nu,
	const std::vector<REAL> &w,
	const REAL &delta, const REAL &eps
) {
	REAL sum(0);
	assert(F.mu() >= 0); /* sonst "eh alles vorbei" :) */

	REAL abs_w(0);
	for (unsigned j = 0; j < w.size(); j++) {
		abs_w += square(w[j]);
	}
	abs_w = sqrt(abs_w);
	REAL abs_wd = abs_w + delta + eps;

	for (auto idx = F.iterator(nu, F.mu()); idx; ++idx) {
		unsigned ik_sum = 0;
		for (unsigned j=0; j<=F.dimension(); j++)
			ik_sum += idx->ik[j];
		if (ik_sum == 0)
			continue;
		sum += abs(F(nu, idx)) * power(abs_wd, ik_sum - 1) * (int)ik_sum;
	}

	return sum;
}

static REAL toREAL(const char *s)
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
	int n, ret;
	FILE *f;
	const char *msg = NULL;

	if (!(f = fopen(fname, "r"))) {
		perror(fname);
		exit(1);
	}
	ret = fscanf(f, "%u", &dimension);			/* d */

	POLYNOMIAL_FLOW F(dimension);
	std::vector<unsigned> ik(dimension+1);
	char *s = NULL;

	if (ret < 0)
		{ msg = "d"; goto err; }

	while (fscanf(f, "%u", &nu) > 0) {			/* nu */
		if (fscanf(f, "%u", &ik[dimension]) < 0)	/* k */
			{ msg = "k"; goto err; }	
		for (unsigned j=0; j<dimension; j++)
			if (fscanf(f, "%u", &ik[j]) < 0)	/* i1,...,id */
				{ msg = "i_j"; goto err; }
		if (fscanf(f, "%*s%n", &n) < 0 || !n)
			{ msg = "real number c_{nu,k,i_1,...,i_d}"; goto err; }
		if (!(s = (char *)realloc(s, n+1))) {
			perror("realloc");
			exit(1);
		}
		fseek(f, -n, SEEK_CUR);
		if (fscanf(f, "%s", s) < 0)			/* c_{nu,k,i1,...,id} */
			{ msg = "real number c_{nu,k,i_1,...,i_d}"; goto err; }

		F.add_coeff(nu, VI(toREAL(s), ik));
	}

	free(s);
	fclose(f);

	return F;
err:
	fprintf(stderr,
		"invalid coefficient input file: expected %s at offset %lu\n",
		msg, ftell(f));
	free(s);
	fclose(f);
	exit(1);
}



#if 0
template <unsigned B>
class Combination {

	const unsigned n, k;
	unsigned long c[B];

	unsigned idx(unsigned bit)
	{
		return bit / (8*sizeof(*c));
	}

	unsigned mask1(unsigned bit)
	{
		return 1UL << (bit % (8*sizeof(*c)));
	}

public:
	Combination(unsigned n, unsigned k)
	: n(k), k(k)
	{
		unsigned i = idx(k);
		memset(c, 0, sizeof(c));
		c[i] = mask1(k) - 1UL;
		while (i--)
			c[i] = ~0UL;
	}

	explicit operator bool() const
	{
		return (c[idx(k)] & mask1(k)) ? true : false;
	}

/* See http://en.wikipedia.org/wiki/Combinatorial_number_system#Applications */
// find next k-combination
void next_combination(INTEGER& x) // assume x has form x'01^a10^b in binary
{
	// extract rightmost bit 1; u =  0'00^a10^b
	INTEGER u = INTEGER_and(x, -x);
	// set last non-trailing bit 0, and clear to the right; v=x'10^a00^b
	INTEGER v = u +  x;
	// v^x = 0'11^a10^b, (v^x)/u = 0'0^b1^{a+2}, and x ← x'100^b1^a
	x = v + ((INTEGER_xor(v,x) / u) >> 2);
}

	void operator++()
	{
		unsigned i;
		unsigned u_idx = 0;
		while (!c[u_idx])
			u_idx++;
		unsigned long u_bit = c[u_idx] & -c[u_idx];
		
	}

	template <typename R>
	void get(R &r) const
	{
		unsigned long bit_idx0 = -1;
		unsigned long bit_idx1;
		for (unsigned i=1; i<n; i++, bit_idx0 = bit_idx1) {
			bit_idx1 = mpz_scan0(combination.value, bit_idx0 + 1);
			r[i-1]   = bit_idx1 - bit_idx0 - 1;
		}
		r[n-1] = (n+k-1) - bit_idx0 - 1;
	}
};
#endif

class Combination {

	const unsigned n, k;
	INTEGER combination;

public:
	Combination(unsigned n, unsigned k)
	: n(n), k(k), combination((INTEGER(1) << (int)k) - 1)
	{}

	explicit operator bool() const
	{
		return mpz_tstbit(combination.value, n+k-1) == 0;
	}

	void operator++() { next_combination(combination); }

	const INTEGER & operator*() const { return combination; }

	template <typename R>
	void get(R &r) const
	{
#if 0
		INTEGER bit(1);
		for (unsigned i=0; i<n; i++) {
			unsigned k;
			for (k=0; INTEGER_and(combination, bit) != 0; k++)
				bit *= 2;
			bit *= 2;
			r[i] = k;
		}
#else
		unsigned long bit_idx0 = -1;
		unsigned long bit_idx1;
		for (unsigned i=1; i<n; i++, bit_idx0 = bit_idx1) {
			bit_idx1 = mpz_scan0(combination.value, bit_idx0 + 1);
			r[i-1]   = bit_idx1 - bit_idx0 - 1;
		}
		r[n-1] = (n+k-1) - bit_idx0 - 1;
#endif
	}
};

#include <sys/time.h>

struct Timer {
	unsigned long t = 0;

	void start() { struct timeval tv; gettimeofday(&tv, NULL); t -= tv.tv_sec * 1000000UL + tv.tv_usec; }
	void stop() { struct timeval tv; gettimeofday(&tv, NULL); t += tv.tv_sec * 1000000UL + tv.tv_usec; }
} t_auto, t_rest0, t_rest1;

#define dbg(fmt, ...) DEBUG0(2,{fprintf(stderr, (fmt), ##__VA_ARGS__);})

template<typename F>
class FLOW_PROXY
{
	const std::vector<int> v;
	const F flow;
public:
	class Iterator {
	public:
		typename F::iterator_type fit;
		const std::vector<int>    &v;

		bool v_contains(int i) const
		{
			/* TODO: optimize to binary search */
			for (unsigned j=0; j<v.size(); j++)
				if (v[j] == i)
					return true;
			return false;
		}

		void search_fit_next()
		{
			for (; fit; ++fit) {
				unsigned j;
				#if 0
				unsigned non_zeros = 0;
				for (j=0; j<fit.size()-1; j++)
					if (fit[j] != 0)
						non_zeros++;
				if (non_zeros <= 1)
					break;
				#else
				for (j=0; j<fit.size()-1; j++)
					if (!v_contains(j) && fit[j] != 0) {
						// printf("discarding i_%d: fit[%d] = %d != 0\n", j, j, fit[j]);
						break;
					}
				if (!(j < fit.size()-1))
					break;
				#endif
			}
		}

		Iterator(typename F::iterator_type fit, const std::vector<int> &v)
		: fit(fit), v(v) { search_fit_next(); }
		explicit operator bool() const { return fit ? true : false; }
		void     operator++()          { ++fit; search_fit_next(); }
		unsigned size()          const { return v.size() + 1; }
		int      operator[](unsigned j) const
		{
			return fit[j == v.size() ? fit.size()-1 : v[j]];
		}
	};

	typedef Iterator iterator_type;

	FLOW_PROXY(const std::vector<int> &v, const F &flow) : v(v), flow(flow)
	{
		assert(flow.dimension() >= v.size());
	}

	iterator_type iterator(unsigned nu, int l) const
	{
		return Iterator(flow.iterator(v[nu], l), v);
	}

	int mu() const
	{
		return flow.mu(); /* TODO: may be too large */
	}

	unsigned dimension() const {
		return v.size();
	}

	REAL operator()(unsigned nu, const iterator_type &idx) const {
		return flow(v[nu], idx.fit);
	}
};

template <typename F>
class FUNCTIONAL_ivp_solver_auto :public FUNCTIONAL_object<unsigned int,std::vector<REAL> >
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
	{/*
		if (l == 0) {
			for (unsigned k=j; k < n.size(); k++) n[k] = 0;
			sum += a_vector_power(n, idx);
			return;
		}*/

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
			unsigned zero_i = 0;
			for (unsigned i=0; i<_dimension; i++)
				if (idx[i] == 0)
					zero_i++;

			for (Combination comb(_dimension - zero_i, l_); comb; ++comb) {
				comb.get(n);

				if (iv_is_zero)
					for (unsigned i=0, j=0; i<_dimension; i++)
						if (idx[i] > 0)
							n[j++] += idx[i];

				dbg("\t");
				for (unsigned i=0, j=0; i<_dimension; i++)
					dbg("n_%u = %d, ", i+1, idx[i] == 0 ? 0 : n[j++]);
				dbg("\n");

#if 1
				REAL mul(1);
				for (unsigned i=0, j=0; i<_dimension; i++)
					if (idx[i] > 0)
						mul *= a[i][idx[i]][n[j++]];
				sum_a_ni += mul;
#else
				for (unsigned i=0, j=0; i<_dimension; i++)
					if (idx[i] == 0)
						n_[i] = 0;
					else
						n_[i] = n[j++];

				sum_a_ni += a_vector_power(n_, idx);
#endif
			}
#elif 0
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
			DEBUG0(3,{cerr << "a(" << nu << ",1," << n << ") = " << result[nu] << "\n";});
		}
		return result;
	}

	virtual void clear()
	{
		if (this->release_check()) return;
		delete this; 
	}
};

struct F_REAL {
	REAL x;
	bool valid;
	bool is_zero;

	inline F_REAL() : valid(false), is_zero(false) {}
	inline explicit F_REAL(int v) : x(v), valid(true), is_zero(v == 0) {}/*
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
	}*/
};

template <typename F>
class FUNCTIONAL_IVP_SOLVER_RECURSIVE : public FUNCTIONAL_object<unsigned int,std::vector<REAL> > {
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
	inline void combination_recursive3(const REAL &mul, unsigned j, unsigned d, REAL &sum, const V &i, const W &idx, unsigned l)
	{
		if (j == d-1) {
			unsigned k = l;
			sum += mul * a(i[j],idx[i[j]],(iv_is_zero ? idx[i[j]] : 0) + k);/*
		} else if (l == 0) {
			unsigned k = 0; 
			combination_recursive3(mul*a(i[j],idx[i[j]],(iv_is_zero ? idx[i[j]] : 0) + k), j+1, d, sum, i, idx, l-k);*/
		} else {
			for (unsigned k = 0; k <= l; k++)
				combination_recursive3(mul*a(i[j],idx[i[j]],(iv_is_zero ? idx[i[j]] : 0) + k), j+1, d, sum, i, idx, l-k);
		}
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
/*
		std::vector<int> n;
		std::vector<int> n_;
		n.resize(_dimension);
		n_.resize(_dimension);
*/
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
#if 0
			if (l_ == 0) {
				REAL mul(1);
				for (unsigned j=0; j<_dimension; j++)
					mul *= a(j,idx[j],iv_is_zero ? idx[j] : 0);
				sum += REAL(c) * mul;
				continue;
			}
#if 0
			REAL sum_a_ni = 0;
			unsigned zero_i = 0;
			for (unsigned i=0; i<_dimension; i++)
				if (idx[i] == 0)
					zero_i++;

			for (Combination comb(_dimension - zero_i, l_); comb; ++comb) {
				comb.get(n);

				if (iv_is_zero)
					for (unsigned i=0, j=0; i<_dimension; i++)
						if (idx[i] > 0)
							n[j++] += idx[i];

				dbg("\t");
				for (unsigned i=0, j=0; i<_dimension; i++)
					dbg("n_%u = %d, ", i+1, idx[i] == 0 ? 0 : n[j++]);
				dbg("\n");


				REAL mul(1);
				for (unsigned i=0, j=0; i<_dimension; i++)
					if (idx[i] > 0)
						mul *= a(i,idx[i],n[j++]);
				sum_a_ni += mul;
			}
			sum += c * sum_a_ni;
#elif 0
			REAL sum_a_ni = 0;
			unsigned j=0;
			for (unsigned i=0; i<_dimension; i++)
				if (idx[i] > 0)
					n[j++] = i;
			if (j > 0) {
				combination_recursive3(REAL(1), 0, j, sum_a_ni, n, idx, l_);
				sum += c * sum_a_ni;
			}
#elif 0
			const typename F::I_type &I = *idx;
			if (I.ni_gt0 > 0) {
				REAL sum_a_ni(0);
				combination_recursive3(REAL(1), 0, I.ni_gt0, sum_a_ni, I.idx_i_gt0, I.ik, l_);
				sum += c * sum_a_ni;
			}
#endif
#else
			const typename F::I_type &I = *idx;
			if (I.ni_gt0 > 0) {
				sum += c * combination_recursive4(0, I.ni_gt0-1, I.idx_i_gt0, I.ik, l_);
			} else if (l_ == 0) {
				sum += c;
			} else {
				/* some a_{\nu,n}^{(i)} is 0 whereas n >= 1:
				 * nothing to add */
			}
#endif
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
			DEBUG0(3,{cerr << "a(" << nu << ",1," << n << ") = " << result[nu] << "\n";});
		}

		//t.stop();
		//fprintf(stderr, "n: %u,\tt: %luµs   \r", n, t.t);
		//fflush(stderr);

		return result;
	}

	virtual void clear()
	{
		if (this->release_check()) return;
		delete this; 
	}
};

class FUNCTIONAL_IVP_SOLVER_PICARD : public FUNCTIONAL_object<unsigned int,std::vector<REAL>> {
public:
	const POLYNOMIAL_FLOW F;
	/* p[n][nu] */
	std::vector<std::vector<POLYNOMIAL<REAL>>> p;
	unsigned n = 0;

	const std::vector<REAL> w;
	const bool iv_is_zero;

	// const unsigned int _dimension;

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

	virtual void clear()
	{
		if (this->release_check()) return;
		delete this; 
	}
};

template <typename F>
inline FUNCTION<unsigned int, std::vector<REAL> > ivp_solver_recursive(
	const F &flow,
	const std::vector<REAL> &w,
	bool iv_is_zero
) {
	return new FUNCTIONAL_IVP_SOLVER_RECURSIVE<F>(flow, w, iv_is_zero);
}

inline FUNCTION<unsigned int, std::vector<REAL> > ivp_solver_picard(
	const POLYNOMIAL_FLOW &flow,
	const std::vector<REAL> &w,
	bool iv_is_zero
) {
	return new FUNCTIONAL_IVP_SOLVER_PICARD(flow, w, iv_is_zero);
}

template <typename F>
inline FUNCTION<unsigned int, std::vector<REAL> > ivp_solver_auto(
	const F &flow,
	const std::vector<REAL> &w,
	bool iv_is_zero
) {
	return new FUNCTIONAL_ivp_solver_auto<F>(flow, w, iv_is_zero);
}


template <typename T>
POLYNOMIAL<T> poly(const T &t)
{
	return { t };
}

void test()
{
	auto p1 = poly(REAL(0));
	auto p2 = poly(p1);
	auto p3 = poly(p2);
	auto p4 = poly(p3);

	//unsigned d = 5;
	// POLYNOMIAL_FLOW F(d);
	//std::vector<POLYNOMIAL2> p = {
	//	POLYNOMIAL2(d), POLYNOMIAL2(d), POLYNOMIAL2(d), POLYNOMIAL2(d), POLYNOMIAL2(d) };
	// F(p);
}

template <typename KP,typename Flow>
POLYNOMIAL<REAL> picard_iteration_step(
	const REAL &t0, const REAL &x0,
	const Flow &F, const POLYNOMIAL<KP> &phi
) {
	const auto integral = F(phi).primitive();
	return integral + POLYNOMIAL<REAL>(x0 - integral(t0));
}

static unsigned last_phi_n;

template <typename Flow>
class SOLVER_PICARD : public FUNCTIONAL_object<REAL,REAL> {

	const Flow F;
	const REAL t0, x0;
	const REAL L, M;

	unsigned n;
	POLYNOMIAL<REAL> phi;

public:
	SOLVER_PICARD(
		const Flow &F,
		const REAL &t0, const REAL &x0,
		const REAL &delta, const REAL &eps
	) : F(F), t0(t0), x0(x0),
	    L(abs(F.derivative())(abs(x0) + eps)),
	    M(abs(F)(abs(x0) + eps)),
	    n(0), phi(x0)
	{}

/****************************************************************************************/
/* Evaluation:										*/
/* Compute coefficients until the truncation error term is smaller than actual_prec	*/
/* or until the accumulated error is already of the order of the truncation error	*/
/* (whichever comes first)								*/
/****************************************************************************************/

	REAL eval(const REAL &t)
	{
		ITERATION_STACK SAVED_STACK;
		ACTUAL_STACK.inlimit+=1;
		REAL sum=0;
		REAL error = M/L;
		REAL errorfactor = L*abs(t-t0);
		DEBUG0(2,{cerr << "TAN_PL starting with precision "<<ACTUAL_STACK.actual_prec
			<< " at ratio "<< errorfactor.vsize.mantissa*pow(2,errorfactor.vsize.exponent)<<"/(n+1)\n";});

		REAL best=0;
		sizetype sum_error,trunc_error,best_error;

		sizetype_add(best_error,error.vsize,error.error);
		int best_index=-1;

		for (unsigned int i=0;;i++){
			sum = phi(t);
			error= error*errorfactor / (int)(i+1);

			sum.geterror(sum_error);
			sizetype_add(trunc_error,error.vsize,error.error);

			if (i == n) {
				cerr << "phi_" << i;
				// cerr << " = " << phi;
				cerr << ", trunc_err: " << trunc_error.mantissa
				     << "*2^(" << trunc_error.exponent
				     << ")   \n";
			}

			sizetype local_error;
			sizetype_add(local_error,trunc_error,sum_error);
			if (sizetype_less(local_error,best_error)) {
				best=sum;
				best_error=local_error;
				best.seterror(best_error);
				best_index=i;
			}
#if 0 /* TODO: ??? */
			while ((trunc_error.mantissa >>= 1) != 0)
				trunc_error.exponent++;
#endif
			if (trunc_error.exponent < ACTUAL_STACK.actual_prec ||
			    sizetype_less(trunc_error,sum_error)) {
				DEBUG0(2,{cerr << "TAN_PL: stop at " <<i<< " with best index at " << best_index << "\n";});
				last_phi_n = i;
				break;
			}
			if (i == n) {
				phi = picard_iteration_step(t0, x0, F, phi);
				n = i+1;
			}
		}
		return best;
	}

	virtual void clear()
	{
		if (this->release_check()) return;
		delete this; 
	}
};

template <typename Flow>
inline FUNCTION<REAL,REAL> solver_picard(
	const Flow &F,
	const REAL &t0,
	const REAL &x0,
	const REAL &delta,
	const REAL &eps
) {
	return new SOLVER_PICARD<Flow>(F, t0, x0, delta, eps);
}


struct sol_desc {
	REAL (*solution)(const REAL &x);
	const char *name;
};
#define sol_desc0(fun)	{ [](const REAL &x){ return fun; }, #fun }

template <typename F>
std::vector<sol_desc> solution() { return { { nullptr, "<unknown>" } }; }

template <> std::vector<sol_desc> solution<TAN_AUTONOMOUS_FLOW_nD>()
	{ return { sol_desc0(tan(x)) }; }
template <> std::vector<sol_desc> solution<TAN_AUTONOMOUS_FLOW>()
	{ return { sol_desc0(tan(x)) }; }
template <> std::vector<sol_desc> solution<SOLUTION_AUTONOMOUS_FLOW_COS>()
	{ return { sol_desc0(2*atan(tanh(x/2))) }; }
template <> std::vector<sol_desc> solution<SIN_COS>()
	{ return { sol_desc0(sin(x)), sol_desc0(cos(x)) }; }
template <> std::vector<sol_desc> solution<EXP>()
	{ return { sol_desc0(exp(x)) }; }

struct Input {
	unsigned p, d, n;
	REAL radius, x;
	std::vector<REAL> w;
	bool w_is_zero;

	void set_w(const std::vector<REAL> &_w)
	{
		w = _w;
		w_is_zero = false;
	}

	void clear_w()
	{
		w = std::vector<REAL>(d);
		w_is_zero = true;
	}
};

Input read_input()
{
	Input in;
	cout << "input: p, d, radius, x, n\n";
	cin >> in.p;
	cin >> in.d;
	cin >> in.radius;
	cin >> in.x;
	cin >> in.n;
	in.clear_w();
	return in;
}

template <typename Flow>
void write_output(const Input &in, const Flow &F, const REAL &R, const REAL &M)
{
	FUNCTION<unsigned int,std::vector<REAL>> a = ivp_solver_auto(F, in.w, in.w_is_zero);
	FUNCTION<REAL,std::vector<REAL>> taylor;

	cout << "M = " << M << "\n";
	/*
	int n = (REAL(-(p-7))/(log(REAL(2))/log(REAL(10))) * log(REAL(2))+log(R-abs(x))-log(M*R))/(log(abs(x))-log(R))-1;
	cout << "n            = " << n << "\n";*/
	taylor = taylor_sum(a, R, M);

	int plot=100;
	for (int m=0; m<= plot; m++){
	  REAL z=REAL(in.x)*m/plot;
	  std::vector<REAL> y = taylor(z);
		cout << "taylor( " << z << " ) = ( ";
  	  cout << y[0];
	  for (unsigned nu=1; nu<y.size(); nu++)
		cout << " , " << y[nu];
	  cout << " )\n";
	}
	for (unsigned int j = 0; j < in.n; j++) {
		std::vector<REAL> aj = a(j);
		for (unsigned nu=0; nu<aj.size(); nu++)
			cout << "a" << j << "[" << nu << "] = " << aj[nu] << std::endl;
	}
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

template <bool autonomous,bool picard,typename Flow>
void plot_output(
	std::vector<REAL> w,
	const REAL &final_t,
	Flow F,
	const REAL &delta,
	const REAL &eps,
	const REAL &R_scale,
	const Smallstep_Control &smallsteps
) {
	FUNCTION<unsigned int,std::vector<REAL>> a;
	FUNCTION<REAL,std::vector<REAL>> taylor;

	const int cmp_p = -10;
	const int delta_t_p = -53;
	const DYADIC end_t = approx(final_t, cmp_p) + scale(INTEGER(2), cmp_p);

	DYADIC_precision dyadic_prec(delta_t_p);

	DYADIC current_t = 0;
	DYADIC small_t = 0;

	Timer t;
	t.start();

	unsigned bigsteps;
	// for (bigsteps = 0; !positive(REAL(current_t) - final_t, cmp_p); bigsteps++) {
	for (bigsteps = 0; current_t <= end_t; bigsteps++) {
		REAL R, R_test, M, M_test;
		REAL R2, M2;
		DYADIC delta_t;
		sizetype err, max_err;

		t.stop();
		cout << "# --------------------------------------------\n";
		cout << "# big step " << bigsteps
		     << ", last max coeff: " << last_taylor_coeff
		     << ", " << (t.t / 1e3) << "ms\n";
		t.t = 0;
		t.start();

		F.get_RM(w, 0, delta, eps, R, M);
		F.get_RM2(w, 0, delta, eps, R2, M2);
		cout << "# (R ,M ) = (" << R << ", " << M << ")\n";
		cout << "# (R2,M2) = (" << R2 << ", " << M2 << ")\n";

		/* TODO: while (INTEGER(R2 * R_scale * 2**n)-1 <= 0) n++;
		 *  und  DYADIC_precision(n) */
		delta_t = approx(R2 * R_scale, delta_t_p) - scale(DYADIC(1), delta_t_p);
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
		cout << ")\n";/*
		cout << "# L = (" << flow_lipschitz_C(F, 0, w, REAL(delta_t), eps);
		for (unsigned k=1; k<F.dimension(); k++)
			cout << ", " << flow_lipschitz_C(F, k, w, REAL(delta_t), eps);
		cout << ")\n";*/
		/*
		cout << "# no. bits lost: (";
		for (unsigned k=0; k<F.dimension(); k++)
			cout << (k == 0 ? "" : ", ")
			     << (flow_lipschitz_C(F, k, w, REAL(delta_t), eps) * REAL(delta_t) * REAL(1)/log(REAL(2)) +
			         err.exponent * log((w[k].geterror(err), 
			             scale(REAL((double)err.mantissa)))));
		cout << ")\n";*/

		a = picard ? ivp_solver_picard(F, w, false) : ivp_solver_recursive(F, w, false);
		taylor = taylor_sum(a, R2, M2);

		while (max_err.mantissa) {
			max_err.mantissa >>= 1;
			max_err.exponent++;
		}

		if (smallsteps.n_not_delta) {
			for (unsigned j=0; j<smallsteps.n_smallsteps; j++) {
				DYADIC t = (delta_t / (int)smallsteps.n_smallsteps) * (int)j;
				std::vector<REAL> y = taylor(REAL(t));
				cout << "taylor( " << (t+current_t) << " ) = ( ";
				cout << y[0];
				for (unsigned k=1; k<y.size(); k++)
					cout << " , " << y[k];
				cout << " ), max_coeff = " << last_taylor_coeff << "\n" << std::flush;
			}
		} else {
			for (; small_t < current_t + delta_t; small_t = small_t + DYADIC(smallsteps.small_delta_t)) {
				std::vector<REAL> y = taylor(REAL(small_t - current_t));
				cout << "taylor( " << small_t << " ) = ( ";
				cout << y[0];
				for (unsigned k=1; k<y.size(); k++)
					cout << " , " << y[k];
				cout << " ), max_coeff = " << last_taylor_coeff << "\n" << std::flush;
			}
		}

		w = taylor(delta_t);
		for (REAL &wj : w)
			wj = approx(wj, -24);
		current_t = current_t + delta_t;

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


#if 0
struct tan_flow {
	template <typename T>
	T operator()(const T &p) const { return p*p + T(1); }
	POLYNOMIAL<int> derivative() const { return POLYNOMIAL<int>({0,2}); }

	friend POLYNOMIAL<int> abs(const tan_flow &c) { return abs(POLYNOMIAL<int>({1,0,1})); }
};

void compute()
{
	int p;
	std::string ssteps;
	REAL delta, eps, final_t;

	cout << "# input: p, { small-steps | delta_t_nom/delta_t_den }, delta, eps, x\n";

	cin >> p;
	cin >> ssteps;
	cin >> delta;
	cin >> eps;
	cin >> final_t;

	cout << setRwidth(p);

	const int delta_t_p = -53;

	//POLYNOMIAL<int> F({1,0,1});
	tan_flow F;
	FUNCTION<REAL,REAL> phi = solver_picard(F, 0, 0, delta, eps);
	Smallstep_Control smallsteps(ssteps.c_str());
	DYADIC delta_t = approx(final_t, delta_t_p) - scale(DYADIC(1), delta_t_p);
	DYADIC current_t = 0;
	DYADIC small_t = 0;

	if (smallsteps.n_not_delta) {
		for (unsigned j=1; j<=smallsteps.n_smallsteps; j++) {
			DYADIC t = (delta_t / (int)smallsteps.n_smallsteps) * (int)j;
			REAL y = phi(REAL(t));
			cout << "taylor( " << (t+current_t)
			     << " ) = ( " << y
			     << " ), last_phi_n = " << last_phi_n << "\n";
		}
	} else {
		for (; small_t < current_t + delta_t; small_t = small_t + DYADIC(smallsteps.small_delta_t)) {
			REAL y = phi(REAL(small_t - current_t));
			cout << "taylor( " << small_t << " ) = ( "
			     << " ) = ( " << y
			     << " ), last_phi_n = " << last_phi_n << "\n";
		}
	}
}
#elif 1
#include <cstdarg>
static void die(int exitval, const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	exit(exitval);
}

static const char usage[] =
"echo <p> { <small-steps> | <delta_t_num>/<delta_t_den> } <delta> <eps>"
	" <R-scale> <x> <path/to/coeffs> <iv1> ... <ivd> | ivp [-iRRAM-OPTS]\n"
"\n"
"  p            precision of output values in decimal places\n"
"  small-steps  # of steps to divide each big-step into\n"
"  delta_t_*    fraction defining the evenly spaced step-size for printing results\n"
"  delta        t-radius for estimation of (R,M)\n"
"  eps          solution-space-radius for estimation of (R,M)\n"
"  R-scale      proceed with each big-step only (<R-scale>)*R in t\n"
"  x            final t to stop the iteration at (or around)\n"
"  path/to/coeffs  coefficients of the flow function in the format '<d> (<nu> <k> <i1> ... <id>)*'\n"
"  iv*          initial values at t=0\n"
;

static void input_error(const char *which)
{
	die(1, "input error: %s invalid\n\n%s", which, usage);
}

void compute()
{
	int p;
	std::string ssteps;
	REAL delta, eps, final_x, alpha, R_scale;
	std::string fname;

	cout << "# input: p, { small-steps | delta_t_num/delta_t_den }, delta, eps, R-scale, x, fname, w1, ..., wd\n";
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
		if (!(cin >> s) || !s.length()) input_error("<iv%d>");
		w[i] = toREAL(s.c_str());
	}

	cout << setRwidth(p);

	print_iterator(F);

#ifndef METHOD_PICARD
# define METHOD_PICARD	0
#endif

	if (F.is_autonomous())
		plot_output<true,!!(METHOD_PICARD-0)>(w, final_x, F, delta, eps, R_scale, Smallstep_Control(ssteps.c_str()));
	else
		plot_output<false,!!(METHOD_PICARD-0)>(w, final_x, F, delta, eps, R_scale, Smallstep_Control(ssteps.c_str()));
}
#elif 1
void compute()
{
	int p;
	int ssteps;
	REAL delta, eps, final_x, alpha, R_scale;
	std::vector<REAL> w(2);

	cout << "# input: p, small-steps, delta, eps, R-scale, x\n";
	cout << "# van-der-pol input: y(0) y'(0) alpha\n";

	cin >> p;
	cin >> ssteps;
	cin >> delta;
	cin >> eps;
	cin >> R_scale;
	cin >> final_x;

	cin >> w[0];
	cin >> w[1];
	cin >> alpha;

	cout << setRwidth(p);

	VANDERPOL_FLOW F(alpha);
	// w = {1,0};

	print_iterator(F);

	plot_output(w, final_x, F, delta, eps, R_scale, ssteps);
}

#else
void compute()
{
	Input in = read_input();
	REAL x = in.x;

	REAL R = in.radius;
#if 0
	EXP F;
	in.set_w({1});
	REAL M = F.bound_of_solution(in.w, R);
#elif 0
	SOLUTION_AUTONOMOUS_FLOW_COS F;
	// PARALLEL_FLOW<SOLUTION_AUTONOMOUS_FLOW_COS> F(F_, in.d);
	REAL M = F.bound_of_solution(in.w, R);
	// RICCATI_FLOW F;
#elif 0
	TAN_AUTONOMOUS_FLOW_nD F(in.d);
	//std::vector<int> v = { 42 };
	//FLOW_PROXY<TAN_AUTONOMOUS_FLOW_nD> F(v, F_);
	REAL M = F.bound_of_solution(in.w, R);
#elif 0
	RICCATI_FLOW_2 F;
	REAL M;
	in.set_w({1});
	F.get_RM(in.w, in.x, 0, 1, 1, R, M);
	cout << "R: " << R << "\n";
#elif 1
	VANDERPOL_FLOW F(0.5);
	REAL M;
	in.set_w({1,0});
	F.get_RM(in.w, in.x, 0, 100, 0.75, R, M);
	cout << "R: " << R << "\n";
#else
	SIN_COS F;
	in.set_w({0, 1, 0, -1});
	REAL M = F.bound_of_solution(in.w, R);/*
	write_output(in, FLOW_PROXY<SIN_COS>(std::vector<int>{0}, F), R, M);
	write_output(in, FLOW_PROXY<SIN_COS>(std::vector<int>{1}, F), R, M);
	write_output(in, FLOW_PROXY<SIN_COS>(std::vector<int>{2}, F), R, M);
	write_output(in, FLOW_PROXY<SIN_COS>(std::vector<int>{3}, F), R, M);*/
#endif

#if 1
	print_iterator(F);
#endif

	cout << setRwidth(in.p);
#if 1
	for (const sol_desc &fd : solution<decltype(F)>()) {
		cout << setiosflags(std::ios::left) << setw(17)
		     << fd.name << " = ";
		if (fd.solution != nullptr)
			cout << (*fd.solution)(x);
		else
			cout << "<unknown>";
		cout << "\n";
	}
#else
	// cout << "y(x)         = " << y(x) << "\n";
	cout << "tan(x)       = " << tan(x) << "\n";
	cout << "sin(x)       = " << sin(x) << "\n";
	cout << "cos(x)       = " << cos(x) << "\n";
	cout << "exp(x)       = " << exp(x) << "\n";
#endif
	write_output(in, F, R, M);
}
#endif
