#include "iRRAM.h"
#include <vector>
#include <cassert>

namespace iRRAM {
    
  
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
		KR r(K(0));
		K2 xi(K(1));
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

	POLYNOMIAL2 & operator-=(const POLYNOMIAL2 &p)
	{
		for (const I &i : p.c) {
			for (I &j : c) {
				if (i.i == j.i) {
					j.ci -= i.ci;
					goto ok;
				}
			}
			c.push_back({i.i,-i.ci});
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
		POLYNOMIAL2 f(d);
		for (const I &i : p.c) {
			f+=(*this)*i;
		}
		(*this)=f;
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
	auto operator()(const std::vector<K2> &x) const  -> decltype(REAL(0) * x[0])
	{
		using KR = decltype(REAL(0) * x[0]);
		assert(x.size() == d);
		KR r = KR(REAL(0));
		for (const I &i : c) {
			K2 s = K2(REAL(1));
			for (unsigned j=0; j<d; j++)
				s *= power(x[j],int( i.i[j]));
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


} // namespace iRRAM