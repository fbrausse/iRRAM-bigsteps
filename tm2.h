
#ifndef TAYLOR_MODEL_H
#define TAYLOR_MODEL_H

#include <map>
#include <array>
#include <tuple>
#include <memory>

#include <iRRAM.h>

namespace iRRAM {

template <unsigned N_MAX,typename sweep_mul>
class TM2 {
public:
	typedef REAL c0_t, cn_t;
	typedef unsigned lambda_refcnt_t;
	typedef unsigned multiplicity_T;
	typedef std::shared_ptr<lambda_refcnt_t>                      lambda_t;
	typedef std::array<std::tuple<lambda_t,multiplicity_T>,N_MAX> idx_t;

	friend sweep_mul;

	// friend idx_t operator+(const idx_t &a, const idx_t &b) { return a; }

	TM2(c0_t &&r) : c0(r) {}
	TM2(const c0_t &r = c0_t()) : c0(r) {}
	TM2(TM2 &&t) = default;
	TM2(const TM2 &t) = default;

	explicit operator REAL() const;

	static void polish(TM2 &t...);

	TM2 & operator=(const TM2 &a) = default;

	TM2 & operator+=(const TM2 &b)
	{
		c0 += b.c0;
		for (auto bi = b.c.cbegin(); bi != b.c.end(); ++bi)
			c[bi->first] += bi->second;
		return *this;
	}
	TM2 & operator-=(const TM2 &b)
	{
		c0 -= b.c0;
		for (auto bi = b.c.cbegin(); bi != b.c.end(); ++bi)
			c[bi->first] -= bi->second;
		return *this;
	}
	TM2 & operator*=(const TM2 &b) { return *this = (*this) * b; }

	friend TM2 operator-(const TM2 &a)
	{
		TM2 r = a.c0;
		for (auto ai = a.c.cbegin(); ai != a.c.end(); ++ai)
			r.c.emplace(ai->first, -ai->second);
		return r;
	}
	friend TM2 operator+(const TM2 &a, const TM2 &b) { return TM2(a) += b; }
	friend TM2 operator-(const TM2 &a, const TM2 &b) { return TM2(a) -= b; }
	friend TM2 operator*(const TM2 &a, const TM2 &b)
	{
		TM2 r = a.c0 * b.c0;

		for (auto ai = a.c.cbegin(); ai != a.c.end(); ++ai) {
			r.c[swp(ai->first,idx_t())] += ai->second * b.c0;
			for (auto bi = b.c.cbegin(); bi != b.c.end(); ++bi)
				r.c[swp(ai->first, bi->first)] += ai->second * bi->second;
		}

		return r;
	}

private:
	c0_t c0;
	std::map<idx_t,cn_t> c;

	static constexpr sweep_mul swp = sweep_mul();

	friend void bla();
};

struct max_swp {
	template <typename idx_t>
	const idx_t & operator()(const idx_t &a, const idx_t &b) const
	{ return std::max(a,b); }
};

void bla()
{
	REAL r = 0;
	TM2<1,max_swp> x = r, y = r, z;
	z += x -= y;
	z *= x;
	z = x + y;
	z = x - y;
	z = -y;
	z = x * y;
	cout << z.c.begin()->second;
}

}

#endif
