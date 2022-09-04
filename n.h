
#ifdef N_MT
#include "threadpool.h"
#endif

#ifdef N_SMALLVECTOR
#include "smallvector.h"
#else
#endif

#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <vector>
#include <queue>
#include <tuple>


template <int static_digits = 100>
class N
{
	private:

		typedef unsigned char D;
		typedef signed char SD;

#ifdef N_DEBUG
		std::string sr;
		void srx()
		{
			sr = s(b, true);
		}
#endif

#ifdef N_SMALLVECTOR
#include "smallvector.h"
		using n_v = llvm_vecsmall::SmallVector<D,static_digits>;
		template <typename X> using n_vv = llvm_vecsmall::SmallVector<X,100>;
#else
		using n_v = std:vector<D>;
		template <typename X> using n_vv = std::vector<X>;
#endif


		D b = 10;
#ifdef N_SMALLVECTOR
		n_v n;
#else
		std::vector<D> n;
#endif
		bool m = false;
		unsigned long long dot = 0;
#ifdef _WIN64
		typedef signed long long ssize_t;
#else
		typedef signed long ssize_t;
#endif

public:

	void clear()
	{
		n.clear();
		m = false;
		dot = 0;
	}

	void RemoveZeroes()
	{
		while (n.size() > 0)
		{
			if (n[n.size() - 1] != 0)
				break;
			n.erase(n.end() - 1);
		}
#ifdef N_DEBUG
		srx();
#endif
	}

	static N<static_digits> rand(unsigned long long to = 100000000LL)
	{
		static std::random_device rd;
		static std::mt19937 mt(rd());
		std::uniform_int_distribution<unsigned long long> dist(0LL,to);
		return dist(mt);
	}

	bool IsEven() const
	{
		if (NumDigits() == 0)
			return true;
		if ((n[0] % 2) == 0)
			return true;
		return false;
	}

	bool IsZero() const
	{
		return NumDigits() == 0;
	}

	size_t NumDigits() const
	{
		return n.size();
	}

	D operator[](size_t idx) const
	{
		if (n.size() <= idx)
			return 0;
		return n[idx];
	}

	void shl2(size_t t2)
	{
		while (t2 > 0)
		{
			*this *= 2LL;
			t2--;
		}
	}

	void shr2(size_t t2)
	{
		while (t2 > 0)
		{
			*this /= 2LL;
			t2--;
		}
	}


	N& operator <<=(size_t t2)
	{
		shl2(t2);
		return *this;
	}

	N& operator >>=(size_t t2)
	{
		shr2(t2);
		return *this;
	}

	N& operator &=(const N& t2)
	{
		*this = w_logical(*this, t2, 0);
		return *this;
	}

	N& operator |=(const N& t2)
	{
		*this = w_logical(*this, t2, 1);
		return *this;
	}
	N& operator ^=(const N& t2)
	{
		if (b == 2)
			*this = w_logical(*this, t2, 2);
		else
			*this = w_pow(*this, t2);
		return *this;
	}
	

	N<static_digits> cshl(int t) const
	{
		N<static_digits> a = *this;
		a.shl(t);
		return a;
	}

	N<static_digits> cshr(int t) const
	{
		N<static_digits> a = *this;
		a.shr(t);
		return a;
	}

	void shr(size_t t)
	{
		if (NumDigits() <= t)
			clear();
		else
		{
			for (size_t i = 0; i < t; i++)
				n.erase(n.end() - 1);
		}
#ifdef N_DEBUG
		srx();
#endif
	}
	void shl(size_t t)
	{
		D d = 0;
		for (size_t i = 0; i < t; i++)
			n.insert(n.begin(), d);
#ifdef N_DEBUG
		srx();
#endif
	}



	static bool Less(const N& n1, const N& n2)
	{
		if (n1.m != n2.m)
		{
			if (n1.m == true)
				return true;
			return false;
		}
		if (n1.n.size() < n2.n.size())
			return true;
		if (n1.n.size() > n2.n.size())
			return false;

		for (ssize_t j = n1.n.size() - 1; j >= 0; j--)
		{
			D d1 = n1.n[j];
			D d2 = n2.n[j];
			if (d1 < d2)
				return true;
			if (d1 > d2)
				return false;
		}
		return false;
	}

	static bool Equal(const N& n1, const N& n2)
	{
		if (n1.m == n2.m && n1.n == n2.n)
			return true;
		return false;
	}


	bool operator <(const N<static_digits>& n2) const
	{
		return Less(*this, n2);
	}
	bool operator >(const N<static_digits>& n2) const
	{
		return Less(n2, *this);
	}
	bool operator <=(const N<static_digits>& n2) const
	{
		return !Less(n2, *this);
	}
	bool operator >=(const N<static_digits>& n2) const
	{
		return !Less(*this, n2);
	}
	bool operator !=(const N<static_digits>& n2) const
	{
		return !Equal(*this, n2);
	}
	bool operator ==(const N<static_digits>& n2) const
	{
		return Equal(*this, n2);
	}
	N<static_digits>& operator ++()
	{
		N<static_digits> nx((signed long long)1);
		return operator +=(nx);
	}
	N<static_digits>& operator --()
	{
		N<static_digits> nx((signed long long)1);
		return operator -=(nx);
	}

	static N<static_digits> w_subx(const N<static_digits>& n1, const N<static_digits>& n2)
	{
		if (n1.m != n2.m)
			return w_add(n1, n2.negative());

		if (n1.absolute() < n2.absolute())
			return w_subx(n2, n1).negative();

		if (n2.IsZero())
			return n1;
		if (n1.IsZero())
			return n2.negative();

		N<static_digits> n;
		n.ChangeInternalBase(n1.b);
		n.m = n1.m;
		int carry = 0;

		for (size_t i = 0; i < n1.NumDigits() || i < n2.NumDigits(); i++)
		{
			SD sum = n1[i] - n2[i] + carry;
			carry = 0;
			if (sum < 0)
			{
				sum = n1.b + sum;
				carry = -1;
			}
			n.n.push_back(sum);
#ifdef N_DEBUG
			n.srx();
#endif
		}
		n.n.push_back(carry);
		n.RemoveZeroes();
		return n;
	}


	static N<static_digits> w_pow(const N& n1, const N& n2)
	{
		N<static_digits> z = n1;
		if (n2 == 0ll)
			return 1ll;
		if (n2 == 1ll)
			return n1;
		if (n1 == 1ll)
			return 1ll;

		for (N<static_digits> j = 1ll; j < n2; ++j)
			z *= n1;

#ifdef N_DEBUG
		z.srx();
#endif
		return z;
	}


#ifdef N_MT
	static N<static_digits> w_add2(n_vv<N<static_digits>>& n, ThreadPool& pool)
	{
		if (n.size() == 0)
			return 0LL;
		if (n.size() == 1)
			return n[0];
		if (n.size() == 2)
		{
			N res = n[0];
			res += n[1];
			return res;
		}
		struct Z
		{
			N<static_digits>* n1;
			N<static_digits>* n2;
			N<static_digits>* res;
			std::future<void> v;
		};
		n_vv<Z> z(n.size());
		n_vv<N<static_digits>> res(n.size() / 2);

		for (size_t i = 0; i < n.size(); i += 2)
		{
			if (i == (n.size() - 1))
				break; // odd number of additions

			auto a = [](Z* z)
			{
				*z->res = w_add(*z->n1, *z->n2);
			};
			Z& zz = z[i];
			zz.n1 = &n[i];
			zz.n2 = &n[i + 1];
			zz.res = &res[i / 2];

			zz.v = pool.enqueue([&](Z* z) { a(z); }, &zz);
		}
		for (size_t i = 0; i < n.size(); i += 2)
		{
			if (z[i].v.valid())
				z[i].v.get();
		}
		if (n.size() % 2)
			res.push_back(n[n.size() - 1]);
		return w_add2(res, pool);
	}


	static N<static_digits> w_pow2(const N<static_digits>& n1, const N<static_digits>& n2, ThreadPool& pool)
	{
		N<static_digits> z = n1;
		if (n2 == 0ll)
			return 1ll;
		if (n2 == 1ll)
			return n1;
		if (n1 == 1ll)
			return 1ll;

		for (N<static_digits> j = 1ll; j < n2; ++j)
			z = w_mul2(z, n1, pool);
		return z;
	}

	static N<static_digits> w_mul2(const N<static_digits>& n1, const N<static_digits>& n2, ThreadPool& pool)
	{
		size_t muls = n1.NumDigits() * n2.NumDigits();
		n_vv<N<static_digits>> a;
		a.reserve(muls);
		for (size_t i = 0; i < n1.NumDigits(); i++)
		{
			for (size_t ii = 0; ii < n2.NumDigits(); ii++)
			{
				N<static_digits> rr;
				D d1 = n1[i];
				D d2 = n2[ii];
				unsigned long long r = d1 * d2;
				rr = r;
				rr.shl(ii + i);
				a.push_back(rr);
			}
		}
		return w_add2(a, pool);
	}
#endif

	static std::tuple<N, N> w_div(const N& n1, const N& n2, bool NoChangeBase = false)
	{
		if (n1.b != n2.b && NoChangeBase == false)
			return w_div(n1.b, n2.tobase(n1.b));
		if (n2 > n1)
		{
			N res = n1;
			return std::make_tuple<N, N>(0LL, std::forward<N>(res));
		}
		if (n2 == n1)
			return std::make_tuple<N, N>(1LL, 0LL);

		N rem = n1;
		N res;
		res.ChangeInternalBase(n1.b);

		for (;;)
		{
			auto nd2 = n2.NumDigits();
			auto upper = rem.upperpart(nd2);
			if (upper < n2)
			{
				nd2++;
				upper = rem.upperpart(nd2);
				if (upper < n2)
				{
					// End...
					return std::make_tuple<N, N>(std::forward<N>(res), std::forward<N>(rem));
				}
			}

			D js = 9;
			N m1;
			for (; js >= 1; js--)
			{
				m1 = w_mul(n2, js);
				if (m1 < upper)
					break;
			}

			res.n.insert(res.n.begin(), js);
#ifdef N_DEBUG
			res.srx();
#endif
			upper -= m1;
			upper.shl(rem.NumDigits() - nd2);
			upper += rem.lowerpart(rem.NumDigits() - nd2);
			rem = upper;
		}
	}

	static N<static_digits> w_mul(const N& n1, const N& n2)
	{
		if (n1.b != n2.b)
			return w_mul(n1, n2.tobase(n1.b));
		N<static_digits> n;
		n.ChangeInternalBase(n1.b);
		n_vv<N<static_digits>> addiz;
		for (size_t i = 0; i < n1.n.size(); i++)
		{
			D d1 = n1.n[i];
			N<static_digits> addi;
			addi.n.reserve(i + n2.n.size());
			for (size_t j = 0; j < i; j++)
				addi.n.push_back(0);
			D carry = 0;
			for (size_t y = 0; y < n2.n.size(); y++)
			{
				D d2 = n2.n[y];
				D dm = (d1 * d2) + carry;
				carry = 0;
				carry = dm / n1.b;
				dm %= n1.b;
				addi.n.push_back(dm);
#ifdef N_DEBUG
				addi.srx();
#endif
			}
			addi.n.push_back(carry);
			addi.RemoveZeroes();
			addiz.push_back(addi);
		}
		for (auto& a : addiz)
			n += a;
		if (n1.m != n2.m)
			n.m = true;
		return n;
	}

	static N<static_digits> w_logical(const N<static_digits>& n1, const N<static_digits>& n2, int x)
	{
		if (n1.b != 2)
			return w_logical(n1.tobase(2), n2, x);
		if (n2.b != 2)
			return w_logical(n1, n2.tobase(2), x);

		N<static_digits> n;
		n.ChangeInternalBase(2);
		n.n.reserve(std::max(n1.NumDigits(), n2.NumDigits()));

		for (size_t i = 0; i < n1.NumDigits() || i < n2.NumDigits(); i++)
		{
			D sum = 0;
			if (x == 0) sum = n1[i] & n2[i];
			if (x == 1) sum = n1[i] | n2[i];
			if (x == 2) sum = n1[i] ^ n2[i];
			n.n.push_back(sum);
#ifdef N_DEBUG
			n.srx();
#endif
		}
		n.RemoveZeroes();
		return n;
	}


	static N<static_digits> w_add(const N<static_digits>& n1, const N<static_digits>& n2)
	{
		if (n1.b != n2.b)
			return w_add(n1, n2.tobase(n1.b));
		if (n1.m != n2.m)
		{
			if (n1.m)
				return w_subx(n2, n1.negative());
			return w_subx(n1, n2.negative());
		}
		if (n1.n.empty()) return n2;
		if (n2.n.empty()) return n1;

		N<static_digits> n;
		n.ChangeInternalBase(n1.b);
		n.n.reserve(std::max(n1.NumDigits(), n2.NumDigits()));
		D carry = 0;

		if (n1.m && n2.m)
			n.m = true;

		size_t j = 0;
		for (size_t i = 0; i < n1.NumDigits() || i < n2.NumDigits(); i++)
		{
			j = i;
			D sum = n1[i] + n2[i] + carry;
			carry = 0;
			if (sum >= n1.b)
			{
				carry = 1;
				sum -= n1.b;
			}
			n.n.push_back(sum);
#ifdef N_DEBUG
			n.srx();
#endif
		}
		n.n.push_back(carry);
		n.RemoveZeroes();
		return n;
	}


	N<static_digits>& operator += (const N& nn)
	{
		*this = w_add(*this, nn);
		return *this;
	}
	N<static_digits>& operator -= (const N& nn)
	{
		*this = w_add(*this, nn.negative());
		return *this;
	}
	N<static_digits>& operator *= (const N& nn)
	{
		*this = w_mul(*this, nn);
		return *this;
	}
	N<static_digits>& operator /= (const N& nn)
	{
		*this = std::get<0>(w_div(*this, nn));
		return *this;
	}
	N<static_digits>& operator %= (const N& nn)
	{
		*this = std::get<1>(w_div(*this, nn));
		return *this;
	}

	void ChangeInternalBase(D nb = 10)
	{
		b = nb;
	}

	N<static_digits> tobase(D nb) const
	{
		N<static_digits> n2(s(nb).c_str(), nb);
		n2.m = m;
		return n2;
	}

	N<static_digits>& ParseBase(const char* a1, unsigned long long B = 16)
	{
		N<static_digits> res;
		res.b = b;
		unsigned long long  k = 0;
		for (ssize_t i = strlen(a1) - 1; i >= 0; i--)
		{
			unsigned long long j = a1[i];
			if (j >= 'a')
				j -= ('a' - 10);
			else
				if (j >= 'A')
					j -= ('A' - 10);
				else
					if (j >= '0')
						j -= '0';

			if (j >= B)
				break; // duh
			auto p = w_pow(B, k);
			p *= j;
			res += p;
			k++;
		}
		operator =(res);
		return *this;
	}


	N<static_digits> lowerpart(size_t digs)
	{
		if (digs > NumDigits())
			digs = NumDigits();
		N<static_digits> n2 = *this;
		auto b1 = n.begin();
		auto e1 = n.begin() + digs;
		n_v nv(b1, e1);
		n2.n = nv;
#ifdef N_DEBUG
		n2.srx();
#endif
		return n2;
	}

	D upperdigit() const
	{
		return n[n.size() - 1];
	}

	N<static_digits> upperpart(size_t digs)
	{
		if (digs > NumDigits())
			digs = NumDigits();
		N<static_digits> n2 = *this;
		auto b1 = n.end();
		auto e1 = n.end() - digs;
		n_v nv(e1, b1);
		n2.n = nv;
#ifdef N_DEBUG
		n2.srx();
#endif
		return n2;
	}

	N<static_digits> negative() const
	{
		N<static_digits> nx = *this;
		nx.m = !nx.m;
#ifdef N_DEBUG
		nx.srx();
#endif
		return nx;
	}

	N<static_digits> absolute() const
	{
		N nn = *this;
		nn.m = false;
#ifdef N_DEBUG
		nn.srx();
#endif
		return nn;
	}

	N<static_digits> HexNeg(int bits) const
	{
		unsigned long long r = (unsigned long long)pow(2LL, bits);
		N<static_digits> a(r);
		//		a += 1ll;
		a -= absolute().tobase(10);
		return a.tobase(16);
	}

	std::string s(unsigned long long base = -1, bool Debug = false) const
	{
		if (base == (unsigned long long) - 1)
			base = b;
		if (base == 0)
			return "";
		std::string a;
		if (base == b)
		{
			if (n.empty())
				return "0";

			auto cdot = dot;
			if (m)
				a += "-";
			for (ssize_t j = (n.size() - 1); j >= 0; j--)
			{
				auto d = n[j];
				if (d >= 10)
					a += (char)(d + 'A' - 10);
				else
					a += (char)(d + 0x30u);
			}
			return a;
		}
		if (Debug)
			return "";

		N<static_digits> e = *this;
		if (m)
			a += "-";
		N<static_digits> rem;
		for (;;)
		{
			auto rx = w_div(e, base, true);
			auto d = abs(atoi(std::get<1>(rx).s().c_str()));
			if (d >= 10)
				d = (d - 10) + 'A';
			else
				d += 0x30;
			a += (char)d;
			e = std::get<0>(rx);
			rem = std::get<1>(rx);
			if (e.IsZero())
				break;
		}
		std::reverse(a.begin(), a.end());
		return a;
	}

	void Set(unsigned long long a = 0)
	{
		clear();
		for (;;)
		{
			D d = a % b;
			n.push_back(d);
			a /= b;
			if (!a)
				break;
		}
		RemoveZeroes();
#ifdef N_DEBUG
		srx();
#endif
	}

	void Set(signed long long a = 0)
	{
		clear();
		if (a < 0)
		{
			Set((unsigned long long) - a);
			m = true;
		}
		else
			Set((unsigned long long)a);
#ifdef N_DEBUG
		srx();
#endif
	}


	void Set(const char* a)
	{
		clear();
		if (!a)
			return;
		bool bdot = false;
		for (size_t i = 0; i < strlen(a); i++)
		{
			char a1 = a[i];
			if (i == 0 && a1 == '-')
			{
				m = true;
				continue;
			}
			if (i == 0 && a1 == '+')
			{
				m = false;
				continue;
			}
			if (a1 == '.')
			{
				bdot = true;
				continue;
			}
			int jn = 0;

			if (a1 < '0')
				break;
			if (a1 >= 'a')
				jn = a1 - ('a' - 10);
			else
				if (a1 >= 'A')
					jn = a1 - ('A' - 10);
				else
					if (a1 >= '0')
						jn = a1 - ('0');
			if (jn >= b)
				break;

			n.insert(n.begin(), jn);
			if (bdot)
				dot++;
		}
		RemoveZeroes();
#ifdef N_DEBUG
		srx();
#endif
	}


	N(const char* a, D ba = 10)
	{
		b = ba;
		Set(a);
	}
	N(unsigned long long a)
	{
		Set(a);
	}
	N(signed long long a = 0)
	{
		Set(a);
	}
	N(unsigned char a)
	{
		Set((unsigned long long)a);
	}
	N<static_digits>& operator=(unsigned long long a)
	{
		Set(a);
		return *this;
	}
	N<static_digits>& operator=(unsigned char a)
	{
		return operator=((unsigned long long)a);
	}
	N<static_digits>& operator=(signed long long a)
	{
		Set(a);
		return *this;
	}
	N<static_digits>& operator=(const char* a)
	{
		Set(a);
		return *this;
	}


};

template <int static_digits = 100>
N<static_digits> operator + (const N<static_digits>& lhs, const N<static_digits>& rhs)
{
	N<static_digits> n = lhs;
	n += rhs;
	return n;
}

template <int static_digits = 100>
N<static_digits> operator - (const N<static_digits>& lhs, const N<static_digits>& rhs)
{
	N<static_digits> n = lhs;
	n -= rhs;
	return n;
}

template <int static_digits = 100>
N<static_digits> operator * (const N<static_digits>& lhs, const N<static_digits>& rhs)
{
	N<static_digits> n = lhs;
	n *= rhs;
	return n;
}
template <int static_digits = 100>
N<static_digits> operator / (const N<static_digits>& lhs, const N<static_digits>& rhs)
{
	N<static_digits> n = lhs;
	n /= rhs;
	return n;
}

template <int static_digits = 100>
N<static_digits> operator % (const N<static_digits>& lhs, const N<static_digits>& rhs)
{
	N<static_digits> n = lhs;
	n %= rhs;
	return n;
}

template <int static_digits = 100>
N<static_digits> operator | (const N<static_digits>& lhs, const N<static_digits>& rhs)
{
	N<static_digits> n = lhs;
	n |= rhs;
	return n;
}

template <int static_digits = 100>
N<static_digits> operator & (const N<static_digits>& lhs, const N<static_digits>& rhs)
{
	N<static_digits> n = lhs;
	n &= rhs;
	return n;
}

template <int static_digits = 100>
N<static_digits> operator ^ (const N<static_digits>& lhs, const N<static_digits>& rhs)
{
	N<static_digits> n = lhs;
	n ^= rhs;
	return n;
}
