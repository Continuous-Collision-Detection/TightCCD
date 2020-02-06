#ifndef BAC_TIGHTBOUND_H
#define BAC_TIGHTBOUND_H

#include <vec.h>

namespace bsc_tightbound {

//#define EPSILON 0.00000000000000022204460492503131
#define BOUND_UNIT 0.00000000000000022204460492503136

	enum ReturnValue{
		RETURN_ZERO,
		RETURN_FALSE,
		RETURN_TRUE
	};


	template<typename T>
	class floating{
	public:
		T v;
		T sigma;
	public:
		inline floating() :v(T(0)), sigma(T(0)){}
		inline floating(T _v) : v(_v), sigma(T(0)){}
		inline floating(T _v, T _sigma) :v(_v), sigma(_sigma){
			assert(sigma >= 0);
		}
		inline int sign(){
			if (v > sigma)
				return 1;
			else if (v < (-sigma))
				return -1;
			else
				return 0;
			/*if (v == 0){
				return 0;
			}
			else{
				if (v >= sigma)
					return 1;
				else if (v <= (-sigma))
					return -1;
				else
					return 0;
			}*/
		}
		inline int sign() const {
			if (v > sigma)
				return 1;
			else if (v < (-sigma))
				return -1;
			else
				return 0;
			/*if (v == 0){
				return 0;
			}
			else{
				if (v >= sigma)
					return 1;
				else if (v <= (-sigma))
					return -1;
				else
					return 0;
			}*/
		}
		inline bool operator == (const T t){
			if (v == t)
				return true;
			return false;
		}
		inline floating& operator += (const floating& rhs){
			v += rhs.v;
			/*if (sizeof(T) == 4)
				sigma = sigma + rhs.sigma + fabs(v)*FLT_EPSILON / (1 - FLT_EPSILON);
			else if (sizeof(T) == 8)*/
			sigma = sigma + rhs.sigma + fabs(v)*BOUND_UNIT;//DBL_EPSILON / (1 - DBL_EPSILON);
			return *this;
		}
		inline floating& operator -= (const floating& rhs){
			v -= rhs.v;
			/*if (sizeof(T) == 4)
				sigma = sigma + rhs.sigma + fabs(v)*FLT_EPSILON / (1 - FLT_EPSILON);
			else if (sizeof(T) == 8)*/
			sigma = sigma + rhs.sigma + fabs(v)*BOUND_UNIT;// DBL_EPSILON / (1 - DBL_EPSILON);
			return *this;
		}
		inline floating& operator *= (const floating& rhs){
			sigma = sigma*rhs.sigma + fabs(v)*rhs.sigma + fabs(rhs.v)*sigma;
			v *= rhs.v;
			/*if (sizeof(T) == 4)
				sigma += fabs(v)*FLT_EPSILON / (1 - FLT_EPSILON);
			else if (sizeof(T) == 8)*/
			sigma += fabs(v)*BOUND_UNIT;// DBL_EPSILON / (1 - DBL_EPSILON);
			return *this;
		}
		inline floating operator + (const floating& other) const{
			T _v = v + other.v;
			T _sigma = sigma + other.sigma + fabs(_v)*BOUND_UNIT;// DBL_EPSILON / (1 - DBL_EPSILON);
			/*if (sizeof(T) == 4){
				_sigma += fabs(_v)*FLT_EPSILON / (1 - FLT_EPSILON);
			}
			else if (sizeof(T) == 8){
				_sigma += fabs(_v)*DBL_EPSILON / (1 - DBL_EPSILON);
			}*/
			return floating(_v, _sigma);
		}
		inline floating operator - (const floating& other) const{
			T _v = v - other.v;
			T _sigma = sigma + other.sigma + fabs(_v)*BOUND_UNIT;//DBL_EPSILON / (1 - DBL_EPSILON);
			/*if (sizeof(T) == 4){
				_sigma += fabs(_v)*FLT_EPSILON / (1 - FLT_EPSILON);
			}
			else if (sizeof(T) == 8){
				_sigma += fabs(_v)*DBL_EPSILON / (1 - DBL_EPSILON);
			}*/
			return floating(_v, _sigma);
		}
		inline floating operator * (const floating& other) const{
			T _v = v * other.v;
			T _sigma = sigma*other.sigma + fabs(v)*other.sigma + fabs(other.v)*sigma + fabs(_v)*BOUND_UNIT;// DBL_EPSILON / (1 - DBL_EPSILON);
			/*if (sizeof(T) == 4){
				_sigma += fabs(_v)*FLT_EPSILON / (1 - FLT_EPSILON);
			}
			else if (sizeof(T) == 8){
				_sigma += fabs(_v)*DBL_EPSILON / (1 - DBL_EPSILON);
			}*/
			return floating(_v, _sigma);
		}
		inline floating operator - () const
		{
			return floating(-v, sigma);
		}
	};


	template <class T>
	class bcrv {
	public:
		floating<T> k0, k1, k2, k3;
		floating<T> kk0, kk1, kk2;
		int ct;
		int sign_k0, sign_k1, sign_k2, sign_k3;
		int sign_kk0, sign_kk1, sign_kk2;
	public:
		bcrv(const floating<T>& _k0, const floating<T>& _k1, const floating<T>& _k2, const floating<T>& _k3,
			const floating<T>& _kk0, const floating<T>& _kk1, const floating<T>& _kk2, int _ct)
		{
			this->k0 = _k0;
			this->k1 = _k1;
			this->k2 = _k2;
			this->k3 = _k3;
			this->kk0 = _kk0;
			this->kk1 = _kk1;
			this->kk2 = _kk2;
			this->ct = _ct;

			sign_k0 = k0.sign();
			sign_k1 = k1.sign();
			sign_k2 = k2.sign();
			sign_k3 = k3.sign();
			sign_kk0 = kk0.sign();
			sign_kk1 = kk1.sign();
			sign_kk2 = kk2.sign();
		}
	};


	template <class T>
	inline ReturnValue diffSign(const floating<T>& a, const floating<T>& b)
	{
		if (a.sign() == 0 || b.sign() == 0){
			return RETURN_ZERO;
		}
		else if ((a.sign() > 0 && b.sign() < 0) || (a.sign() < 0 && b.sign() > 0)){
			return RETURN_TRUE;
		}
		else{
			return RETURN_FALSE;
		}
	}


	template <class T>
	inline ReturnValue sameSign(const floating<T>& a, const floating<T>& b)
	{
		// || (a == 0 && b == 0);
		if (a.sign() == 0 || b.sign() == 0){
			return RETURN_ZERO;
		}
		else if ((a.sign() > 0 && b.sign() > 0) || (a.sign() < 0 && b.sign() < 0)){
			return RETURN_TRUE;
		}
		else{
			return RETURN_FALSE;
		}
	}


	template<unsigned int N, class T>
	void make_vector(const Vec<N, T> v, Vec<N, floating<T>>& out){
		for (int i = 0; i < N; i++){
			out[i] = floating<T>(v[i], 0);
		}
	}

	//					| a  b |
	//  return		| c  d | = a*d - b*c
	template <class T>
	inline floating<T> det2x2(const floating<T> &a, const floating<T> &b, const floating<T> &c, const floating<T> &d)
	{
		return a*d - b*c;
	}


	template <class T>
	inline floating<T> _evaluateBezier2(const floating<T> &p0, const floating<T> &p1, const floating<T> &p2,
		const floating<T> &t, const floating<T> &s)
	{
		floating<T> s2 = s*s;
		floating<T> t2 = t*t;

		return p0*s2 + p1*floating<T>(2.0, 0)*s*t + p2*t2;
	}


	template <class T>
	inline floating<T> evaluateBezier1(const floating<T> &p0, const floating<T> &p1,
		const floating<T> &t)
	{
		floating<T> s = floating<T>(1.0, 0) - t;
		return p0*s + p1*t;
	}


	template <class T>
	inline floating<T> evaluateBezier2(const floating<T> &p0, const floating<T> &p1, const floating<T> &p2,
		const floating<T> &t)
	{
		floating<T> s = floating<T>(1.0, 0) - t;
		return _evaluateBezier2(p0, p1, p2, t, s);
	}


	template <class T>
	inline floating<T> _evaluateBezier(const floating<T> &p0, const floating<T> &p1, const floating<T> &p2, const floating<T> &p3,
		const floating<T> &t, const floating<T> &s)
	{
		floating<T> s2 = s*s;
		floating<T> s3 = s2*s;
		floating<T> t2 = t*t;
		floating<T> t3 = t2*t;

		return p0*s3 + p1*floating<T>(3.0, 0)*s2*t + p2*floating<T>(3.0, 0)*s*t2 + p3*t3;
	}


	template <class T>
	inline floating<T> evaluateBezier(const floating<T> &p0, const floating<T> &p1, const floating<T> &p2, const floating<T> &p3,
		const floating<T> &t)
	{
		floating<T> s = floating<T>(1.0, 0) - t;
		return _evaluateBezier(p0, p1, p2, p3, t, s);
	}


	template <class T>
	inline floating<T> _evaluateBezier4(const floating<T> &p0, const floating<T> &p1, const floating<T> &p2, const floating<T> &p3, const floating<T> &p4,
		const floating<T> &t, const floating<T> &s)
	{
		floating<T> s2 = s*s;
		floating<T> s3 = s2*s;
		floating<T> s4 = s2*s2;
		floating<T> t2 = t*t;
		floating<T> t3 = t2*t;
		floating<T> t4 = t2*t2;

		return p0*s4 + p1 * 4 * s3*t + p2 * 6 * s2*t2 + p3 * 4 * s*t3 + p4*t4;
	}


	template <class T>
	inline floating<T> evaluateBezier4(const floating<T> &p0, const floating<T> &p1, const floating<T> &p2, const floating<T> &p3, const floating<T> &p4,
		const floating<T> &t)
	{
		floating<T> s = floating<T>(1.0, 0) - t;
		return _evaluateBezier4(p0, p1, p2, p3, p4, t, s);
	}

	// from polynomial decomposition theorem
	template <class T>
	inline bool bezierDecomposition(const floating<T> &k0, const floating<T> &k1, const floating<T> &k2, const floating<T> &k3,
		const floating<T> &j0, const floating<T> &j1, const floating<T> &j2,
		floating<T> &m0, floating<T> &m1, floating<T> &n0, floating<T> &n1)
	{
		floating<T> A = (j1 - j2)*floating<T>(2.0, 0);
		floating<T> B = j0 - j2;
		floating<T> C = k2*floating<T>(3.0, 0) - k3*floating<T>(2.0, 0) - k0;
		floating<T> D = k1*floating<T>(3.0, 0) - k0*floating<T>(2.0, 0) - k3;
		floating<T> E = j2 - j0;
		floating<T> F = (j1 - j0)*floating<T>(2.0, 0);

		floating<T> tt = det2x2(A, B, E, F);

		//here has some problem
		if (tt == 0 || (j0 + j2 - j1*floating<T>(2.0, 0) == 0)){//(fabs(tt)<epsilon){//(tt.is_certainly_zero()) {
			//printf("bscfloat-2@@@det = 0.\n");
			return false;
		}

		/*
		m0 = det2x2(A, B, C, D)/tt;
		m1 = det2x2(F, E, D, C)/tt;
		n0 = k0-m0*j0;
		n1 = k3-m1*j2;
		*/
		m0 = det2x2(A, B, C, D);
		m1 = det2x2(F, E, D, C);
		n0 = k0*tt - m0*j0;
		n1 = k3*tt - m1*j2;

		return true;
	}


	template <class T>
	inline Vec<3, floating<T>> norm(const Vec<3, floating<T>> &p1, const Vec<3, floating<T>> &p2, const Vec<3, floating<T>> &p3)
	{
		return cross(p2 - p1, p3 - p1);
	}


	template <class T>
	inline bool DNF_Culling(
		const Vec<3, floating<T>> &a0, const Vec<3, floating<T>> &b0, const Vec<3, floating<T>> &c0, const Vec<3, floating<T>> &d0,
		const Vec<3, floating<T>> &a1, const Vec<3, floating<T>> &b1, const Vec<3, floating<T>> &c1, const Vec<3, floating<T>> &d1)
	{
		Vec<3, floating<T>> n0 = norm(a0, b0, c0);
		Vec<3, floating<T>> n1 = norm(a1, b1, c1);
		Vec<3, floating<T>> delta = norm(a1 - a0, b1 - b0, c1 - c0);
		Vec<3, floating<T>> nX = (n0 + n1 - delta)*floating<T>(0.5, 0);

		Vec<3, floating<T>> pa0 = d0 - a0;
		Vec<3, floating<T>> pa1 = d1 - a1;

		floating<T> A = dot(n0, pa0);
		floating<T> B = dot(n1, pa1);
		floating<T> C = dot(nX, pa0);
		floating<T> D = dot(nX, pa1);
		floating<T> E = dot(n1, pa0);
		floating<T> F = dot(n0, pa1);

		floating<T> p0 = A;
		floating<T> p1 = C*floating<T>(2.0, 0) + F;
		floating<T> p2 = D*floating<T>(2.0, 0) + E;
		floating<T> p3 = B;

		/*int sign_p0 = p0.sign();
		int sign_p1 = p1.sign();
		int sign_p2 = p2.sign();
		int sign_p3 = p3.sign();*/

		if (p0.sign() > 0 && p1.sign() > 0 && p2.sign() > 0 && p3.sign() > 0)
			return false;

		if (p0.sign() < 0 && p1.sign() < 0 && p2.sign() < 0 && p3.sign() < 0)
			return false;

		return true;
	}


	template <class T>
	inline bool getBezier(
		const Vec<3, floating<T>> &a0, const Vec<3, floating<T>> &b0, const Vec<3, floating<T>> &c0, const Vec<3, floating<T>> &d0,
		const Vec<3, floating<T>> &a1, const Vec<3, floating<T>> &b1, const Vec<3, floating<T>> &c1, const Vec<3, floating<T>> &d1,
		floating<T> &p0, floating<T> &p1, floating<T> &p2, floating<T> &p3,
		Vec<3, floating<T>> &n0, Vec<3, floating<T>> &n1, Vec<3, floating<T>> &delta, Vec<3, floating<T>> &nX)
	{
		n0 = norm(a0, b0, c0);
		n1 = norm(a1, b1, c1);
		delta = norm(a1 - a0, b1 - b0, c1 - c0);
		nX = (n0 + n1 - delta)*floating<T>(0.5, 0);

		Vec<3, floating<T>> pa0 = d0 - a0;
		Vec<3, floating<T>> pa1 = d1 - a1;

		/*floating<T> A = dot(n0, pa0);
		floating<T> B = dot(n1, pa1);
		floating<T> C = dot(nX, pa0);
		floating<T> D = dot(nX, pa1);
		floating<T> E = dot(n1, pa0);
		floating<T> F = dot(n0, pa1);*/

		p0 = dot(n0, pa0)*floating<T>(3.0, 0);//A*floating<T>(3.0, 0);
		p1 = dot(nX, pa0)*floating<T>(2.0, 0) + dot(n0, pa1);//C*floating<T>(2.0, 0) + F;
		p2 = dot(nX, pa1)*floating<T>(2.0, 0) + dot(n1, pa0);//D*floating<T>(2.0, 0) + E;
		p3 = dot(n1, pa1)*floating<T>(3.0, 0);//B*floating<T>(3.0, 0);

		if (p0.sign() > 0 && p1.sign() > 0 && p2.sign() > 0 && p3.sign() > 0)
			return false;
		//if (sign(p0)>0 && sign(p1)>0 && sign(p2)>0 && sign(p3)>0)
		/*if (p0.is_certainly_positive() && p1.is_certainly_positive() &&
		p2.is_certainly_positive() && p3.is_certainly_positive())
		return false;*/

		if (p0.sign() < 0 && p1.sign() < 0 && p2.sign() < 0 && p3.sign() < 0)
			return false;
		//if (sign(p0)<0 && sign(p1)<0 && sign(p2)<0 && sign(p3)<0)
		/*if (p0.is_certainly_negative() && p1.is_certainly_negative() &&
		p2.is_certainly_negative() && p3.is_certainly_negative())
		return false;*/

		return true;
	}

	template <class T>
	inline void getBezier4(
		const Vec<3, floating<T>> &a0, const Vec<3, floating<T>> &b0, const Vec<3, floating<T>> &c0, const Vec<3, floating<T>> &d0,
		const Vec<3, floating<T>> &a1, const Vec<3, floating<T>> &b1, const Vec<3, floating<T>> &c1, const Vec<3, floating<T>> &d1,
		const Vec<3, floating<T>> &n0, const Vec<3, floating<T>> &n1, const Vec<3, floating<T>> &deltaN, const Vec<3, floating<T>> &nX,
		floating<T> &l0, floating<T> &l1, floating<T> &l2, floating<T> &l3, floating<T> &l4, int which, bool ee_test)
	{
		/*Vec<3, floating<T>> n0 = norm(a0, b0, c0);
		Vec<3, floating<T>> n1 = norm(a1, b1, c1);
		Vec<3, floating<T>> deltaN = norm(a1 - a0, b1 - b0, c1 - c0);
		Vec<3, floating<T>> nX = (n0 + n1 - deltaN)*floating<T>(0.5,0);*/

		Vec<3, floating<T>> m0, m1, deltaM, mX;

		if (which == 0) { // (bt-pt) x (ct-pt) . nt
			m0 = norm(d0, b0, c0);
			m1 = norm(d1, b1, c1);
			deltaM = norm(d1 - d0, b1 - b0, c1 - c0);
		}
		else if (which == 1) { // ct-pt x at-pt . nt
			m0 = norm(d0, c0, a0);
			m1 = norm(d1, c1, a1);
			deltaM = norm(d1 - d0, c1 - c0, a1 - a0);
		}
		else if (which == 2) {// at-pt x bt-pt .nt
			m0 = norm(d0, a0, b0);
			m1 = norm(d1, a1, b1);
			deltaM = norm(d1 - d0, a1 - a0, b1 - b0);
		}
		else
			printf("bscfloat@@@Imposible be here!");

		mX = (m0 + m1 - deltaM)*floating<T>(0.5, 0);

		l0 = dot(m0, n0)*floating<T>(6.0, 0);
		l1 = (dot(m0, nX) + dot(mX, n0))*floating<T>(3.0, 0);
		l2 = dot(m0, n1) + dot(mX, nX)*floating<T>(4.0, 0) + dot(m1, n0);
		l3 = (dot(mX, n1) + dot(m1, nX))*floating<T>(3.0, 0);
		l4 = dot(m1, n1)*floating<T>(6.0, 0);

		if (which == 2 && ee_test) {
			l0 = -l0, l1 = -l1, l2 = -l2, l3 = -l3, l4 = -l4;
		}
	}

	template <class T>
	inline int bezierClassification(const floating<T>& k0, const floating<T>& k1, const floating<T>& k2, const floating<T>& k3,
		floating<T> &kk0, floating<T> &kk1, floating<T> &kk2)
	{
		if (k0.sign() > 0 && k1.sign() < 0 && k2.sign() < 0 && k3.sign() < 0){
			return 0;
		}
		//if (sign(k0) > 0 && sign(k1) < 0 && sign(k2) < 0 && sign(k3) < 0)
		/*if (k0.is_certainly_positive() &&
		k1.is_certainly_negative() &&
		k2.is_certainly_negative() &&
		k3.is_certainly_negative())
		return 0;*/

		if (k0.sign() < 0 && k1.sign() > 0 && k2.sign() > 0 && k3.sign() > 0){
			return 0;
		}
		//if (sign(k0) < 0 && sign(k1) > 0 && sign(k2) > 0 && sign(k3) > 0)
		/*if (k0.is_certainly_negative() &&
		k1.is_certainly_positive() &&
		k2.is_certainly_positive() &&
		k3.is_certainly_positive())
		return 0;*/

		if (k3.sign() > 0 && k1.sign() < 0 && k2.sign() < 0 && k0.sign() < 0){
			return 0;
		}
		//if (sign(k3) > 0 && sign(k1) < 0 && sign(k2) < 0 && sign(k0) < 0)
		/*if (k3.is_certainly_positive() &&
		k1.is_certainly_negative() &&
		k2.is_certainly_negative() &&
		k0.is_certainly_negative())
		return 0;*/

		if (k3.sign() < 0 && k1.sign() > 0 && k2.sign() > 0 && k0.sign() > 0){
			return 0;
		}
		//if (sign(k3) < 0 && sign(k1) > 0 && sign(k2) > 0 && sign(k0) > 0)
		/*if (k3.is_certainly_negative() &&
		k1.is_certainly_positive() &&
		k2.is_certainly_positive() &&
		k0.is_certainly_positive())
		return 0;*/

		// f'' = 6*(k2-2*k1+k0)*B^0_1 + 6*(k3-2*k2+k1)*B^1_1
		floating<T> a = k2 - k1*floating<T>(2.0, 0) + k0;
		floating<T> b = k3 - k2*floating<T>(2.0, 0) + k1;

		// f = 3*(k1-k0) B^2_0 + 3*(k2-k1)*B^2_1 + 3*(k3-k2)*B^2_2
		kk0 = k1 - k0;
		kk1 = k2 - k1;
		kk2 = k3 - k2;

		//printf("%.32f %.32f %.32f %.32f %.32f %.32f %.32f %.32f %.32f\n",
			//k0.sigma, k1.sigma, k2.sigma, k3.sigma, kk0.sigma, kk1.sigma, kk2.sigma, a.sigma, b.sigma);
		/*printf("%.32f %.32f %.32f %.32f %.32f %.32f",
			k0.sigma, k3.sigma, kk0.sigma, kk2.sigma, a.sigma, b.sigma);*/

		if (diffSign(a, b) != RETURN_FALSE) {
			return 2; // 1 inflexion
		}

		/*// f = 3*(k1-k0) B^2_0 + 3*(k2-k1)*B^2_1 + 3*(k3-k2)*B^2_2
		kk0 = k1 - k0;
		kk1 = k2 - k1;
		kk2 = k3 - k2;*/

		if (diffSign(kk0, kk2) != RETURN_FALSE)
			return 1; // no inflexion, 1 extreme
		else
			return 0; // no inflexion, no extreme
	}


	template <class T>
	inline int coplanarTest(bcrv<T> &c)
	{
		//printf("T-bound: k0:%.32f k1:%.32f k2:%.32f k3:%.32f kk0:%.32f kk1:%.32f kk2:%.32f\n", c.k0.sigma, c.k1.sigma, c.k2.sigma, c.k3.sigma, c.kk0.sigma, c.kk1.sigma, c.kk2.sigma);

		if (c.k0.sign() == 0 && c.k3.sign() == 0){
			return 2; //conservative operation
		}

		if (c.ct == 0) {// we only need to make sure sign(k0) != sign(k3)
			if (diffSign(c.k0, c.k3) != RETURN_FALSE){

				if (c.k0.sign() == 0 && c.k3.sign() != 0)
					c.sign_k0 = -c.sign_k3;
				else if (c.k0.sign() != 0 && c.k3.sign() == 0)
					c.sign_k3 = -c.sign_k0;

				return 1;
			}
			else
				return 0;
		}
		else {
			if (diffSign(c.k0, c.k3) == RETURN_TRUE)
				return 1;

			//ensure Y(0) and Y(1) have the same nonzero sign
			//int sign_k0, sign_k3;
			bool flag = false;

			if (c.k0.sign() == 0){
				c.sign_k0 = c.k3.sign();
				c.sign_k3 = c.k3.sign();
				flag = true;
			}
			else if (c.k3.sign() == 0){
				c.sign_k0 = c.k0.sign();
				c.sign_k3 = c.k0.sign();
				flag = true;
			}

			/*assert(c.sign_k0 == c.sign_k3);
			assert(c.sign_k0 != 0);*/

			/*if (c.kk0 + c.kk2 - c.kk1*T(2.0) == 0){//(fabs(c.kk0 + c.kk2 - c.kk1*T(2.0))<epsilon){//.is_certainly_zero()) {
				//printf("bscfloat-3@@@degenerated ...\n");
				//T t = lineRoot(c.kk0, c.kk2);
				//T fk = evaluateBezier(c.k0, c.k1, c.k2, c.k3, t);
				//T fk = _evaluateBezier(c.k0, c.k1, c.k2, c.k3, c.kk0,-c.kk2);

				floating<T> fk = _evaluateBezier(c.k0, c.k1, c.k2, c.k3, c.kk0, -(c.kk1*T(2.0) - c.kk0));//fk is the value of extreme point

				if ((c.kk0 - c.kk1).sign() == 0){
					return 2;//conservative operation
				}
				else if ((c.kk0 - c.kk1).sign() < 0){
						fk = -fk;
						if (sign_k0 == fk.sign())
							return 0;
						else
							return 2;//conservative operation
				}
			}*/

			floating<T> s0, s1;
			floating<T> t0, t1;

			//degenerated case
			if (!bezierDecomposition(c.k0, c.k1, c.k2, c.k3, c.kk0, c.kk1, c.kk2, s0, s1, t0, t1))//which indicates Y'(t) is a linear polynomial
			{
				//printf("bscfloat-4@@@degenerated ...\n");
				//T t = lineRoot(c.kk0, c.kk2);
				//T fk = evaluateBezier(c.k0, c.k1, c.k2, c.k3, t);
				//T fk = _evaluateBezier(c.k0, c.k1, c.k2, c.k3, c.kk0, -c.kk2);
				floating<T> fk = _evaluateBezier(c.k0, c.k1, c.k2, c.k3, c.kk0, -(c.kk1*floating<T>(2.0, 0) - c.kk0));//fk is the value of extreme point

				if ((c.kk0 - c.kk1).sign() == 0)
					return 2;//conservative operation

				if ((c.kk0 - c.kk1).sign() < 0)//if (c.kk0<0)//.is_certainly_negative())
					fk = -fk;

				if (c.sign_k0 == fk.sign())
					return 0;
				else
					return 2;//conservative operation
			}

			//ensure Y''(0) and Y''(1) have the same sign, which indicate Y'(t) is monotonic
			floating<T> kkk0 = c.k2 - c.k1*floating<T>(2.0, 0) + c.k0;
			floating<T> kkk1 = c.k3 - c.k2*floating<T>(2.0, 0) + c.k1;

			//printf("T-bound: kkk0:%.32f kkk1:%.32f\n", kkk0.sigma, kkk1.sigma);

			int sign_kkk0, sign_kkk1;
			if (kkk0.sign() == 0 && kkk1.sign() != 0){
				sign_kkk0 = kkk1.sign();
				sign_kkk1 = kkk1.sign();
			}
			else if (kkk0.sign() != 0 && kkk1.sign() == 0){
				sign_kkk0 = kkk0.sign();
				sign_kkk1 = kkk0.sign();
			}
			/*assert(sign_kkk0 == sign_kkk1);*/

			//ensure Y'(0) and Y'(1) have the different nonzero signs, which indicate the existing of extreme point
			//int sign_kk0, sign_kk1;
			if (c.kk0.sign() == 0 && c.kk2.sign() == 0){
				if (sign_kkk0 == 0){
					return 2;//conservative operation
				}
				else if (sign_kkk0 == 1){
					c.sign_kk0 = -1;
					c.sign_kk2 = 1;
				}
				else{
					c.sign_kk0 = 1;
					c.sign_kk2 = -1;
				}
			}
			else if (c.kk0.sign() == 0 && c.kk2.sign() != 0){
				c.sign_kk0 = -c.kk2.sign();
				//sign_kk1 = c.kk1.sign();
			}
			else if (c.kk0.sign() != 0 && c.kk2.sign() == 0){
				//sign_kk0 = c.kk0.sign();
				c.sign_kk2 = -c.kk0.sign();
			}
			/*assert(c.sign_kk0 != 0);
			assert(c.sign_kk0 == -c.sign_kk2);*/

			//T(0) == T(1)
			if (t0.sign() == t1.sign()){//(sameSign(t0, t1) == RETURN_TRUE) {
				if (c.sign_k0 == t0.sign()){
					if (flag)
						return 1;//conservative operation
					return 0;
				}
				else
					return 2;//conservative operation if(t0.sign()==0)
			}

			//T(0) != T(1)
			//T t = lineRoot(t0, t1);
			//T fk = evaluateBezier2(c.kk0, c.kk1, c.kk2, t);
			//T fk = _evaluateBezier2(c.k0, c.k1, c.k2, t0, -t1);
			floating<T> fk = _evaluateBezier2(c.kk0, c.kk1, c.kk2, t0, -t1);//fk = Y'(t_T)

			if (fk.sign() == 0)
				return 2;//conservative operation

			if (c.sign_kk0 == fk.sign()) {
				if (c.sign_k0 == t1.sign()){//(sameSign(c.k0, t1))
					if (flag)
						return 1;//conservative operation
					return 0;
				}
				else
					return 2; //conservative operation if(t1.sign()==0)
			}
			else if (c.sign_kk2 = fk.sign()){
				if (c.sign_k0 == t0.sign()){
					if (flag)
						return 1;//conservative operation
					return 0;
				}
				else
					return 2;//conservative operation if(t0.sign()==0)
			}
		}
		throw "TightCCD: coplanarTest failed";
	}


	template <class T>
	inline bool getSimplifyed(
		const floating<T>& k0, const floating<T>& k1, const floating<T>& k2, const floating<T>& k3,
		const floating<T>& l0, const floating<T>& l1, const floating<T>& l2, const floating<T>& l3, const floating<T>& l4,
		floating<T>& j0, floating<T>& j1, floating<T>& j2)
	{
		floating<T> kk0 = k0*floating<T>(4.0, 0);
		floating<T> kk1 = k0 + k1*floating<T>(3.0, 0);
		floating<T> kk2 = (k1 + k2)*floating<T>(2.0, 0);
		floating<T> kk3 = k2*floating<T>(3.0, 0) + k3;
		floating<T> kk4 = k3*floating<T>(4.0, 0);

		floating<T> s0 = (l1*kk0 - l0*kk1)*floating<T>(12.0, 0);
		floating<T> s1 = (l2*kk0 - l0*kk2)*floating<T>(6.0, 0);
		floating<T> s2 = (l3*kk0 - l0*kk3)*floating<T>(4.0, 0);
		floating<T> s3 = (l4*kk0 - l0*kk4)*floating<T>(3.0, 0);

		j0 = (s1*k0 - s0*k1)*floating<T>(6.0, 0);
		j1 = (s2*k0 - s0*k2)*floating<T>(3.0, 0);
		j2 = (s3*k0 - s0*k3)*floating<T>(2.0, 0);

		return true;
	}


	template <class T>
	inline bool getSigns(const floating<T>& t0, const floating<T>& t1, bcrv<T> &c, floating<T> &lt0, floating<T> &lt1, int root_nums)
	{
		/*assert(root_nums > 0);*/

		if (sameSign(t0, t1) == RETURN_TRUE) {
			lt0 = t0;
			lt1 = t0;
			return true;
		}

		if ((t0 - t1).sign() == 0){
			return false;//here have some problem
		}

		//one root of Y(t)=0
		if (root_nums == 1){
			//((c.ct == 0) || (c.ct == 1 && diffSign(c.k0, c.k3)) == RETURN_TRUE) {
			//T t = lineRoot(t0, t1);
			//T ft = evaluateBezier(c.k0, c.k1, c.k2, c.k3, t);
			floating<T> ft = _evaluateBezier(c.k0, c.k1, c.k2, c.k3, t0, -t1);//ft = Y(t_L)

			if (ft.sign() == 0){
				lt0 = floating<T>(0, 1);//So lt0.sign()==0
				lt1 = lt0;//So lt1.sign()==0
				return true;
			}

			if ((t0 - t1).sign() < 0)//if (t0<0)//.is_certainly_negative())
				ft = -ft;

			if (ft.sign() == c.sign_k0){//(sameSign(ft, c.k0) == RETURN_TRUE) {
				lt0 = t1;
				lt1 = t1;
			}
			else {
				lt0 = t0;
				lt1 = t0;
			}
			return true;
		}

		//two roots of Y(t)=0
		if (root_nums == 2) {
			//				T t = lineRoot(t0, t1);
			//				T ft = evaluateBezier(c.k0, c.k1, c.k2, c.k3, t);
			floating<T> ft = _evaluateBezier(c.k0, c.k1, c.k2, c.k3, t0, -t1);//ft = Y(t_L)

			if (ft.sign() == 0){
				floating<T> fk = _evaluateBezier2(c.kk0, c.kk1, c.kk2, t0, -t1);//fk = Y'(t_L)
				if (c.sign_kk0 == 0){
					/*assert(c.sign_kk2 == 0);*/
					lt0 = floating<T>(0, 1);//So lt0.sign()==0
					lt1 = lt0;//So lt1.sign()==0
					return true;
				}
				if (fk.sign() == c.sign_kk0){
					lt0 = floating<T>(0, 1);//So lt0.sign()==0
					lt1 = t1;
					return true;
				}
				if (fk.sign() == c.sign_kk2){
					lt0 = t0;
					lt1 = floating<T>(0, 1);//So lt1.sign()==0
					return true;
				}
				if (fk.sign() == 0){
					lt0 = t0;
					lt1 = t1;
					return true;
				}
			}

			if ((t0 - t1).sign()<0)//if (t0<0)//.is_certainly_negative())
				ft = -ft;

			if (ft.sign() != c.sign_k0){//(diffSign(ft, c.k0)) {
				lt0 = t0;
				lt1 = t1;
				return true;
			}

			//T fk = evaluateBezier2(c.kk0, c.kk1, c.kk2, t);
			floating<T> fk = _evaluateBezier2(c.kk0, c.kk1, c.kk2, t0, -t1);//fk = Y'(t_L)
			if (c.sign_kk0 == 0){
				/*assert(c.sign_kk2 == 0);*/
				lt0 = floating<T>(0, 1);//So lt0.sign()==0
				lt1 = lt0;//So lt1.sign()==0
				return true;
			}

			if (fk.sign() == c.sign_kk0){
				lt0 = t1;
				lt1 = t1;
				return true;
			}
			if (fk.sign() == c.sign_kk2){
				lt0 = t0;
				lt1 = t0;
				return true;
			}
			if (fk.sign() == 0){
				lt0 = t0;
				lt1 = t1;
				//printf("Here is unreasonable!\n");
				return true;
			}
		}

		//printf("Impossible to be here!\n");
		return false;
	}


	/*template <class T>
	inline int insideTest(
		const Vec<3, floating<T>> &a0, const Vec<3, floating<T>> &b0, const Vec<3, floating<T>> &c0, const Vec<3, floating<T>> &d0,
		const Vec<3, floating<T>> &a1, const Vec<3, floating<T>> &b1, const Vec<3, floating<T>> &c1, const Vec<3, floating<T>> &d1,
		bcrv<T> &c, bool ee_test, int root_nums)
	{
		floating<T> l0, l1, l2, l3, l4;
		floating<T> j0, j1, j2;
		floating<T> s0, s1; // for L(t)
		floating<T> t0, t1; // for K(t)

		floating<T> lt0, lt1, kt0, kt1; // for signs of lt and kt

		bool bt0[3], bt1[3];

		for (int i = 0; i<3; i++) {
			bt0[i] = true;
			bt1[i] = true;

			getBezier4(a0, b0, c0, d0, a1, b1, c1, d1, l0, l1, l2, l3, l4, i, ee_test);
			getSimplifyed(c.k0, c.k1, c.k2, c.k3, l0, l1, l2, l3, l4, j0, j1, j2);

			if ((j0 + j2 - j1*floating<T>(2.0, 0)) == 0) {// degenerate j0, j1, j2
				if ((j1 - j0) == 0){
					if (j0.sign() < 0)
						return false;
					else
						continue;
				}

				//printf("@@@Degenerated...\n");
				//getSigns(j0, j2, c, lt0, lt1);//wzd modify 2015.1.28
				getSigns(j0, j1*floating<T>(2.0, 0) - j0, c, lt0, lt1, root_nums);

				if (root_nums == 1){//(c.ct == 0 || (c.ct == 1 && diffSign(c.k0, c.k3)) == RETURN_TRUE) {//(c.ct == 0){////wzd modify 2015.1.28
					if (lt0.sign() < 0)
						return false;
				}
				else if (root_nums == 2){ //(c.ct == 1){
					//if (lt0 < 0)
					if (lt0.sign() < 0)
						bt0[i] = false;//bt0 = false;//wzd modify 2015.1.28

					//if (lt1 < 0)
					if (lt1.sign() < 0)
						bt1[i] = false;//bt1 = false;//wzd modify 2015.1.28

					if (!bt0[i] && !bt1[i])//(!bt0 && !bt1)//wzd modify 2015.1.28
						return false;
				}
				else if (c.ct == 2){
					printf("Impossible to be here!\n");
				}

				continue;
			}

			if (false == bezierDecomposition(c.k0, c.k1, c.k2, c.k3, j0, j1, j2, s0, s1, t0, t1)) {
				getSigns(j0, j1*floating<T>(2.0, 0) - j0, c, lt0, lt1, root_nums);//getSigns(j0, j2, c, lt0, lt1);//wzd modify 2015.1.28

				if (root_nums == 1){//(c.ct == 0 || (c.ct == 1 && diffSign(c.k0, c.k3)) == RETURN_TRUE) {//(c.ct == 0){////wzd modify 2015.1.28
					if (lt0.sign() < 0)
						return false;
				}
				else if (root_nums == 2){ //(c.ct == 1){
					if (lt0.sign() < 0)
						bt0[i] = false;//bt0 = false;//wzd modify 2015.1.28

					if (lt1.sign() < 0)
						bt1[i] = false;//bt1 = false;//wzd modify 2015.1.28

					if (!bt0[i] && !bt1[i])//(!bt0 && !bt1)//wzd modify 2015.1.28
						return false;
				}
				else if (c.ct == 2){
					printf("Impossible to be here!\n");
				}

				continue;
			}

			bool dg_l = getSigns(s0, s1, c, lt0, lt1, root_nums);//for L(t)
			bool dg_k = getSigns(t0, t1, c, kt0, kt1, root_nums);//for K(t)

			if (!dg_l || !dg_k){
				continue;
			}

			if (root_nums == 1){//(c.ct == 0 || (c.ct == 1 && diffSign(c.k0, c.k3)) == RETURN_TRUE) {//(c.ct == 0){////wzd modify 2015.1.28
				if (sameSign(lt0, kt0) == RETURN_TRUE)
					return false;
				continue;
			}

			// kill an possiblity
			if (sameSign(lt0, kt0) == RETURN_TRUE)
				bt0[i] = false;//bt0 = false;//wzd modify 2015.1.28

			// kill an possiblity
			if (sameSign(lt1, kt1) == RETURN_TRUE)
				bt1[i] = false;//bt1 = false;//wzd modify 2015.1.28

			//if no possiblity left, return false ...
			if (!bt0[i] && !bt1[i])//(!bt0 && !bt1)//wzd modify 2015.1.28
				return false;
		}
		//wzd modify 2015.1.28
		if (root_nums == 1)//(c.ct == 0 || (c.ct == 1 && diffSign(c.k0, c.k3)) == RETURN_TRUE)
			return true;
		else if (root_nums == 2){ //(c.ct == 1){
			bool bb0 = (bt0[0] && bt0[1] && bt0[2]);
			bool bb1 = (bt1[0] && bt1[1] && bt1[2]);
			if (!bb0 && !bb1)
				return false;
		}
		//end modify
		return true;
	}*/


	template <class T>
	inline int insideTest(bcrv<T> c, floating<T>* g, int root_nums)
	{
		floating<T> l0, l1, l2, l3, l4;
		floating<T> j0, j1, j2;
		floating<T> s0, s1; // for L(t)
		floating<T> t0, t1; // for K(t)

		floating<T> lt0, lt1, kt0, kt1; // for signs of lt and kt

		bool bt0[3], bt1[3];

		for (int i = 0; i<3; i++) {
			bt0[i] = true;
			bt1[i] = true;

			l0 = g[i * 5 + 0]; l1 = g[i * 5 + 1]; l2 = g[i * 5 + 2];
			l3 = g[i * 5 + 3]; l4 = g[i * 5 + 4];

			getSimplifyed(c.k0, c.k1, c.k2, c.k3, l0, l1, l2, l3, l4, j0, j1, j2);

			/*if ((j0 + j2 - j1*floating<T>(2.0, 0)) == 0) {// degenerate j0, j1, j2
				if ((j1 - j0) == 0){
					if (j0.sign() < 0)
						return false;
					else
						continue;
				}

				//printf("@@@Degenerated...\n");
				//getSigns(j0, j2, c, lt0, lt1);//wzd modify 2015.1.28
				getSigns(j0, j1*floating<T>(2.0, 0) - j0, c, lt0, lt1, root_nums);

				if (root_nums == 1){//(c.ct == 0 || (c.ct == 1 && diffSign(c.k0, c.k3)) == RETURN_TRUE) {//(c.ct == 0){////wzd modify 2015.1.28
					if (lt0.sign()<0)
						return false;
				}
				else if (root_nums == 2){ //(c.ct == 1){
					//if (lt0 < 0)
					if (lt0.sign()<0)
						bt0[i] = false;//bt0 = false;//wzd modify 2015.1.28

					//if (lt1 < 0)
					if (lt1.sign()<0)
						bt1[i] = false;//bt1 = false;//wzd modify 2015.1.28

					if (!bt0[i] && !bt1[i])//(!bt0 && !bt1)//wzd modify 2015.1.28
						return false;
				}
				else if (c.ct == 2){
					printf("Impossible to be here!\n");
				}

				continue;
			}*/

			if (false == bezierDecomposition(c.k0, c.k1, c.k2, c.k3, j0, j1, j2, s0, s1, t0, t1)) {
				if ((j1 - j0) == 0){
					if (j0.sign() < 0)
						return false;
					else
						continue;
				}

				getSigns(j0, j1*floating<T>(2.0, 0) - j0, c, lt0, lt1, root_nums);//getSigns(j0, j2, c, lt0, lt1);//wzd modify 2015.1.28

				if (root_nums == 1){//(c.ct == 0 || (c.ct == 1 && diffSign(c.k0, c.k3)) == RETURN_TRUE) {//(c.ct == 0){////wzd modify 2015.1.28
					if (lt0.sign() < 0)
						return false;
				}
				else if (root_nums == 2){ //(c.ct == 1){
					if (lt0.sign()<0)
						bt0[i] = false;//bt0 = false;//wzd modify 2015.1.28

					if (lt1.sign()<0)
						bt1[i] = false;//bt1 = false;//wzd modify 2015.1.28

					if (!bt0[i] && !bt1[i])//(!bt0 && !bt1)//wzd modify 2015.1.28
						return false;
				}
				else if (c.ct == 2){
					printf("bscwzd-Impossible to be here!\n");
				}

				continue;
			}

			bool dg_l = getSigns(s0, s1, c, lt0, lt1, root_nums);//for L(t)
			bool dg_k = getSigns(t0, t1, c, kt0, kt1, root_nums);//for K(t)

			if (!dg_l || !dg_k){
				continue;
			}

			if (root_nums == 1){//(c.ct == 0 || (c.ct == 1 && diffSign(c.k0, c.k3)) == RETURN_TRUE) {//(c.ct == 0){////wzd modify 2015.1.28
				if (sameSign(lt0, kt0) == RETURN_TRUE)
					return false;
				continue;
			}

			// kill an possiblity
			if (sameSign(lt0, kt0) == RETURN_TRUE)
				bt0[i] = false;//bt0 = false;//wzd modify 2015.1.28

			// kill an possiblity
			if (sameSign(lt1, kt1) == RETURN_TRUE)
				bt1[i] = false;//bt1 = false;//wzd modify 2015.1.28

			//if no possiblity left, return false ...
			if (!bt0[i] && !bt1[i])//(!bt0 && !bt1)//wzd modify 2015.1.28
				return false;
		}
		/*//wzd modify 2015.1.28
		if (root_nums == 1)//(c.ct == 0 || (c.ct == 1 && diffSign(c.k0, c.k3)) == RETURN_TRUE)
			return true;
		else*/
		if (root_nums == 2){ //(c.ct == 1){
			bool bb0 = (bt0[0] && bt0[1] && bt0[2]);
			bool bb1 = (bt1[0] && bt1[1] && bt1[2]);
			if (!bb0 && !bb1)
				return false;
		}
		//end modify
		return true;
	}


	template <class T>
	inline int subdivide(floating<T> k0, floating<T> k1, floating<T> k2, floating<T> k3, floating<T>* g){
		// f'' = 6*(k2-2*k1+k0)*B^0_1 + 6*(k3-2*k2+k1)*B^1_1
		// inflection = f''(0)/(f''(0) - f''(1))
		floating<T> t = k0 - k1*floating<T>(2.0, 0) + k2;
		floating<T> division = k0 - k1*floating<T>(3.0, 0) + k2*floating<T>(3.0, 0) - k3;

		//k0, (1-t)*k0 + t*k1, (1-t)^2*k0 + 2*(1-t)*t*k1 + t^2*k2, (1-t)^3*k0 + 3*(1-t)^2*t*k1 + 3*(1-t)*t^2*k2 + t^3*k3
		floating<T> l0 = k0 * division * division * division;
		floating<T> l1 = ((division - t) * k0 + t * k1) * division * division;
		floating<T> l2 = ((division - t)*(division - t)*k0 + floating<T>(2.0, 0) * t*(division - t)*k1 + t*t*k2)*division;
		floating<T> l3 = (division - t)*(division - t)*(division - t)*k0 + floating<T>(3.0, 0)*t*(division - t)*(division - t)*k1
			+ floating<T>(3.0, 0)*t*t*(division - t)*k2 + t*t*t*k3;

		//(1-t)^3*k0 + 3*(1-t)^2*t*k1 + 3*(1-t)*t^2*k2 + t^3*k3, (1-t)^2*k1 + 2*(1-t)*t*k2 + t^2*k3, (1-t)*k2 + t*k3, k3
		floating<T> r0 = l3;
		floating<T> r1 = ((division - t)*(division - t)*k1 + floating<T>(2.0, 0) * t*(division - t)*k2 + t*t*k3)*division;
		floating<T> r2 = ((division - t) * k2 + t * k3) * division * division;
		floating<T> r3 = k3 * division * division * division;

		if (division.sign() == 0)
			return true; //conservative operation

		if (division.sign() < 0){
			l0 = -l0;
			l1 = -l1;
			l2 = -l2;
			l3 = -l3;
			r0 = -r0;
			r1 = -r1;
			r2 = -r2;
			r3 = -r3;
		}

		floating<T> ll0, ll1, ll2;
		floating<T> rr0, rr1, rr2;

		int ct_l = bezierClassification(l0, l1, l2, l3, ll0, ll1, ll2);
		int ct_r = bezierClassification(r0, r1, r2, r3, rr0, rr1, rr2);

		if (ct_l == 2){
			if (diffSign(ll0, ll2) != RETURN_FALSE)
				ct_l = 1; // no inflexion, 1 extreme
			else
				ct_l = 0; // no inflexion, no extreme
		}
		if (ct_r == 2){
			if (diffSign(rr0, rr2) != RETURN_FALSE)
				ct_r = 1; // no inflexion, 1 extreme
			else
				ct_r = 0; // no inflexion, no extreme
		}

		bcrv<T> c_l(l0, l1, l2, l3, ll0, ll1, ll2, ct_l);
		bcrv<T> c_r(r0, r1, r2, r3, rr0, rr1, rr2, ct_r);

		int left = true;
		int right = true;

		left = coplanarTest(c_l);
		right = coplanarTest(c_r);

		if (!left && !right){
			return false;//return -1;//
		}

		if (left){
			floating<T> gl[15];
			for (int i = 0; i < 3; i++){
				gl[i * 5 + 0] = g[i * 5 + 0] * division*division*division*division;
				gl[i * 5 + 1] = ((division - t)*g[i * 5 + 0]
					+ t*g[i * 5 + 1])*division*division*division;
				gl[i * 5 + 2] = ((division - t)*(division - t)*g[i * 5 + 0]
					+ floating<T>(2.0, 0) * t*(division - t)*g[i * 5 + 1]
					+ t*t*g[i * 5 + 2])*division*division;
				gl[i * 5 + 3] = ((division - t)*(division - t)*(division - t)*g[i * 5 + 0]
					+ floating<T>(3.0, 0) * t*(division - t)*(division - t)*g[i * 5 + 1]
					+ floating<T>(3.0, 0) * t*t*(division - t)*g[i * 5 + 2]
					+ t*t*t*g[i * 5 + 3])*division;
				gl[i * 5 + 4] = ((division - t)*(division - t)*(division - t)*(division - t)*g[i * 5 + 0]
					+ floating<T>(4.0, 0) * t*(division - t)*(division - t)*(division - t)*g[i * 5 + 1]
					+ floating<T>(6.0, 0) * t*t*(division - t)*(division - t)*g[i * 5 + 2]
					+ floating<T>(4.0, 0) * t*t*t*(division - t)*g[i * 5 + 3]
					+ t*t*t*t*g[i * 5 + 4]);
			}
			left = insideTest(c_l, gl, left);
		}

		if (right){
			floating<T> gr[15];
			for (int i = 0; i < 3; i++){
				gr[i * 5 + 4] = g[i * 5 + 4] * division*division*division*division;
				gr[i * 5 + 3] = ((division - t)*g[i * 5 + 3]
					+ t*g[i * 5 + 4])*division*division*division;
				gr[i * 5 + 2] = ((division - t)*(division - t)*g[i * 5 + 2]
					+ floating<T>(2.0, 0) * t*(division - t)*g[i * 5 + 3]
					+ t*t*g[i * 5 + 4])*division*division;
				gr[i * 5 + 1] = ((division - t)*(division - t)*(division - t)*g[i * 5 + 1]
					+ floating<T>(3.0, 0) * t*(division - t)*(division - t)*g[i * 5 + 2]
					+ floating<T>(3.0, 0) * t*t*(division - t)*g[i * 5 + 3]
					+ t*t*t*g[i * 5 + 4])*division;
				gr[i * 5 + 0] = ((division - t)*(division - t)*(division - t)*(division - t)*g[i * 5 + 0]
					+ floating<T>(4.0, 0) * t*(division - t)*(division - t)*(division - t)*g[i * 5 + 1]
					+ floating<T>(6.0, 0) * t*t*(division - t)*(division - t)*g[i * 5 + 2]
					+ floating<T>(4.0, 0) * t*t*t*(division - t)*g[i * 5 + 3]
					+ t*t*t*t*g[i * 5 + 4]);
			}
			right = insideTest(c_r, gr, right);
		}

		if (!left && !right)//
			return false;//return -1;//             //here has false negatives

		return true;
	}


	template <class T>
	inline int Intersect_robust(
		const Vec<3, floating<T>> &a0, const Vec<3, floating<T>> &b0, const Vec<3, floating<T>> &c0, const Vec<3, floating<T>> &d0,
		const Vec<3, floating<T>> &a1, const Vec<3, floating<T>> &b1, const Vec<3, floating<T>> &c1, const Vec<3, floating<T>> &d1,
		bool ee_test)
	{
		floating<T> k0, k1, k2, k3;
		Vec<3, floating<T>> n0, n1, deltaN, nX;

		if (!getBezier(a0, b0, c0, d0, a1, b1, c1, d1, k0, k1, k2, k3, n0, n1, deltaN, nX)) {
			return false;//return -1;//
		}

		//printf("T-bound: k0:%.32f k1:%.32f k2:%.32f k3:%.32f ", k0.sigma, k1.sigma, k2.sigma, k3.sigma);

		floating<T> kk0, kk1, kk2;
		int ct = bezierClassification(k0, k1, k2, k3, kk0, kk1, kk2);

		floating<T> kkk00 = k2 - k1*floating<T>(2.0) + k0;
		floating<T> kkk11 = k3 - k2*floating<T>(2.0) + k1;
		//printf("%.32f %.32f %.32f %.32f %.32f %.32f",
			//k0.sigma, k3.sigma, kk0.sigma, kk2.sigma, kkk00.sigma, kkk11.sigma);

		//printf("kk0:%.32f kk1:%.32f kk2:%.32f ", kk0.sigma, kk1.sigma, kk2.sigma);

		if (ct == 2) {

			/*floating<T> t = inflexion;
			Vec<3, floating<T>> at = a0*(floating<T>(1.0) - t) + a1*t;
			Vec<3, floating<T>> bt = b0*(floating<T>(1.0) - t) + b1*t;
			Vec<3, floating<T>> ct = c0*(floating<T>(1.0) - t) + c1*t;
			Vec<3, floating<T>> dt = d0*(floating<T>(1.0) - t) + d1*t;

			bool ret1 = Intersect_robust(a0, b0, c0, d0, at, bt, ct, dt, ee_test);
			bool ret2 = Intersect_robust(at, bt, ct, dt, a1, b1, c1, d1, ee_test);

			std::cout << "Here!" << std::endl;

			return ret1 || ret2;*/

			floating<T> g[15]; // for inside test
			for (int i = 0; i < 3; i++){
				getBezier4(a0, b0, c0, d0, a1, b1, c1, d1, n0, n1, deltaN, nX, g[i * 5 + 0], g[i * 5 + 1], g[i * 5 + 2], g[i * 5 + 3], g[i * 5 + 4], i, ee_test);
			}

			//std::cout << std::endl;

			return subdivide(k0, k1, k2, k3, g);//here has false negatives

			return true;
		}

		bcrv<T> c(k0, k1, k2, k3, kk0, kk1, kk2, ct);

		int root_nums = coplanarTest(c);

		//std::cout << std::endl;

		if (root_nums == 0)
			return false;//return -1;//

		floating<T> g[15]; // for inside test
		for (int i = 0; i < 3; i++){
			getBezier4(a0, b0, c0, d0, a1, b1, c1, d1, n0, n1, deltaN, nX, g[i * 5 + 0], g[i * 5 + 1], g[i * 5 + 2], g[i * 5 + 3], g[i * 5 + 4], i, ee_test);
		}

		return insideTest(c, g, root_nums);
		//insideTest(a0, b0, c0, d0, a1, b1, c1, d1, c, ee_test, root_nums);

		return true;
	}


	int Intersect_VF_robust(
		const Vec3d &a0, const Vec3d &b0, const Vec3d &c0, const Vec3d &d0,
		const Vec3d &a1, const Vec3d &b1, const Vec3d &c1, const Vec3d &d1)
	{
		//DNF culling with conservative bound
		Vec<3,floating<double>> d_a0, d_b0, d_c0, d_d0, d_a1, d_b1, d_c1, d_d1;
		make_vector(a0, d_a0);
		make_vector(b0, d_b0);
		make_vector(c0, d_c0);
		make_vector(d0, d_d0);
		make_vector(a1, d_a1);
		make_vector(b1, d_b1);
		make_vector(c1, d_c1);
		make_vector(d1, d_d1);

		/*if (!DNF_Culling(d_a0, d_b0, d_c0, d_d0, d_a1, d_b1, d_c1, d_d1))
			return false;*/

		int ret = Intersect_robust(d_a0, d_b0, d_c0, d_d0, d_a1, d_b1, d_c1, d_d1, false);

		return ret;
	}


	int Intersect_EE_robust(
		const Vec3d &a0, const Vec3d &b0, const Vec3d &c0, const Vec3d &d0,
		const Vec3d &a1, const Vec3d &b1, const Vec3d &c1, const Vec3d &d1)
	{
		//DNF culling with conservative bound
		Vec<3, floating<double>> d_a0, d_b0, d_c0, d_d0, d_a1, d_b1, d_c1, d_d1;
		make_vector(a0, d_a0);
		make_vector(b0, d_b0);
		make_vector(c0, d_c0);
		make_vector(d0, d_d0);
		make_vector(a1, d_a1);
		make_vector(b1, d_b1);
		make_vector(c1, d_c1);
		make_vector(d1, d_d1);

		/*if (!DNF_Culling(d_a0, d_b0, d_c0, d_d0, d_a1, d_b1, d_c1, d_d1))
			return false;*/

		int ret = Intersect_robust(d_a0, d_b0, d_c0, d_d0, d_a1, d_b1, d_c1, d_d1, true);

		return ret;
	}
}

#endif
