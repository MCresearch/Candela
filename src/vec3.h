#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

template <class T>
class Vector3
{
public:
	T x;
	T y;
	T z;
	Vector3(const T &x1 = 0,const T &y1 = 0,const T &z1 = 0);

	Vector3<T>& operator=(const Vector3<T> &u);
	Vector3<T>& operator+=(const Vector3<T> &u);
	Vector3<T>& operator-=(const Vector3<T> &u);
	Vector3<T>& operator*=(const Vector3<T> &u);
	Vector3<T>& operator*=(const T &s);
	Vector3<T>& operator/=(const Vector3<T> &u);
	Vector3<T>& operator/=(const T &s);
	//qianrui add
    Vector3<T> cos(Vector3<T> v);
    Vector3<T> sin(Vector3<T> v);

	void set(const T &x1, const T &y1,const T &z1);

	inline T norm(void) const
	{ return sqrt(x*x + y*y + z*z); }

	T norm2(void)const;
	void normalize(void);
	void reverse(void);

	// mohan add 2009-11-29
	void print(void)const ;
};

template <class T>
double vector_cos(Vector3<T> &u1, Vector3<T> &v1)
{
	double dot_mult = u1.x*v1.x + u1.y*v1.y + u1.z*v1.z;
	//cout << "dot_mult = " << dot_mult << endl;
	double vec_mag_u = sqrt(u1.x*u1.x + u1.y*u1.y + u1.z*u1.z);
	double vec_mag_v = sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);
	double acos_uv = acos(dot_mult / vec_mag_u / vec_mag_v) * 180 / 3.1415926535897932384626;
	return acos_uv;
}

template <class T>
Vector3<T> operator+(const Vector3<T> &u,const Vector3<T> &v);
template <class T>
Vector3<T> operator-(const Vector3<T> &u,const Vector3<T> &v);
// | i  j  k  |
// | ux uy uz |
// | vx vy vz |
//u.v=(uy*vz-uz*vy)i+(-ux*vz+uz*vx)j+(ux*vy-uy*vx)k
template <class T>
Vector3<T> operator^(const Vector3<T> &u,const Vector3<T> &v);
//u.v=(ux*vx)+(uy*vy)+(uz*vz)
template <class T>
T operator*(const Vector3<T> &u,const Vector3<T> &v);
template <class T>
Vector3<T> operator*(const T &s,const Vector3<T> &u);
template <class T>
Vector3<T> operator/(const Vector3<T> &u,const T &s);
template <class T>
Vector3<T> cross(const Vector3<T> &u,const Vector3<T> &v);
//u.v=(ux*vx)+(uy*vy)+(uz*vz)
template <class T>
T dot(const Vector3<T> &u,const Vector3<T> &v);
//template <class T>
//T TripleSalarProduct(Vector3<T> u, Vector3<T> v, Vector3<T> w);

//==========================================================
// Define Constructor
//==========================================================
template <class T>
Vector3<T>::Vector3(const T &x1,const T &y1, const T &z1)
:x(x1),y(y1),z(z1)
{}

template <class T>
void Vector3<T>::set(const T &x1,const T &y1,const T &z1)
{
	x = x1; y = y1; z = z1;
}

template <class T>
T Vector3<T>::norm2(void)const
{
	return x*x + y*y + z*z;
}

//|v| = sqrt(x*x+y*y+z*z)=1
template <class T>
void Vector3<T>::normalize(void)
{
	T m = sqrt(x*x + y*y + z*z);
	x /= m;y /= m;z /= m;
}

template <class T>
void Vector3<T>::reverse(void)
{
	x = -x;y = -y;z = -z;
}

template <class T>
Vector3<T>& Vector3<T>::operator=(const Vector3<T> &u)
{
	x = u.x;y = u.y;z = u.z;
	return *this;
}

template <class T>
Vector3<T>& Vector3<T>::operator+=(const Vector3<T> &u)
{
	x += u.x;y += u.y;z += u.z;
	return *this;
}

template <class T>
Vector3<T>& Vector3<T>::operator-=(const Vector3<T> &u)
{
	x -= u.x;y -= u.y;z -= u.z;
	return *this;
}

template <class T>
Vector3<T>& Vector3<T>::operator*=(const T &s)
{
	x *= s;y *= s;z *= s;
	return *this;
}

template <class T>
Vector3<T>& Vector3<T>::operator/=(const T &s)
{
	x /= s;y /= s;z /= s;
	return *this;
}

template <class T>
Vector3<T> operator+(const Vector3<T> &u,const Vector3<T> &v)
{
	return Vector3<T>(u.x + v.x, u.y + v.y, u.z + v.z);
}

template <class T>
Vector3<T> operator-(const Vector3<T> &u,const Vector3<T> &v)
{
	return Vector3<T>(u.x - v.x, u.y - v.y, u.z - v.z);
}

//	u.v=(uy*vz-uz*vy)i+(-ux*vz+uz*vx)j+(ux*vy-uy*vzx)k
template <class T>
Vector3<T> operator^(const Vector3<T> &u,const Vector3<T> &v)
{
	return Vector3<T> (u.y * v.z - u.z * v.y,
	                   -u.x * v.z + u.z * v.x,
	                   u.x * v.y - u.y * v.x);
}

//u.v=(ux*vx)+(uy*vy)+(uz*vz)
template <class T>
T operator*(const Vector3<T> &u,const Vector3<T> &v)
{
	return (u.x * v.x + u.y * v.y + u.z * v.z);
}

template <class T>// mohan add 2009-5-10
Vector3<T> operator*(const Vector3<T> &u, const T &s)
{
	return Vector3<T>(u.x * s, u.y * s, u.z * s);
}

template <class T>
Vector3<T> operator*(const T &s,const Vector3<T> &u)
{
	return Vector3<T>(u.x * s, u.y * s, u.z * s);
}

template <class T>
Vector3<T> operator /(const Vector3<T> &u,const T &s)
{
	return Vector3<T>(u.x / s, u.y / s, u.z / s);
}

//	u.v=(uy*vz-uz*vy)i+(-ux*vz+uz*vx)j+(ux*vy-uy*vzx)k
template <class T>
Vector3<T> cross(const Vector3<T> &u, const Vector3<T> &v)
{
	return Vector3<T> (u.y * v.z - u.z * v.y,
	                   -u.x * v.z + u.z * v.x,
	                   u.x * v.y - u.y * v.x);
}

template <class T>
T dot(const Vector3<T> &u,const Vector3<T> &v)
{
	return (u.x * v.x + u.y * v.y + u.z * v.z);
}

//whether u == v
template <class T>
bool operator ==(const Vector3<T> &u, const Vector3<T> &v)
{
	if(u.x == v.x && u.y == v.y && u.z == v.z)
	{
		return true;
	}
	return false;
}

//whether m1 != m2
template <class T>
bool operator !=(const Vector3<T> &u, const Vector3<T> &v)
{
    return !(u == v);
}

template <class T>
void Vector3<T>::print(void)const
{
	cout.precision(5) ;
	cout << "(" << setw(10) << x << "," << setw(10) << y << ","
	<< setw(10) << z  << ")"  << endl ;
	return ;
}
template<class T>//qianrui add
Vector3<T> cos(Vector3<T> v)
{
        Vector3<T> y;
        y.x=cos(v.x);
        y.y=cos(v.y);
        y.z=cos(v.z);
        return y;
}
template<class T>//qianrui add
Vector3<T> sin(Vector3<T> v)
{
        Vector3<T> y;
        y.x=sin(v.x);
        y.y=sin(v.y);
        y.z=sin(v.z);
        return y;
}

//s = u.(v x w)
//template <class T>
//T TripleScalarProduct(Vector3<T> u, Vector3<T> v, Vector3<T> w)
//{
//	return T((u.x * (v.y * w.z - v.z * w.y)) +
//	         (u.y * (-v.x * w.z + v.z * w.x)) +
//	         (u.z * (v.x * w.y - v.y * w.x)));
//}

#endif
