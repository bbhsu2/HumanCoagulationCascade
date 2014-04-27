/*
	Human Coagulation Cascade
	http://github.com/bbhsu2/HumanCoagulationCascade
	C++ Implementation by Bernard Hsu and anonymous
	3/27/2014
*/

#ifndef RungeKutta_h
#define RungeKutta_h

template < size_t N >
inline void rk2(double (&finish)[N], double (&x)[N], double dt, void (&f)(double (&)[N], double (&)[N]))
{
	double k1[N], k2[N], x1[N];
	double halfdt = 0.5 * dt;

	f(k1, x);

	for (size_t i = 0; i < N; ++i)
	{
		x1[i] = x[i] + k1[i] * halfdt;
	}

	f(k2, x1);

	for (size_t i = 0; i < N; ++i)
	{
		finish[i] = x[i] + k2[i] * dt;
	}
}

#endif
