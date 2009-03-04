#ifndef POINTWISE_ONE_HH_
#define POINTWISE_ONE_HH_

#include "pointwise_dummy.hh"

namespace gmes
{
class OneEx: public DummyEx
{
public:
	OneEx(const int * const idx, int size, double epsilon_r = 1) :
		DummyEx(idx, size, epsilon_r)
	{
	}

	void update(double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const hz, int hz_x_size, int hz_y_size,
			int hz_z_size, const double * const hy, int hy_x_size,
			int hy_y_size, int hy_z_size, double dt, double dy, double dz);
};

class OneEy: public DummyEy
{
public:
	OneEy(const int * const idx, int size, double epsilon_r = 1) :
		DummyEy(idx, size, epsilon_r)
	{
	}

	void update(double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const hx, int hx_x_size, int hx_y_size,
			int hx_z_size, const double * const hz, int hz_x_size,
			int hz_y_size, int hz_z_size, double dt, double dz, double dx);
};

class OneEz: public DummyEz
{
public:
	OneEz(const int * const idx, int size, double epsilon_r = 1) :
		DummyEz(idx, size, epsilon_r)
	{
	}

	void update(double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const hy, int hy_x_size, int hy_y_size,
			int hy_z_size, const double * const hx, int hx_x_size,
			int hx_y_size, int hx_z_size, double dt, double dx, double dy);
};

class OneHx: public DummyHx
{
public:
	OneHx(const int * const idx, int size, double mu_r = 1) :
		DummyHx(idx, size, mu_r)
	{
	}

	void update(double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const ez, int ez_x_size, int ez_y_size,
			int ez_z_size, const double * const ey, int ey_x_size,
			int ey_y_size, int ey_z_size, double dt, double dy, double dz);
};

class OneHy: public DummyHy
{
public:
	OneHy(const int * const idx, int size, double mu_r = 1) :
		DummyHy(idx, size, mu_r)
	{
	}

	void update(double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const ex, int ex_x_size, int ex_y_size,
			int ex_z_size, const double * const ez, int ez_x_size,
			int ez_y_size, int ez_z_size, double dt, double dz, double dx);
};

class OneHz: public DummyHz
{
public:
	OneHz(const int * const idx, int size, double mu_r = 1) :
		DummyHz(idx, size, mu_r)
	{
	}

	void update(double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const ey, int ey_x_size, int ey_y_size,
			int ey_z_size, const double * const ex, int ex_x_size,
			int ex_y_size, int ex_z_size, double dt, double dx, double dy);
};
}

#endif /*POINTWISE_ONE_HH_*/
