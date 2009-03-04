#ifndef POINTWISE_ZERO_HH_
#define POINTWISE_ZERO_HH_

#include "pointwise_dummy.hh"

namespace gmes
{
class ZeroEx: public DummyEx
{
public:
	ZeroEx(const int * const idx, int size, double epsilon_r = 1) :
		DummyEx(idx, size, epsilon_r)
	{
	}

	void update(double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const hz, int hz_x_size, int hz_y_size,
			int hz_z_size, const double * const hy, int hy_x_size,
			int hy_y_size, int hy_z_size, double dt, double dy, double dz);
};

class ZeroEy: public DummyEy
{
public:
	ZeroEy(const int * const idx, int size, double epsilon_r = 1) :
		DummyEy(idx, size, epsilon_r)
	{
	}

	void update(double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const hx, int hx_x_size, int hx_y_size,
			int hx_z_size, const double * const hz, int hz_x_size,
			int hz_y_size, int hz_z_size, double dt, double dz, double dx);
};

class ZeroEz: public DummyEz
{
public:
	ZeroEz(const int * const idx, int size, double epsilon_r = 1) :
		DummyEz(idx, size, epsilon_r)
	{
	}

	void update(double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const hy, int hy_x_size, int hy_y_size,
			int hy_z_size, const double * const hx, int hx_x_size,
			int hx_y_size, int hx_z_size, double dt, double dx, double dy);
};

class ZeroHx: public DummyHx
{
public:
	ZeroHx(const int * const idx, int size, double mu_r = 1) :
		DummyHx(idx, size, mu_r)
	{
	}

	void update(double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const ez, int ez_x_size, int ez_y_size,
			int ez_z_size, const double * const ey, int ey_x_size,
			int ey_y_size, int ey_z_size, double dt, double dy, double dz);
};

class ZeroHy: public DummyHy
{
public:
	ZeroHy(const int * const idx, int size, double mu_r = 1) :
		DummyHy(idx, size, mu_r)
	{
	}

	void update(double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const ex, int ex_x_size, int ex_y_size,
			int ex_z_size, const double * const ez, int ez_x_size,
			int ez_y_size, int ez_z_size, double dt, double dz, double dx);
};

class ZeroHz: public DummyHz
{
public:
	ZeroHz(const int * const idx, int size, double mu_r = 1) :
		DummyHz(idx, size, mu_r)
	{
	}

	void update(double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const ey, int ey_x_size, int ey_y_size,
			int ey_z_size, const double * const ex, int ex_x_size,
			int ex_y_size, int ex_z_size, double dt, double dx, double dy);
};
}

#endif /*POINTWISE_ZERO_HH_*/
