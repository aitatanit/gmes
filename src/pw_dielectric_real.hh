#ifndef PW_DIELECTRIC_REAL_HH_
#define PW_DIELECTRIC_REAL_HH_

#include "pw_dummy_real.hh"

namespace gmes
{
class DielectricExReal: public DummyExReal
{
public:
	DielectricExReal(const int * const idx, int size, double epsilon_r = 1) :
		DummyExReal(idx, size, epsilon_r)
	{
	}

	void update(double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			double dt, double dy, double dz);
};

class DielectricEyReal: public DummyEyReal
{
public:
	DielectricEyReal(const int * const idx, int size, double epsilon_r = 1) :
		DummyEyReal(idx, size, epsilon_r)
	{
	}

	void update(double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			double dt, double dz, double dx);
};

class DielectricEzReal: public DummyEzReal
{
public:
	DielectricEzReal(const int * const idx, int size, double epsilon_r = 1) :
		DummyEzReal(idx, size, epsilon_r)
	{
	}

	void update(double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			double dt, double dx, double dy);
};

class DielectricHxReal: public DummyHxReal
{
public:
	DielectricHxReal(const int * const idx, int size, double mu_r = 1) :
		DummyHxReal(idx, size, mu_r)
	{
	}

	void update(double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			double dt, double dy, double dz);
};

class DielectricHyReal: public DummyHyReal
{
public:
	DielectricHyReal(const int * const idx, int size, double mu_r = 1) :
		DummyHyReal(idx, size, mu_r)
	{
	}

	void update(double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			double dt, double dz, double dx);
};

class DielectricHzReal: public DummyHzReal
{
public:
	DielectricHzReal(const int * const idx, int size, double mu_r = 1) :
		DummyHzReal(idx, size, mu_r)
	{
	}

	void update(double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			double dt, double dx, double dy);
};
}

#endif /*PW_DIELECTRIC_REAL_HH_*/
