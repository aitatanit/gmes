#ifndef PW_DUMMY_REAL_HH_
#define PW_DUMMY_REAL_HH_

#include "pw_material_real.hh"
#include "constants.hh"

namespace gmes
{
class DummyElectricReal: public MaterialElectricReal
{
public:
	DummyElectricReal(const int * const idx, int size, double epsilon_r) :
		MaterialElectricReal(idx, size), epsilon(epsilon_r * epsilon0)
	{
	}
	double get_epsilon()
	{
		return epsilon;
	}
	void set_epsilon(double epsilon_r)
	{
		epsilon = epsilon_r * epsilon0;
	}

protected:
	double epsilon;
};

class DummyExReal: public DummyElectricReal
{
public:
	DummyExReal(const int * const idx, int size, double epsilon_r = 1) :
		DummyElectricReal(idx, size, epsilon_r)
	{
	}

	void update(double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			double dt, double dy, double dz);
};

class DummyEyReal: public DummyElectricReal
{
public:
	DummyEyReal(const int * const idx, int size, double epsilon_r = 1) :
		DummyElectricReal(idx, size, epsilon_r)
	{
	}

	void update(double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			double dt, double dz, double dx);
};

class DummyEzReal: public DummyElectricReal
{
public:
	DummyEzReal(const int * const idx, int size, double epsilon_r = 1) :
		DummyElectricReal(idx, size, epsilon_r)
	{
	}

	void update(double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			double dt, double dx, double dy);
};

class DummyMagneticReal: public MaterialMagneticReal
{
public:
	DummyMagneticReal(const int * const idx, int size, double mu_r) :
		MaterialMagneticReal(idx, size), mu(mu_r * mu0)
	{
	}
	double get_mu()
	{
		return mu;
	}
	void set_mu(double mu_r)
	{
		mu = mu_r * mu0;
	}

protected:
	double mu;
};

class DummyHxReal: public DummyMagneticReal
{
public:
	DummyHxReal(const int * const idx, int size, double mu_r = 1) :
		DummyMagneticReal(idx, size, mu_r)
	{
	}

	void update(double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			double dt, double dy, double dz);
};

class DummyHyReal: public DummyMagneticReal
{
public:
	DummyHyReal(const int * const idx, int size, double mu_r = 1) :
		DummyMagneticReal(idx, size, mu_r)
	{
	}

	void update(double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			double dt, double dz, double dx);
};

class DummyHzReal: public DummyMagneticReal
{
public:
	DummyHzReal(const int * const idx, int size, double mu_r = 1) :
		DummyMagneticReal(idx, size, mu_r)
	{
	}

	void update(double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			double dt, double dx, double dy);
};
}

#endif /*PW_DUMMY_REAL_HH_*/
