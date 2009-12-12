#ifndef PW_DIELECTRIC_REAL_HH_
#define PW_DIELECTRIC_REAL_HH_

#include "pw_material_real.hh"
#include "constants.hh"

namespace gmes
{
class DielectricElectricReal: public MaterialElectricReal
{
public:
	DielectricElectricReal(const int * const idx, int size, double epsilon_r = 1) :
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

class DielectricExReal: public DielectricElectricReal
{
public:
	DielectricExReal(const int * const idx, int size, double epsilon_r = 1) :
			DielectricElectricReal(idx, size, epsilon_r)
	{
	}

	void update(double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			double dy, double dz, double dt, double t);
};

class DielectricEyReal: public DielectricElectricReal
{
public:
	DielectricEyReal(const int * const idx, int size, double epsilon_r = 1) :
			DielectricElectricReal(idx, size, epsilon_r)
	{
	}

	void update(double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			double dz, double dx, double dt, double t);
};

class DielectricEzReal: public DielectricElectricReal
{
public:
	DielectricEzReal(const int * const idx, int size, double epsilon_r = 1) :
			DielectricElectricReal(idx, size, epsilon_r)
	{
	}

	void update(double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			double dx, double dy, double dt, double t);
};

class DielectricMagneticReal: public MaterialMagneticReal
{
public:
	DielectricMagneticReal(const int * const idx, int size, double mu_r = 1) :
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

class DielectricHxReal: public DielectricMagneticReal
{
public:
	DielectricHxReal(const int * const idx, int size, double mu_r = 1) :
			DielectricMagneticReal(idx, size, mu_r)
	{
	}

	void update(double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			double dy, double dz, double dt, double t);
};

class DielectricHyReal: public DielectricMagneticReal
{
public:
	DielectricHyReal(const int * const idx, int size, double mu_r = 1) :
			DielectricMagneticReal(idx, size, mu_r)
	{
	}

	void update(double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			double dz, double dx, double dt, double t);
};

class DielectricHzReal: public DielectricMagneticReal
{
public:
	DielectricHzReal(const int * const idx, int size, double mu_r = 1) :
			DielectricMagneticReal(idx, size, mu_r)
	{
	}

	void update(double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			double dx, double dy, double dt, double t);
};
}

#endif /*PW_DIELECTRIC_REAL_HH_*/
