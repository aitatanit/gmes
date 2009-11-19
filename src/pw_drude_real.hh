/* This implementation is based on the following article.
 * M. Okoniewski and E. Okoniewska, "Drude dispersion in ADE FDTD revisited,"
 * Electron. Lett., vol. 42, no. 9, pp. 503-504, 2006. */

#ifndef PW_DRUDE_REAL_HH_
#define PW_DRUDE_REAL_HH_

#include <vector>

#include "pw_dummy_real.hh"
#include "pw_dielectric_real.hh"

namespace gmes
{
class DrudeElectricReal: public DummyElectricReal
{
public:
	DrudeElectricReal(const int * const idx, int size, double epsilon_inf,
			const double * const omega_p, int omega_p_size,
			const double * const gamma_p, int gamma_p_size) :
		DummyElectricReal(idx, size, epsilon_inf), omega_p(omega_p, omega_p
				+ omega_p_size), gamma_p(gamma_p, gamma_p + gamma_p_size),
				q_new(omega_p_size, 0.), q_old(omega_p_size, 0.)
	{
	}

protected:
	std::vector<double> omega_p;
	std::vector<double> gamma_p;
	std::vector<double> q_new;
	std::vector<double> q_old;
};

class DrudeExReal: public DrudeElectricReal
{
public:
	DrudeExReal(const int * const idx, int size, double epsilon_inf,
			const double * const omega_p, int omega_p_size,
			const double * const gamma_p, int gamma_p_size) :
		DrudeElectricReal(idx, size, epsilon_inf, omega_p, omega_p_size, gamma_p,
				gamma_p_size)
	{
	}

	void update(double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			double dt, double dy, double dz);
};

class DrudeEyReal: public DrudeElectricReal
{
public:
	DrudeEyReal(const int * const idx, int size, double epsilon_inf,
			const double * const omega_p, int omega_p_size,
			const double * const gamma_p, int gamma_p_size) :
		DrudeElectricReal(idx, size, epsilon_inf, omega_p, omega_p_size, gamma_p,
				gamma_p_size)
	{
	}

	void update(double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			double dt, double dz, double dx);
};

class DrudeEzReal: public DrudeElectricReal
{
public:
	DrudeEzReal(const int * const idx, int size, double epsilon_inf,
			const double * const omega_p, int omega_p_size,
			const double * const gamma_p, int gamma_p_size) :
		DrudeElectricReal(idx, size, epsilon_inf, omega_p, omega_p_size, gamma_p,
				gamma_p_size)
	{
	}

	void update(double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			double dt, double dx, double dy);
};

class DrudeHxReal: public DielectricHxReal
{
public:
	DrudeHxReal(const int * const idx, int size, double mu_r = 1) :
		DielectricHxReal(idx, size, mu_r)
	{
	}
};

class DrudeHyReal: public DielectricHyReal
{
public:
	DrudeHyReal(const int * const idx, int size, double mu_r = 1) :
		DielectricHyReal(idx, size, mu_r)
	{
	}
};

class DrudeHzReal: public DielectricHzReal
{
public:
	DrudeHzReal(const int * const idx, int size, double mu_r = 1) :
		DielectricHzReal(idx, size, mu_r)
	{
	}
};

}

#endif /*PW_DRUDE_REAL_HH_*/
