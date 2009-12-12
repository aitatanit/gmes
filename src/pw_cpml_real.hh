/* This implementation is based on the following article.
 * S. Gedney, “Perfectly Matched Layer Absorbing Boundary Conditions,”
 * Computational Electrodynamics: The Finite-Difference Time-Domain Method,
 * Third Edition, A. Taflove and S.C. Hagness, eds., Artech House Publishers,
 *  2005, pp. 273-328.
 */

#ifndef PW_CPML_REAL_HH_
#define PW_CPML_REAL_HH_

#include "pw_material_real.hh"
#include "constants.hh"

namespace gmes
{
class CpmlElectricReal: public MaterialElectricReal
{
public:
	CpmlElectricReal(const int * const idx, int size, double epsilon_r,
			double b1_in, double b2_in, double c1_in, double c2_in,
			double kappa1_in, double kappa2_in) :
		MaterialElectricReal(idx, size), epsilon(epsilon_r * epsilon0),
		b1(b1_in), b2(b2_in), c1(c1_in), c2(c2_in),
		kappa1(kappa1_in), kappa2(kappa2_in), psi1(0), psi2(0)
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
	double b1, b2;
	double c1, c2;
	double kappa1, kappa2;
	double psi1, psi2;
};

class CpmlExReal: public CpmlElectricReal
{
public:
	CpmlExReal(const int * const idx, int size, double epsilon_r, double by,
			double bz, double cy, double cz, double kappay, double kappaz) :
		CpmlElectricReal(idx, size, epsilon_r, by, bz, cy, cz, kappay, kappaz)
	{
	}

	void update(double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			double dy, double dz, double dt, double t);
};

class CpmlEyReal: public CpmlElectricReal
{
public:
	CpmlEyReal(const int * const idx, int size, double epsilon_r, double bz,
			double bx, double cz, double cx, double kappaz, double kappax) :
		CpmlElectricReal(idx, size, epsilon_r, bz, bx, cz, cx, kappaz, kappax)
	{
	}

	void update(double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			double dz, double dx, double dt, double t);
};

class CpmlEzReal: public CpmlElectricReal
{
public:
	CpmlEzReal(const int * const idx, int size, double epsilon_r, double bx,
			double by, double cx, double cy, double kappax, double kappay) :
		CpmlElectricReal(idx, size, epsilon_r, bx, by, cx, cy, kappax, kappay)
	{
	}

	void update(double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			double dx, double dy, double dt, double t);
};

class CpmlMagneticReal: public MaterialMagneticReal
{
public:
	CpmlMagneticReal(const int * const idx, int size, double mu_r, double b1_in,
			double b2_in, double c1_in, double c2_in, double kappa1_in,
			double kappa2_in) :
		MaterialMagneticReal(idx, size), mu(mu_r * mu0),
		b1(b1_in), b2(b2_in), c1(c1_in), c2(c2_in),
		kappa1(kappa1_in), kappa2(kappa2_in), psi1(0), psi2(0)
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
	double b1, b2;
	double c1, c2;
	double kappa1, kappa2;
	double psi1, psi2;
};

class CpmlHxReal: public CpmlMagneticReal
{
public:
	CpmlHxReal(const int * const idx, int size, double mu_r, double by, double bz,
			double cy, double cz, double kappay, double kappaz) :
		CpmlMagneticReal(idx, size, mu_r, by, bz, cy, cz, kappay, kappaz)
	{
	}

	void update(double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			double dy, double dz, double dt, double t);
};

class CpmlHyReal: public CpmlMagneticReal
{
public:
	CpmlHyReal(const int * const idx, int size, double mu_r, double bz, double bx,
			double cz, double cx, double kappaz, double kappax) :
		CpmlMagneticReal(idx, size, mu_r, bz, bx, cz, cx, kappaz, kappax)
	{
	}

	void update(double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			double dz, double dx, double dt, double t);
};

class CpmlHzReal: public CpmlMagneticReal
{
public:
	CpmlHzReal(const int * const idx, int size, double mu_r, double bx, double by,
			double cx, double cy, double kappax, double kappay) :
		CpmlMagneticReal(idx, size, mu_r, bx, by, cx, cy, kappax, kappay)
	{
	}

	void update(double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			double dx, double dy, double dt, double t);
};
}

#endif /*PW_CPML_REAL_HH_*/
