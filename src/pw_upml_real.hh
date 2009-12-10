/* This implementation is based on the following article.
 * S. Gedney, “Perfectly Matched Layer Absorbing Boundary Conditions,”
 * Computational Electrodynamics: The Finite-Difference Time-Domain Method,
 * Third Edition, A. Taflove and S.C. Hagness, eds., Artech House Publishers,
 *  2005, pp. 273-328.
 */

#ifndef PW_UPML_REAL_HH_
#define PW_UPML_REAL_HH_

#include "pw_dummy_real.hh"

namespace gmes
{
class UpmlElectricReal: public DummyElectricReal
{
public:
	UpmlElectricReal(const int * const idx, int size, double epsilon_r,
			double c1_in, double c2_in, double c3_in, double c4_in,
			double c5_in, double c6_in) :
		DummyElectricReal(idx, size, epsilon_r), c1(c1_in), c2(c2_in), c3(c3_in),
				c4(c4_in), c5(c5_in), c6(c6_in), d(0)
	{
	}

protected:
	double c1, c2, c3, c4, c5, c6;
	double d;
};

class UpmlExReal: public UpmlElectricReal
{
public:
	UpmlExReal(const int * const idx, int size, double epsilon_r, double c1,
			double c2, double c3, double c4, double c5, double c6) :
		UpmlElectricReal(idx, size, epsilon_r, c1, c2, c3, c4, c5, c6)
	{
	}

	void update(double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			double dy, double dz, double dt, double t);
};

class UpmlEyReal: public UpmlElectricReal
{
public:
	UpmlEyReal(const int * const idx, int size, double epsilon_r, double c1,
			double c2, double c3, double c4, double c5, double c6) :
		UpmlElectricReal(idx, size, epsilon_r, c1, c2, c3, c4, c5, c6)
	{
	}

	void update(double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			double dz, double dx, double dt, double t);
};

class UpmlEzReal: public UpmlElectricReal
{
public:
	UpmlEzReal(const int * const idx, int size, double epsilon_r, double c1,
			double c2, double c3, double c4, double c5, double c6) :
		UpmlElectricReal(idx, size, epsilon_r, c1, c2, c3, c4, c5, c6)
	{
	}

	void update(double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			double dx, double dy, double dt, double t);
};

class UpmlMagneticReal: public DummyMagneticReal
{
public:
	UpmlMagneticReal(const int * const idx, int size, double mu_r, double c1_in,
			double c2_in, double c3_in, double c4_in, double c5_in,
			double c6_in) :
		DummyMagneticReal(idx, size, mu_r), c1(c1_in), c2(c2_in), c3(c3_in), c4(
				c4_in), c5(c5_in), c6(c6_in), b(0)
	{
	}

protected:
	double c1, c2, c3, c4, c5, c6;
	double b;
};

class UpmlHxReal: public UpmlMagneticReal
{
public:
	UpmlHxReal(const int * const idx, int size, double muR, double c1, double c2,
			double c3, double c4, double c5, double c6) :
		UpmlMagneticReal(idx, size, muR, c1, c2, c3, c4, c5, c6)
	{
	}

	void update(double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			double dy, double dz, double dt, double t);
};

class UpmlHyReal: public UpmlMagneticReal
{
public:
	UpmlHyReal(const int * const idx, int size, double muR, double c1, double c2,
			double c3, double c4, double c5, double c6) :
		UpmlMagneticReal(idx, size, muR, c1, c2, c3, c4, c5, c6)
	{
	}

	void update(double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			double dz, double dx, double dt, double t);
};

class UpmlHzReal: public UpmlMagneticReal
{
public:
	UpmlHzReal(const int * const idx, int size, double muR, double c1, double c2,
			double c3, double c4, double c5, double c6) :
		UpmlMagneticReal(idx, size, muR, c1, c2, c3, c4, c5, c6)
	{
	}

	void update(double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			double dx, double dy, double dt, double t);
};
}

#endif /*PW_UPML_REAL_HH_*/
