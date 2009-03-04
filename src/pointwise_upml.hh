#ifndef POINTWISE_UPML_HH_
#define POINTWISE_UPML_HH_

#include "pointwise_dummy.hh"

namespace gmes
{
class UPMLElectric: public DummyElectric
{
public:
	UPMLElectric(const int * const idx, int size, double epsilon_r,
			double c1_in, double c2_in, double c3_in, double c4_in,
			double c5_in, double c6_in) :
		DummyElectric(idx, size, epsilon_r), c1(c1_in), c2(c2_in), c3(c3_in),
				c4(c4_in), c5(c5_in), c6(c6_in), d(0)
	{
	}

protected:
	double c1, c2, c3, c4, c5, c6;
	double d;
};

class UPMLEx: public UPMLElectric
{
public:
	UPMLEx(const int * const idx, int size, double epsilon_r, double c1,
			double c2, double c3, double c4, double c5, double c6) :
		UPMLElectric(idx, size, epsilon_r, c1, c2, c3, c4, c5, c6)
	{
	}

	void update(double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const hz, int hz_x_size, int hz_y_size,
			int hz_z_size, const double * const hy, int hy_x_size,
			int hy_y_size, int hy_z_size, double dt, double dy, double dz);
};

class UPMLEy: public UPMLElectric
{
public:
	UPMLEy(const int * const idx, int size, double epsilon_r, double c1,
			double c2, double c3, double c4, double c5, double c6) :
		UPMLElectric(idx, size, epsilon_r, c1, c2, c3, c4, c5, c6)
	{
	}

	void update(double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const hx, int hx_x_size, int hx_y_size,
			int hx_z_size, const double * const hz, int hz_x_size,
			int hz_y_size, int hz_z_size, double dt, double dz, double dx);
};

class UPMLEz: public UPMLElectric
{
public:
	UPMLEz(const int * const idx, int size, double epsilon_r, double c1,
			double c2, double c3, double c4, double c5, double c6) :
		UPMLElectric(idx, size, epsilon_r, c1, c2, c3, c4, c5, c6)
	{
	}

	void update(double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const hy, int hy_x_size, int hy_y_size,
			int hy_z_size, const double * const hx, int hx_x_size,
			int hx_y_size, int hx_z_size, double dt, double dx, double dy);
};

class UPMLMagnetic: public DummyMagnetic
{
public:
	UPMLMagnetic(const int * const idx, int size, double mu_r, double c1_in,
			double c2_in, double c3_in, double c4_in, double c5_in,
			double c6_in) :
		DummyMagnetic(idx, size, mu_r), c1(c1_in), c2(c2_in), c3(c3_in), c4(
				c4_in), c5(c5_in), c6(c6_in), b(0)
	{
	}

protected:
	double c1, c2, c3, c4, c5, c6;
	double b;
};

class UPMLHx: public UPMLMagnetic
{
public:
	UPMLHx(const int * const idx, int size, double muR, double c1, double c2,
			double c3, double c4, double c5, double c6) :
		UPMLMagnetic(idx, size, muR, c1, c2, c3, c4, c5, c6)
	{
	}

	void update(double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const ez, int ez_x_size, int ez_y_size,
			int ez_z_size, const double * const ey, int ey_x_size,
			int ey_y_size, int ey_z_size, double dt, double dy, double dz);
};

class UPMLHy: public UPMLMagnetic
{
public:
	UPMLHy(const int * const idx, int size, double muR, double c1, double c2,
			double c3, double c4, double c5, double c6) :
		UPMLMagnetic(idx, size, muR, c1, c2, c3, c4, c5, c6)
	{
	}

	void update(double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const ex, int ex_x_size, int ex_y_size,
			int ex_z_size, const double * const ez, int ez_x_size,
			int ez_y_size, int ez_z_size, double dt, double dz, double dx);
};

class UPMLHz: public UPMLMagnetic
{
public:
	UPMLHz(const int * const idx, int size, double muR, double c1, double c2,
			double c3, double c4, double c5, double c6) :
		UPMLMagnetic(idx, size, muR, c1, c2, c3, c4, c5, c6)
	{
	}

	void update(double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const ey, int ey_x_size, int ey_y_size,
			int ey_z_size, const double * const ex, int ex_x_size,
			int ex_y_size, int ex_z_size, double dt, double dx, double dy);
};
}

#endif /*POINTWISE_UPML_HH_*/
