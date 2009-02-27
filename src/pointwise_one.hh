#ifndef POINTWISE_ONE_HH_
#define POINTWISE_ONE_HH_

#include "pointwise_material.hh"

namespace gmes
{
class OneElectric: public PointwiseMaterial
{
public:
	OneElectric(const int * const idx, int size, double epsilon_r);
	double get_epsilon() { return epsilon; }
	void set_epsilon(double epsilon_r);

protected:
	double epsilon;
};


class OneEx: public OneElectric
{
public:
	OneEx(const int * const idx, int size, double epsilon_r=1)
	: OneElectric(idx, size, epsilon_r) {}

	void update(double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			double dt, double dy, double dz);
};


class OneEy: public OneElectric
{
public:
	OneEy(const int * const idx, int size, double epsilon_r=1)
		: OneElectric(idx, size, epsilon_r) {}

	void update(double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			double dt, double dz, double dx);
};


class OneEz: public OneElectric
{
public:
	OneEz(const int * const idx, int size, double epsilon_r=1)
		: OneElectric(idx, size, epsilon_r) {}

	void update(double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			double dt, double dx, double dy);
};


class OneMagnetic: public PointwiseMaterial
{
public:
	OneMagnetic(const int * const idx, int size, double mu_r);
	double get_mu() { return mu; }
	void set_mu(double mu_r);

protected:
	double mu;
};


class OneHx: public OneMagnetic
{
public:
	OneHx(const int * const idx, int size, double mu_r=1)
		: OneMagnetic(idx, size, mu_r) {}

	void update(double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			double dt, double dy, double dz);
};


class OneHy: public OneMagnetic
{
public:
	OneHy(const int * const idx, int size, double mu_r=1)
		: OneMagnetic(idx, size, mu_r) {}

	void update(double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			double dt, double dz, double dx);
};


class OneHz: public OneMagnetic
{
public:
	OneHz(const int * const idx, int size, double mu_r=1)
		: OneMagnetic(idx, size, mu_r) {}

	void update(double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			double dt, double dx, double dy);
};
}

#endif /*POINTWISE_ONE_HH_*/
