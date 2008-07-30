#ifndef POINTWISE_DIELECTRIC_HH_
#define POINTWISE_DIELECTRIC_HH_

#include "pointwise_material.hh"

namespace gmes 
{
class DielectricElectric: public PointwiseMaterial
{
public:
	DielectricElectric(const int * const idx, int size, double epsilon_r);
	double get_epsilon() { return epsilon; }
	void set_epsilon(double epsilon_r);

protected:
	double epsilon;
};


class DielectricEx: public DielectricElectric
{
public:
	DielectricEx(const int * const idx, int size, double epsilon_r=1)
	: DielectricElectric(idx, size, epsilon_r) {}

	void update(double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			double dt, double dy, double dz);
};


class DielectricEy: public DielectricElectric
{
public:
	DielectricEy(const int * const idx, int size, double epsilon_r=1)
	: DielectricElectric(idx, size, epsilon_r) {}

	void update(double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			double dt, double dz, double dx);
};


class DielectricEz: public DielectricElectric
{
public:
	DielectricEz(const int * const idx, int size, double epsilon_r=1)
	: DielectricElectric(idx, size, epsilon_r) {}

	void update(double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			double dt, double dx, double dy);
};


class DielectricMagnetic: public PointwiseMaterial
{
public:
	DielectricMagnetic(const int * const idx, int size, double mu_r);
	double get_mu() { return mu; }
	void set_mu(double mu_r);

protected:
	double mu;
};


class DielectricHx: public DielectricMagnetic
{
public:
	DielectricHx(const int * const idx, int size, double mu_r=1)
	: DielectricMagnetic(idx, size, mu_r) {}

	void update(double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,	    
			double dt, double dy, double dz);
};


class DielectricHy: public DielectricMagnetic
{
public:
	DielectricHy(const int * const idx, int size, double mu_r=1)
	: DielectricMagnetic(idx, size, mu_r) {}

	void update(double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			double dt, double dz, double dx);
};


class DielectricHz: public DielectricMagnetic
{
public:
	DielectricHz(const int * const idx, int size, double mu_r=1)
	: DielectricMagnetic(idx, size, mu_r) {}

	void update(double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			double dt, double dx, double dy);
};
}

#endif /*POINTWISE_DIELECTRIC_HH_*/
