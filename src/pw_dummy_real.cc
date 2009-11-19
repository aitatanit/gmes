#include "pw_dummy_real.hh"
#include "constants.hh"

#define ex(i,j,k) ex[((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[((i)*hz_y_size+(j))*hz_z_size+(k)]

using namespace gmes;

void DummyExReal::update(double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
		const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
		const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
		double dt, double dy, double dz)
{
	return;
}

void DummyEyReal::update(double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
		const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
		const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
		double dt, double dz, double dx)
{
	return;
}

void DummyEzReal::update(double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
		const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
		const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
		double dt, double dx, double dy)
{
	return;
}

void DummyHxReal::update(double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
		const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
		const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
		double dt, double dy, double dz)
{
	return;
}

void DummyHyReal::update(double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
		const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
		const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
		double dt, double dz, double dx)
{
	return;
}

void DummyHzReal::update(double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
		const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
		const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
		double dt, double dx, double dy)
{
	return;
}
