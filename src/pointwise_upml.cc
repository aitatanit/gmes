#include "pointwise_upml.hh"
#include "constants.hh"

#define ex(i,j,k) ex[((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[((i)*hz_y_size+(j))*hz_z_size+(k)]

using namespace gmes;

void UPMLEx::update(double * const ex, int ex_x_size, int ex_y_size,
		int ex_z_size, const double * const hz, int hz_x_size, int hz_y_size,
		int hz_z_size, const double * const hy, int hy_x_size, int hy_y_size,
		int hy_z_size, double dt, double dy, double dz)
{
	const double dstore = d;

	d = c1 * d + c2 * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - (hy(i+1,j,k+1)
			- hy(i+1,j,k)) / dz);
	ex(i,j,k) = c3 * ex(i,j,k) + c4 * (c5 * d - c6 * dstore) / epsilon;
}

void UPMLEy::update(double * const ey, int ey_x_size, int ey_y_size,
		int ey_z_size, const double * const hx, int hx_x_size, int hx_y_size,
		int hx_z_size, const double * const hz, int hz_x_size, int hz_y_size,
		int hz_z_size, double dt, double dz, double dx)
{
	const double dstore = d;

	d = c1 * d + c2 * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - (hz(i+1,j+1,k)
			- hz(i,j+1,k)) / dx);
	ey(i,j,k) = c3 * ey(i,j,k) + c4 * (c5 * d - c6 * dstore) / epsilon;
}

void UPMLEz::update(double * const ez, int ez_x_size, int ez_y_size,
		int ez_z_size, const double * const hy, int hy_x_size, int hy_y_size,
		int hy_z_size, const double * const hx, int hx_x_size, int hx_y_size,
		int hx_z_size, double dt, double dx, double dy)
{
	const double dstore = d;

	d = c1 * d + c2 * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - (hx(i,j+1,k+1)
			- hx(i,j,k+1)) / dy);
	ez(i,j,k) = c3 * ez(i,j,k) + c4 * (c5 * d - c6 * dstore) / epsilon;
}

void UPMLHx::update(double * const hx, int hx_x_size, int hx_y_size,
		int hx_z_size, const double * const ez, int ez_x_size, int ez_y_size,
		int ez_z_size, const double * const ey, int ey_x_size, int ey_y_size,
		int ey_z_size, double dt, double dy, double dz)
{
	const double bstore = b;

	b = c1 * b - c2 * ((ez(i,j,k-1) - ez(i,j-1,k-1)) / dy - (ey(i,j-1,k)
			- ey(i,j-1,k-1)) / dz);
	hx(i,j,k) = c3 * hx(i,j,k) + c4 * (c5 * b - c6 * bstore) / mu;
}

void UPMLHy::update(double * const hy, int hy_x_size, int hy_y_size,
		int hy_z_size, const double * const ex, int ex_x_size, int ex_y_size,
		int ex_z_size, const double * const ez, int ez_x_size, int ez_y_size,
		int ez_z_size, double dt, double dz, double dx)
{
	const double bstore = b;

	b = c1 * b - c2 * ((ex(i-1,j,k) - ex(i-1,j,k-1)) / dz - (ez(i,j,k-1)
			- ez(i-1,j,k-1)) / dx);
	hy(i,j,k) = c3 * hy(i,j,k) + c4 * (c5 * b - c6 * bstore) / mu;
}

void UPMLHz::update(double * const hz, int hz_x_size, int hz_y_size,
		int hz_z_size, const double * const ey, int ey_x_size, int ey_y_size,
		int ey_z_size, const double * const ex, int ex_x_size, int ex_y_size,
		int ex_z_size, double dt, double dx, double dy)
{
	const double bstore = b;

	b = c1 * b - c2 * ((ey(i,j-1,k) - ey(i-1,j-1,k)) / dx - (ex(i-1,j,k)
			- ex(i-1,j-1,k)) / dy);
	hz(i,j,k) = c3 * hz(i,j,k) + c4 * (c5 * b - c6 * bstore) / mu;
}
