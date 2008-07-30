#include "pointwise_cpml.hh"
#include "constants.hh"

#define ex(i,j,k) ex[((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[((i)*hz_y_size+(j))*hz_z_size+(k)]

using namespace gmes;

CPMLElectric::CPMLElectric(const int * const idx, int size, double epsilon_r, 
			   double b1_in, double b2_in, 
			   double c1_in, double c2_in, 
			   double kappa1_in, double kappa2_in)
    : PointwiseMaterial(idx, size), 
      b1(b1_in), b2(b2_in), c1(c1_in), c2(c2_in), kappa1(kappa1_in), kappa2(kappa2_in), 
      psi1(0), psi2(0)
{
    epsilon = epsilon_r * epsilon0;
}

void 
CPMLElectric::set_epsilon(double epsilon_r)
{
    epsilon = epsilon_r * epsilon0;
}

void 
CPMLEx::update(double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       double dt, double dy, double dz)
{
    psi1 = b1 * psi1 + c1 * (hz(i+1,j+1,k) - hz(i+1,j,k)) / dy;
    psi2 = b2 * psi2 + c2 * (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz;

    ex(i,j,k) += dt / epsilon *
	((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy / kappa1 - 
	 (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz / kappa2 + psi1 - psi2);
}
         
void
CPMLEy::update(double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       double dt, double dz, double dx)
{
    psi1 = b1 * psi1 + c1 * (hx(i,j+1,k+1) - hx(i,j+1,k)) / dz;
    psi2 = b2 * psi2 + c2 * (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx;

    ey(i,j,k) += dt / epsilon * 
	((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz / kappa1 - 
	 (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx / kappa2 + psi1 - psi2);
}

void
CPMLEz::update(double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       double dt, double dx, double dy)
{
    psi1 = b1 * psi1 + c1 * (hy(i+1,j,k+1) - hy(i,j,k+1)) / dx;
    psi2 = b2 * psi2 + c2 * (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy;

    ez(i,j,k) += dt / epsilon * 
	((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx / kappa1 - 
	 (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy / kappa2 + psi1 - psi2);
}

CPMLMagnetic::CPMLMagnetic(const int * const idx, int size, double mu_r, 
			   double b1_in, double b2_in, 
			   double c1_in, double c2_in, 
			   double kappa1_in, double kappa2_in)
    : PointwiseMaterial(idx, size), 
      b1(b1_in), b2(b2_in), c1(c1_in), c2(c2_in), kappa1(kappa1_in), kappa2(kappa2_in), 
      psi1(0), psi2(0)
{
    mu = mu_r * mu0;
}

void 
CPMLMagnetic::set_mu(double mu_r)
{
    mu = mu_r * mu0;
}

void
CPMLHx::update(double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       double dt, double dy, double dz)
{
    psi1 = b1 * psi1 + c1 * (ez(i,j,k-1) - ez(i,j-1,k-1)) / dy;
    psi2 = b2 * psi2 + c2 * (ey(i,j-1,k) - ey(i,j-1,k-1)) / dz;
	
    hx(i,j,k) -= dt / mu * 
	((ez(i,j,k-1) - ez(i,j-1,k-1)) / dy / kappa1 - 
	 (ey(i,j-1,k) - ey(i,j-1,k-1)) / dz / kappa2 + psi1 - psi2);
}

void
CPMLHy::update(double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       double dt, double dz, double dx)
{
    psi1 = b1 * psi1 + c1 * (ex(i-1,j,k) - ex(i-1,j,k-1)) / dz;
    psi2 = b2 * psi2 + c2 * (ez(i,j,k-1) - ez(i-1,j,k-1)) / dx;
   
    hy(i,j,k) -= dt / mu * 
	((ex(i-1,j,k) - ex(i-1,j,k-1)) / dz / kappa1 - 
	 (ez(i,j,k-1) - ez(i-1,j,k-1)) / dx / kappa2 + psi1 - psi2);
}

void
CPMLHz::update(double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       double dt, double dx, double dy)
{
    psi1 = b1 * psi1 + c1 * (ey(i,j-1,k) - ey(i-1,j-1,k)) / dx;
    psi2 = b2 * psi2 + c2 * (ex(i-1,j,k) - ex(i-1,j-1,k)) / dy;

    hz(i,j,k) -= dt / mu * 
	((ey(i,j-1,k) - ey(i-1,j-1,k)) / dx / kappa1 - 
	 (ex(i-1,j,k) - ex(i-1,j-1,k)) / dy / kappa2 + psi1 - psi2);
}
