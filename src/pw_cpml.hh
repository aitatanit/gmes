/* This implementation is based on the following article.
 * S. Gedney, “Perfectly Matched Layer Absorbing Boundary Conditions,”
 * Computational Electrodynamics: The Finite-Difference Time-Domain Method,
 * Third Edition, A. Taflove and S.C. Hagness, eds., Artech House Publishers,
 *  2005, pp. 273-328.
 */ 
#ifndef PW_CPML_HH_
#define PW_CPML_HH_

#include "pw_material.hh"

#define ex(i,j,k) ex[((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[((i)*hz_y_size+(j))*hz_z_size+(k)]

namespace gmes
{
  template <typename T> class CpmlElectric: public MaterialElectric<T>
  {
  public:
    CpmlElectric(double epsilon,
		 double b1_in, double b2_in, double c1_in, double c2_in,
		 double kappa1_in, double kappa2_in):
      eps(epsilon),
      b1(b1_in), b2(b2_in), c1(c1_in), c2(c2_in),
      kappa1(kappa1_in), kappa2(kappa2_in), psi1(0), psi2(0)
    {
    }

    double get_epsilon() const
    {
      return eps;
    }

    void set_epsilon(double epsilon)
    {
      eps = epsilon;
    }

  protected:
    double eps;
    double b1, b2;
    double c1, c2;
    double kappa1, kappa2;
    T psi1, psi2;
  };

  template <typename T> class CpmlEx: public CpmlElectric<T>
  {
  public:
    CpmlEx(double epsilon, double by,
	   double bz, double cy, double cz, double kappay, double kappaz):
      CpmlElectric<T>(epsilon, by, bz, cy, cz, kappay, kappaz)
    {
    }

    void update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
		const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
		const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
		double dy, double dz, double dt, double n, int i, int j, int k)
    {
      psi1 = b1 * psi1 + c1 * (hz(i+1,j+1,k) - hz(i+1,j,k)) / dy;
      psi2 = b2 * psi2 + c2 * (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz;

      ex(i,j,k) += dt / eps * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy / kappa1 -
			       (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz / kappa2 
			       + psi1 - psi2);
    }

  protected:
    using CpmlElectric<T>::eps;
    using CpmlElectric<T>::b1;
    using CpmlElectric<T>::b2;
    using CpmlElectric<T>::c1;
    using CpmlElectric<T>::c2;
    using CpmlElectric<T>::kappa1;
    using CpmlElectric<T>::kappa2;
    using CpmlElectric<T>::psi1;
    using CpmlElectric<T>::psi2;
  };

  template <typename T> class CpmlEy: public CpmlElectric<T>
  {
  public:
    CpmlEy(double epsilon, double bz,
	   double bx, double cz, double cx, double kappaz, double kappax):
      CpmlElectric<T>(epsilon, bz, bx, cz, cx, kappaz, kappax)
    {
    }

    void update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
		const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
		const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
		double dz, double dx, double dt, double n, int i, int j, int k)
    {
      psi1 = b1 * psi1 + c1 * (hx(i,j+1,k+1) - hx(i,j+1,k)) / dz;
      psi2 = b2 * psi2 + c2 * (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx;

      ey(i,j,k) += dt / eps * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz / kappa1 - 
			       (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx / kappa2 
			       + psi1 - psi2);
    }

  protected:
    using CpmlElectric<T>::eps;
    using CpmlElectric<T>::b1;
    using CpmlElectric<T>::b2;
    using CpmlElectric<T>::c1;
    using CpmlElectric<T>::c2;
    using CpmlElectric<T>::kappa1;
    using CpmlElectric<T>::kappa2;
    using CpmlElectric<T>::psi1;
    using CpmlElectric<T>::psi2;
  };

  template <typename T> class CpmlEz: public CpmlElectric<T>
  {
  public:
    CpmlEz(double epsilon, double bx,
	   double by, double cx, double cy, double kappax, double kappay):
      CpmlElectric<T>(epsilon, bx, by, cx, cy, kappax, kappay)
    {
    }

    void update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
		const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
		const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
		double dx, double dy, double dt, double n, int i, int j, int k)
    {
      psi1 = b1 * psi1 + c1 * (hy(i+1,j,k+1) - hy(i,j,k+1)) / dx;
      psi2 = b2 * psi2 + c2 * (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy;

      ez(i,j,k) += dt / eps * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx / kappa1 - 
			       (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy / kappa2 
			       + psi1 - psi2);
    }

  protected:
    using CpmlElectric<T>::eps;
    using CpmlElectric<T>::b1;
    using CpmlElectric<T>::b2;
    using CpmlElectric<T>::c1;
    using CpmlElectric<T>::c2;
    using CpmlElectric<T>::kappa1;
    using CpmlElectric<T>::kappa2;
    using CpmlElectric<T>::psi1;
    using CpmlElectric<T>::psi2;
  };

  template <typename T> class CpmlMagnetic: public MaterialMagnetic<T>
  {
  public:
    CpmlMagnetic(double mu, double b1_in,
		 double b2_in, double c1_in, double c2_in, double kappa1_in,
		 double kappa2_in):
      mu(mu),
      b1(b1_in), b2(b2_in), c1(c1_in), c2(c2_in),
      kappa1(kappa1_in), kappa2(kappa2_in), psi1(0), psi2(0)
    {
    }

    double get_mu() const
    {
      return mu;
    }

    void set_mu(double mu)
    {
      this->mu = mu;
    }

  protected:
    double mu;
    double b1, b2;
    double c1, c2;
    double kappa1, kappa2;
    T psi1, psi2;
  };

  template <typename T> class CpmlHx: public CpmlMagnetic<T>
  {
  public:
    CpmlHx(double mu, double by, double bz,
	   double cy, double cz, double kappay, double kappaz):
      CpmlMagnetic<T>(mu, by, bz, cy, cz, kappay, kappaz)
    {
    }

    void update(T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
		const T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
		const T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
		double dy, double dz, double dt, double n, int i, int j, int k)
    {
      psi1 = b1 * psi1 + c1 * (ez(i,j,k-1) - ez(i,j-1,k-1)) / dy;
      psi2 = b2 * psi2 + c2 * (ey(i,j-1,k) - ey(i,j-1,k-1)) / dz;

      hx(i,j,k) -= dt / mu
	* ((ez(i,j,k-1) - ez(i,j-1,k-1)) / dy / kappa1 
	   - (ey(i,j-1,k) - ey(i,j-1,k-1)) / dz / kappa2 + psi1 - psi2);
    }

  protected:
    using CpmlMagnetic<T>::mu;
    using CpmlMagnetic<T>::b1;
    using CpmlMagnetic<T>::b2;
    using CpmlMagnetic<T>::c1;
    using CpmlMagnetic<T>::c2;
    using CpmlMagnetic<T>::kappa1;
    using CpmlMagnetic<T>::kappa2;
    using CpmlMagnetic<T>::psi1;
    using CpmlMagnetic<T>::psi2;
  };

  template <typename T> class CpmlHy: public CpmlMagnetic<T>
  {
  public:
    CpmlHy(double mu, double bz, double bx,
	   double cz, double cx, double kappaz, double kappax):
      CpmlMagnetic<T>(mu, bz, bx, cz, cx, kappaz, kappax)
    {
    }

    void update(T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
		const T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
		const T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
		double dz, double dx, double dt, double n, int i, int j, int k)
    {
      psi1 = b1 * psi1 + c1 * (ex(i-1,j,k) - ex(i-1,j,k-1)) / dz;
      psi2 = b2 * psi2 + c2 * (ez(i,j,k-1) - ez(i-1,j,k-1)) / dx;

      hy(i,j,k) -= dt / mu
	* ((ex(i-1,j,k) - ex(i-1,j,k-1)) / dz / kappa1 
	   - (ez(i,j,k-1) - ez(i-1,j,k-1)) / dx / kappa2 + psi1 - psi2);
    }

  protected:
    using CpmlMagnetic<T>::mu;
    using CpmlMagnetic<T>::b1;
    using CpmlMagnetic<T>::b2;
    using CpmlMagnetic<T>::c1;
    using CpmlMagnetic<T>::c2;
    using CpmlMagnetic<T>::kappa1;
    using CpmlMagnetic<T>::kappa2;
    using CpmlMagnetic<T>::psi1;
    using CpmlMagnetic<T>::psi2;
  };

  template <typename T> class CpmlHz: public CpmlMagnetic<T>
  {
  public:
    CpmlHz(double mu, double bx, double by,
	   double cx, double cy, double kappax, double kappay):
      CpmlMagnetic<T>(mu, bx, by, cx, cy, kappax, kappay)
    {
    }

    void update(T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
		const T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
		const T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
		double dx, double dy, double dt, double n, int i, int j, int k)
    {
      psi1 = b1 * psi1 + c1 * (ey(i,j-1,k) - ey(i-1,j-1,k)) / dx;
      psi2 = b2 * psi2 + c2 * (ex(i-1,j,k) - ex(i-1,j-1,k)) / dy;

      hz(i,j,k) -= dt / mu
	* ((ey(i,j-1,k) - ey(i-1,j-1,k)) / dx / kappa1 
	   - (ex(i-1,j,k) - ex(i-1,j-1,k)) / dy / kappa2 + psi1 - psi2);
    }

  protected:
    using CpmlMagnetic<T>::mu;
    using CpmlMagnetic<T>::b1;
    using CpmlMagnetic<T>::b2;
    using CpmlMagnetic<T>::c1;
    using CpmlMagnetic<T>::c2;
    using CpmlMagnetic<T>::kappa1;
    using CpmlMagnetic<T>::kappa2;
    using CpmlMagnetic<T>::psi1;
    using CpmlMagnetic<T>::psi2;
  };
}

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif /*PW_CPML_HH_*/
