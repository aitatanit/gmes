#ifndef PW_DIELECTRIC_HH_
#define PW_DIELECTRIC_HH_

#include "pw_material.hh"

#define ex(i,j,k) ex[((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[((i)*hz_y_size+(j))*hz_z_size+(k)]

namespace gmes
{
  template <typename T> class DielectricElectric: public MaterialElectric<T>
  {
  public:
    DielectricElectric(double epsilon = 1):
      eps(epsilon)
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
  };

  template <typename T> class DielectricEx: public DielectricElectric<T>
  {
  public:
    DielectricEx(double epsilon = 1):
      DielectricElectric<T>(epsilon)
    {
    }

    void update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
		const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
		const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
		double dy, double dz, double dt, double n, int i, int j, int k)
    {
      ex(i,j,k) += dt / eps * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
			       (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz);
    }

  protected:
    using DielectricElectric<T>::eps;
  };

  template <typename T> class DielectricEy: public DielectricElectric<T>
  {
  public:
    DielectricEy(double epsilon = 1):
      DielectricElectric<T>(epsilon)
    {
    }

    void update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
		const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
		const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
		double dz, double dx, double dt, double n, int i, int j, int k)
    {
      ey(i,j,k) += dt / eps * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
			       (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx);
    }

  protected:
    using DielectricElectric<T>::eps;
  };

  template <typename T> class DielectricEz: public DielectricElectric<T>
  {
  public:
    DielectricEz(double epsilon = 1):
      DielectricElectric<T>(epsilon)
    {
    }

    void update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
		const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
		const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
		double dx, double dy, double dt, double n, int i, int j, int k)
    {
      ez(i,j,k) += dt / eps * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx -
			       (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy);
    }

  protected:
    using DielectricElectric<T>::eps;
  };

  template <typename T> class DielectricMagnetic: public MaterialMagnetic<T>
  {
  public:
    DielectricMagnetic(double mu = 1):
      mu(mu)
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
  };

  template <typename T> class DielectricHx: public DielectricMagnetic<T>
  {
  public:
    DielectricHx(double mu = 1):
      DielectricMagnetic<T>(mu)
    {
    }

    void update(T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
		const T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
		const T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
		double dy, double dz, double dt, double n, int i, int j, int k)
    {
      hx(i,j,k) += dt / mu * ((ey(i,j-1,k) - ey(i,j-1,k-1)) / dz -
			      (ez(i,j,k-1) - ez(i,j-1,k-1)) / dy);
    }

  protected:
    using DielectricMagnetic<T>::mu;
  };

  template <typename T> class DielectricHy: public DielectricMagnetic<T>
  {
  public:
    DielectricHy(double mu = 1):
      DielectricMagnetic<T>(mu)
    {
    }

    void update(T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
		const T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
		const T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
		double dz, double dx, double dt, double n, int i, int j, int k)
    {
      hy(i,j,k) += dt / mu * ((ez(i,j,k-1) - ez(i-1,j,k-1)) / dx -
			      (ex(i-1,j,k) - ex(i-1,j,k-1)) / dz);
    }

  protected:
    using DielectricMagnetic<T>::mu;
  };

  template <typename T> class DielectricHz: public DielectricMagnetic<T>
  {
  public:
    DielectricHz(double mu = 1):
      DielectricMagnetic<T>(mu)
    {
    }

    void update(T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
		const T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
		const T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
		double dx, double dy, double dt, double n, int i, int j, int k)
    {
      hz(i,j,k) += dt / mu * ((ex(i-1,j,k) - ex(i-1,j-1,k)) / dy -
			      (ey(i,j-1,k) - ey(i-1,j-1,k)) / dx);
    }

  protected:
    using DielectricMagnetic<T>::mu;
  };
}

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif /*PW_DIELECTRIC_HH_*/
