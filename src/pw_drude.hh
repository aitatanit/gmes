/* This implementation is based on the following article.
 * M. Okoniewski and E. Okoniewska, "Drude dispersion in ADE FDTD revisited,"
 * Electron. Lett., vol. 42, no. 9, pp. 503-504, 2006. */

#ifndef PW_DRUDE_HH_
#define PW_DRUDE_HH_

#include <vector>

#include "pw_dielectric.hh"

#define ex(i,j,k) ex[((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[((i)*hz_y_size+(j))*hz_z_size+(k)]

namespace gmes
{
  template <typename T> class DrudeElectric: public MaterialElectric<T>
  {
  public:
    DrudeElectric(double epsilon,
		  const double * const a, int a_i_size, int a_j_size,
		  const double * const c, int c_size):
      eps(epsilon),
      c(c, c + c_size),
      q_new(a_i_size, static_cast<T>(0)),
      q_now(a_i_size, static_cast<T>(0))
    {
      for (int i = 0; i < a_i_size; i++)
	{
	  std::vector<double> tmp(a + i * a_j_size, a + (i + 1) * a_j_size);
	  this->a.push_back(tmp);
	}
    }

    double get_epsilon() const
    {
      return eps;
    }

    void set_epsilon(double epsilon)
    {
      eps = epsilon;
    }

    T dps_sum(const T& init)
    {
      T sum(init);
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i)
	{
	  sum += q_new[i] - q_now[i];
	}

      return sum;
    }

    void update_q(const T& e_now)
    {
      std::vector<T> q_old(q_now);
      std::copy(q_new.begin(), q_new.end(), q_now.begin());
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i)
	{
	  q_new[i] = a[i][0] * q_old[i] + a[i][1] * q_now[i] + a[i][2] * e_now;
	}
    }

  protected:
    double eps;
    std::vector<std::vector<double> > a;
    std::vector<double> c;
    std::vector<T> q_new;
    std::vector<T> q_now;
  };

  template <typename T> class DrudeEx: public DrudeElectric<T>
  {
  public:
    DrudeEx(double epsilon,
	    const double * const a, int a_i_size, int a_j_size,
	    const double * const c, int c_size):
      DrudeElectric<T>(epsilon, a, a_i_size, a_j_size, c, c_size)
    {
    }

    void update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
		const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
		const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
		double dy, double dz, double dt, double n, int i, int j, int k)
    {
      T e_now = ex(i,j,k);
      update_q(e_now);
      ex(i,j,k) = (c[0] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
			   (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz)
		   + c[1] * dps_sum(static_cast<T>(0)) + c[2] * e_now);
    }

  protected:
    using DrudeElectric<T>::eps;
    using DrudeElectric<T>::a;
    using DrudeElectric<T>::c;
    using DrudeElectric<T>::q_new;
    using DrudeElectric<T>::q_now;
  };

  template <typename T> class DrudeEy: public DrudeElectric<T>
  {
  public:
    DrudeEy(double epsilon,
	    const double * const a, int a_i_size, int a_j_size,
	    const double * const c, int c_size):
      DrudeElectric<T>(epsilon, a, a_i_size, a_j_size, c, c_size)
    {
    }

    void update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
		const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
		const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
		double dz, double dx, double dt, double n, int i, int j, int k)
    {
      T e_now = ey(i,j,k);
      update_q(e_now);
      ey(i,j,k) = (c[0] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
			   (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx)
		   + c[1] * dps_sum(static_cast<T>(0)) + c[2] * e_now);
    }

  protected:
    using DrudeElectric<T>::eps;
    using DrudeElectric<T>::a;
    using DrudeElectric<T>::c;
    using DrudeElectric<T>::q_new;
    using DrudeElectric<T>::q_now;
  };

  template <typename T> class DrudeEz: public DrudeElectric<T>
  {
  public:
    DrudeEz(double epsilon,
	    const double * const a, int a_i_size, int a_j_size,
	    const double * const c, int c_size):
      DrudeElectric<T>(epsilon, a, a_i_size, a_j_size, c, c_size)
    {
    }

    void update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
		const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
		const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
		double dx, double dy, double dt, double n, int i, int j, int k)
    {
      T e_now = ez(i,j,k);
      update_q(e_now);
      ez(i,j,k) = (c[0] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - 
			   (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy)
		   + c[1] * dps_sum(static_cast<T>(0)) + c[2] * e_now);
    }

  protected:
    using DrudeElectric<T>::eps;
    using DrudeElectric<T>::a;
    using DrudeElectric<T>::c;
    using DrudeElectric<T>::q_new;
    using DrudeElectric<T>::q_now;
  };

  template <typename T> class DrudeHx: public DielectricHx<T>
  {
  public:
    DrudeHx(double mu = 1):
      DielectricHx<T>(mu)
    {
    }
  };

  template <typename T> class DrudeHy: public DielectricHy<T>
  {
  public:
    DrudeHy(double mu = 1):
      DielectricHy<T>(mu)
    {
    }
  };

  template <typename T> class DrudeHz: public DielectricHz<T>
  {
  public:
    DrudeHz(double mu = 1):
      DielectricHz<T>(mu)
    {
    }
  };
}

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif /*PW_DRUDE_HH_*/
