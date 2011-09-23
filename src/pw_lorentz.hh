#ifndef PW_LORENTZ_HH_
#define PW_LORENTZ_HH_

#include <vector>

#include "pw_dielectric.hh"

#define ex(i,j,k) ex[ex_y_size==1?0:((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[ey_z_size==1?0:((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[ez_x_size==1?0:((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[hx_y_size==1?0:((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[hy_z_size==1?0:((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[hz_x_size==1?0:((i)*hz_y_size+(j))*hz_z_size+(k)]

namespace gmes
{
  template <typename T> class LorentzElectric: public MaterialElectric<T>
  {
  public:
    LorentzElectric(double epsilon,
		    const double * const a, int a_i_size, int a_j_size,
		    const double * const c, int c_size):
      eps(epsilon),
      c(c, c + c_size),
      l_new(a_i_size, static_cast<T>(0)),
      l_now(a_i_size, static_cast<T>(0))
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

    T lps_sum(const T& init)
    {
      T sum(init);
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i)
	{
	  sum += l_new[i] - l_now[i];
	}

      return sum;
    }

    void update_q(const T& e_now)
    {
      std::vector<T> l_old(l_now);
      std::copy(l_new.begin(), l_new.end(), l_now.begin());
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i)
	{
	  l_new[i] = a[i][0] * l_old[i] + a[i][1] * l_now[i] + a[i][2] * e_now;
	}
    }

  protected:
    double eps;
    std::vector<std::vector<double> > a;
    std::vector<double> c;
    std::vector<T> l_new;
    std::vector<T> l_now;
  };

  template <typename T> class LorentzEx: public LorentzElectric<T>
  {
  public:
    LorentzEx(double epsilon,
	      const double * const a, int a_i_size, int a_j_size,
	      const double * const c, int c_size):
      LorentzElectric<T>(epsilon, a, a_i_size, a_j_size, c, c_size)
    {
    }

    T update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	     const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	     const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	     double dy, double dz, double dt, double n, int i, int j, int k)
    {
      T e_now = ex(i,j,k);
      update_q(e_now);
      ex(i,j,k) = (c[0] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
			   (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz)
		   + c[1] * lps_sum(static_cast<T>(0)) + c[2] * e_now);

      return ex(i,j,k);
    }

  protected:
    using LorentzElectric<T>::eps;
    using LorentzElectric<T>::a;
    using LorentzElectric<T>::c;
    using LorentzElectric<T>::l_new;
    using LorentzElectric<T>::l_now;
  };

  template <typename T> class LorentzEy: public LorentzElectric<T>
  {
  public:
    LorentzEy(double epsilon,
	      const double * const a, int a_i_size, int a_j_size,
	      const double * const c, int c_size):
      LorentzElectric<T>(epsilon, a, a_i_size, a_j_size, c, c_size)
    {
    }

    T update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	     const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	     const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	     double dz, double dx, double dt, double n, int i, int j, int k)
    {
      T e_now = ey(i,j,k);
      update_q(e_now);
      ey(i,j,k) = (c[0] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
			   (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx)
		   + c[1] * lps_sum(static_cast<T>(0)) + c[2] * e_now);

      return ey(i,j,k);
    }

  protected:
    using LorentzElectric<T>::eps;
    using LorentzElectric<T>::a;
    using LorentzElectric<T>::c;
    using LorentzElectric<T>::l_new;
    using LorentzElectric<T>::l_now;
  };

  template <typename T> class LorentzEz: public LorentzElectric<T>
  {
  public:
    LorentzEz(double epsilon,
	      const double * const a, int a_i_size, int a_j_size,
	      const double * const c, int c_size):
      LorentzElectric<T>(epsilon, a, a_i_size, a_j_size, c, c_size)
    {
    }

    T update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	     const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	     const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	     double dx, double dy, double dt, double n, int i, int j, int k)
    {
      T e_now = ez(i,j,k);
      update_q(e_now);
      ez(i,j,k) = (c[0] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - 
			   (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy)
		   + c[1] * lps_sum(static_cast<T>(0)) + c[2] * e_now);

      return ez(i,j,k);
    }

  protected:
    using LorentzElectric<T>::eps;
    using LorentzElectric<T>::a;
    using LorentzElectric<T>::c;
    using LorentzElectric<T>::l_new;
    using LorentzElectric<T>::l_now;
  };

  template <typename T> class LorentzHx: public DielectricHx<T>
  {
  public:
    LorentzHx(double mu = 1):
      DielectricHx<T>(mu)
    {
    }
  };

  template <typename T> class LorentzHy: public DielectricHy<T>
  {
  public:
    LorentzHy(double mu = 1):
      DielectricHy<T>(mu)
    {
    }
  };

  template <typename T> class LorentzHz: public DielectricHz<T>
  {
  public:
    LorentzHz(double mu = 1):
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

#endif /*PW_LORENTZ_HH_*/
