/* This implementation is based on the following article.
 * M. Okoniewski and E. Okoniewska, "Drude dispersion in ADE FDTD revisited,"
 * Electron. Lett., vol. 42, no. 9, pp. 503-504, 2006. */

#ifndef PW_DRUDE_HH_
#define PW_DRUDE_HH_

#include <vector>

#include "pw_dielectric.hh"

#define ex(i,j,k) ex[((this->i)*ex_y_size+(this->j))*ex_z_size+(this->k)]
#define ey(i,j,k) ey[((this->i)*ey_y_size+(this->j))*ey_z_size+(this->k)]
#define ez(i,j,k) ez[((this->i)*ez_y_size+(this->j))*ez_z_size+(this->k)]
#define hx(i,j,k) hx[((this->i)*hx_y_size+(this->j))*hx_z_size+(this->k)]
#define hy(i,j,k) hy[((this->i)*hy_y_size+(this->j))*hy_z_size+(this->k)]
#define hz(i,j,k) hz[((this->i)*hz_y_size+(this->j))*hz_z_size+(this->k)]

namespace gmes
{
template <typename T> class DrudeElectric: public MaterialElectric<T>
{
public:
	DrudeElectric(const int * const idx, int size, double epsilon,
			const double * const a, int a_i_size, int a_j_size,
			const double * const c, int c_size) :
		MaterialElectric<T>(idx, size), epsilon(epsilon),
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

	double get_epsilon()
		{
			return epsilon;
		}

	void set_epsilon(double epsilon)
		{
			epsilon = epsilon;
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
	double epsilon;
	std::vector<std::vector<double> > a;
	std::vector<double> c;
	std::vector<T> q_new;
	std::vector<T> q_now;
};

template <typename T> class DrudeEx: public DrudeElectric<T>
{
public:
	DrudeEx(const int * const idx, int size, double epsilon,
			const double * const a, int a_i_size, int a_j_size,
			const double * const c, int c_size) :
		DrudeElectric<T>(idx, size, epsilon, a, a_i_size, a_j_size, c, c_size)
	{
	}

	void update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			double dy, double dz, double dt, double n)
	{
		T e_now = ex(i,j,k);
		update_q(e_now);
		ex(i,j,k) = c[0] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz)
				+ c[1] * dps_sum(static_cast<T>(0)) + c[2] * e_now;
	}

protected:
	using DrudeElectric<T>::epsilon;
	using DrudeElectric<T>::a;
	using DrudeElectric<T>::c;
	using DrudeElectric<T>::q_new;
	using DrudeElectric<T>::q_now;
};

template <typename T> class DrudeEy: public DrudeElectric<T>
{
public:
	DrudeEy(const int * const idx, int size, double epsilon,
			const double * const a, int a_i_size, int a_j_size,
			const double * const c, int c_size) :
		DrudeElectric<T>(idx, size, epsilon, a, a_i_size, a_j_size, c, c_size)
	{
	}

	void update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			double dz, double dx, double dt, double n)
	{
		T e_now = ey(i,j,k);
		update_q(e_now);
		ey(i,j,k) = c[0] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx)
				+ c[1] * dps_sum(static_cast<T>(0)) + c[2] * e_now;
	}

protected:
	using DrudeElectric<T>::epsilon;
	using DrudeElectric<T>::a;
	using DrudeElectric<T>::c;
	using DrudeElectric<T>::q_new;
	using DrudeElectric<T>::q_now;
};

template <typename T> class DrudeEz: public DrudeElectric<T>
{
public:
	DrudeEz(const int * const idx, int size, double epsilon,
			const double * const a, int a_i_size, int a_j_size,
			const double * const c, int c_size) :
		DrudeElectric<T>(idx, size, epsilon, a, a_i_size, a_j_size, c, c_size)
	{
	}

	void update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			double dx, double dy, double dt, double n)
	{
		T e_now = ez(i,j,k);
		update_q(e_now);
		ez(i,j,k) = c[0] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy)
				+ c[1] * dps_sum(static_cast<T>(0)) + c[2] * e_now;
	}

protected:
	using DrudeElectric<T>::epsilon;
	using DrudeElectric<T>::a;
	using DrudeElectric<T>::c;
	using DrudeElectric<T>::q_new;
	using DrudeElectric<T>::q_now;
};

template <typename T> class DrudeHx: public DielectricHx<T>
{
public:
	DrudeHx(const int * const idx, int size, double mu = 1) :
		DielectricHx<T>(idx, size, mu)
	{
	}
};

template <typename T> class DrudeHy: public DielectricHy<T>
{
public:
	DrudeHy(const int * const idx, int size, double mu = 1) :
		DielectricHy<T>(idx, size, mu)
	{
	}
};

template <typename T> class DrudeHz: public DielectricHz<T>
{
public:
	DrudeHz(const int * const idx, int size, double mu = 1) :
		DielectricHz<T>(idx, size, mu)
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
