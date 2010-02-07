/* This implementation is based on the following article.
 * M. Okoniewski and E. Okoniewska, "Drude dispersion in ADE FDTD revisited,"
 * Electron. Lett., vol. 42, no. 9, pp. 503-504, 2006. */

#ifndef PW_DRUDE_HH_
#define PW_DRUDE_HH_

#include <vector>
#include <numeric>

#include "constants.hh"
#include "pw_material.hh"
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
	DrudeElectric(const int * const idx, int size, double epsilon_inf,
			const double * const omega_p, int omega_p_size,
			const double * const gamma_p, int gamma_p_size) :
		MaterialElectric<T>(idx, size), epsilon(epsilon_inf * epsilon0),
		omega_p(omega_p, omega_p + omega_p_size),
		gamma_p(gamma_p, gamma_p + gamma_p_size),
		q_new(omega_p_size, static_cast<T>(0.)),
		q_old(omega_p_size, static_cast<T>(0.))
	{
	}

	double get_epsilon()
		{
			return epsilon;
		}

	void set_epsilon(double epsilon_r)
		{
			epsilon = epsilon_r * epsilon0;
		}

protected:
	double epsilon;
	std::vector<double> omega_p;
	std::vector<double> gamma_p;
	std::vector<T > q_new;
	std::vector<T > q_old;
};

template <typename T> class DrudeEx: public DrudeElectric<T>
{
public:
	DrudeEx(const int * const idx, int size, double epsilon_inf,
			const double * const omega_p, int omega_p_size,
			const double * const gamma_p, int gamma_p_size) :
		DrudeElectric<T>(idx, size, epsilon_inf, omega_p, omega_p_size,
		gamma_p, gamma_p_size)
	{
	}

	void update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			double dy, double dz, double dt, double n)
	{
		std::vector<T > q_tmp(omega_p.size());

		for (unsigned int u = 0; u != q_tmp.size(); ++u)
		{
			q_tmp[u] = (4. * q_new[u] + (gamma_p[u] * dt - 2.) * q_old[u] - (2.
					* dt * dt * epsilon0 * omega_p[u] * omega_p[u]) * ex(i,j,k))
					/ (gamma_p[u] * dt + 2.);
		}

		std::copy(q_new.begin(), q_new.end(), q_old.begin());
		std::copy(q_tmp.begin(), q_tmp.end(), q_new.begin());

		T q_diff_sum = std::accumulate(q_new.begin(), q_new.end(), static_cast<T >(0.))
				- std::accumulate(q_old.begin(), q_old.end(), static_cast<T >(0.));

		ex(i,j,k) += (dt * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy
				- (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz) + q_diff_sum) / epsilon;
	}

protected:
	using DrudeElectric<T>::epsilon;
	using DrudeElectric<T>::omega_p;
	using DrudeElectric<T>::gamma_p;
	using DrudeElectric<T>::q_new;
	using DrudeElectric<T>::q_old;
};

template <typename T> class DrudeEy: public DrudeElectric<T>
{
public:
	DrudeEy(const int * const idx, int size, double epsilon_inf,
			const double * const omega_p, int omega_p_size,
			const double * const gamma_p, int gamma_p_size) :
		DrudeElectric<T>(idx, size, epsilon_inf, omega_p, omega_p_size, gamma_p,
				gamma_p_size)
	{
	}

	void update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			double dz, double dx, double dt, double n)
	{
		std::vector<T > q_tmp(omega_p.size());

		for (unsigned int u = 0; u != q_tmp.size(); ++u)
		{
			q_tmp[u] = (4. * q_new[u] + (gamma_p[u] * dt - 2.) * q_old[u] - (2.
					* dt * dt * epsilon0 * omega_p[u] * omega_p[u]) * ey(i,j,k))
					/ (gamma_p[u] * dt + 2.);
		}

		std::copy(q_new.begin(), q_new.end(), q_old.begin());
		std::copy(q_tmp.begin(), q_tmp.end(), q_new.begin());

		T q_diff_sum = std::accumulate(q_new.begin(), q_new.end(), static_cast<T >(0.))
				- std::accumulate(q_old.begin(), q_old.end(), static_cast<T >(0.));

		ey(i,j,k) += (dt * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz
				- (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx) + q_diff_sum) / epsilon;
	}

protected:
	using DrudeElectric<T>::epsilon;
	using DrudeElectric<T>::omega_p;
	using DrudeElectric<T>::gamma_p;
	using DrudeElectric<T>::q_new;
	using DrudeElectric<T>::q_old;
};

template <typename T> class DrudeEz: public DrudeElectric<T>
{
public:
	DrudeEz(const int * const idx, int size, double epsilon_inf,
			const double * const omega_p, int omega_p_size,
			const double * const gamma_p, int gamma_p_size) :
		DrudeElectric<T>(idx, size, epsilon_inf, omega_p, omega_p_size, gamma_p,
				gamma_p_size)
	{
	}

	void update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			double dx, double dy, double dt, double n)
	{
		std::vector<T > q_tmp(omega_p.size());

		for (unsigned int u = 0; u != q_tmp.size(); ++u)
		{
			q_tmp[u] = (4. * q_new[u] + (gamma_p[u] * dt - 2.) * q_old[u] - (2.
					* dt * dt * epsilon0 * omega_p[u] * omega_p[u]) * ez(i,j,k))
					/ (gamma_p[u] * dt + 2.);
		}

		std::copy(q_new.begin(), q_new.end(), q_old.begin());
		std::copy(q_tmp.begin(), q_tmp.end(), q_new.begin());

		T q_diff_sum = std::accumulate(q_new.begin(), q_new.end(), static_cast<T >(0.))
				- std::accumulate(q_old.begin(), q_old.end(), static_cast<T >(0.));

		ez(i,j,k) += (dt * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx
				- (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy) + q_diff_sum) / epsilon;
	}

protected:
	using DrudeElectric<T>::epsilon;
	using DrudeElectric<T>::omega_p;
	using DrudeElectric<T>::gamma_p;
	using DrudeElectric<T>::q_new;
	using DrudeElectric<T>::q_old;
};

template <typename T> class DrudeHx: public DielectricHx<T>
{
public:
	DrudeHx(const int * const idx, int size, double mu_r = 1) :
		DielectricHx<T>(idx, size, mu_r)
	{
	}
};

template <typename T> class DrudeHy: public DielectricHy<T>
{
public:
	DrudeHy(const int * const idx, int size, double mu_r = 1) :
		DielectricHy<T>(idx, size, mu_r)
	{
	}
};

template <typename T> class DrudeHz: public DielectricHz<T>
{
public:
	DrudeHz(const int * const idx, int size, double mu_r = 1) :
		DielectricHz<T>(idx, size, mu_r)
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
