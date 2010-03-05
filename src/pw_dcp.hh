/* The implementation of Drude-critical points model based on the following article.
 * P. G. Etchegoin, E. C. Le Ru, and M. Meyer, "An analytic model for the
 * optical properties of gold," J. Chem. Phys. 125, 164705, 2001.
 */

#ifndef PW_DCP_HH_
#define PW_DCP_HH_

#include <vector>
#include <numeric>
#include <iostream>

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
template <typename T> class DCPElectric: public MaterialElectric<T>
{
public:
	DCPElectric(const int * const idx, int size, double epsilon,
			const double * const a, int a_i_size, int a_j_size,
			const double * const b, int b_i_size, int b_j_size,
			const double * const c, int c_size) :
		MaterialElectric<T>(idx, size), epsilon(epsilon),
		c(c, c + c_size),
		e_old(static_cast<T >(0.)),
		q_now(a_i_size, static_cast<T >(0.)),
		q_old(a_i_size, static_cast<T >(0.)),
		p_now(b_i_size, static_cast<T >(0.)),
		p_old(b_i_size, static_cast<T >(0.))
	{
		for (int i = 0; i < a_i_size; i++)
		{
			std::vector<double> tmp(a + i * a_j_size, a + (i + 1) * a_j_size);
			this->a.push_back(tmp);
		}

		for (int i = 0; i < b_i_size; i++)
		{
			std::vector<double> tmp(b + i * b_j_size, b + (i + 1) * b_j_size);
			this->b.push_back(tmp);
		}
	}

	double get_epsilon()
		{
			return epsilon;
		}

	void set_epsilon(double epsilon)
		{
			this->epsilon = epsilon;
		}

protected:
	double epsilon;
	std::vector<std::vector<double> > a;
	std::vector<std::vector<double> > b;
	std::vector<double> c;
	T e_old;
	std::vector<T > q_now;
	std::vector<T > q_old;
	std::vector<T > p_now;
	std::vector<T > p_old;
};

template <typename T> class DCPEx: public DCPElectric<T>
{
public:
	DCPEx(const int * const idx, int size, double epsilon,
			const double * const a, int a_i_size, int a_j_size,
			const double * const b, int b_i_size, int b_j_size,
			const double * const c, int c_size) :
				DCPElectric<T>(idx, size, epsilon, a, a_i_size, a_j_size, b, b_i_size, b_j_size, c, c_size)
	{
	}

	void update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			double dy, double dz, double dt, double n)
	{
		T dps_sum = static_cast<T >(0.);
		for (unsigned int i = 0; i < a.size(); ++i)
		{
			dps_sum += a[i][0] * q_old[i] + (a[i][1] - 1) * q_now[i];
		}

		T cps_sum = static_cast<T >(0.);
		for (unsigned int i = 0; i < b.size(); ++i)
		{
			cps_sum += b[i][0] * p_old[i] + (b[i][1] - 1) * p_now[i];
		}

		T e_now = ex(i,j,k);
		T e_new = c[0] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz)
				+ c[1] * (dps_sum - cps_sum) - c[2] * e_old + c[3] * e_now;

		std::vector<T > q_new(a.size());
		for (unsigned int i = 0; i < q_new.size(); ++i)
		{
			q_new[i] = a[i][0] * q_old[i] + a[i][1] * q_now[i] - a[i][2] * e_now;
		}
		std::copy(q_now.begin(), q_now.end(), q_old.begin());
		std::copy(q_new.begin(), q_new.end(), q_now.begin());

		std::vector<T > p_new(b.size());
		for (unsigned int i = 0; i < p_new.size(); ++i)
		{
			p_new[i] = b[i][0] * p_old[i] + b[i][1] * p_now[i] + b[i][2] * (e_old - e_new) + b[i][3] * e_now;
		}
		std::copy(p_now.begin(), p_now.end(), p_old.begin());
		std::copy(p_new.begin(), p_new.end(), p_now.begin());

		e_old = e_now;
		ex(i,j,k) = e_new;
	}

protected:
	using DCPElectric<T >::epsilon;
	using DCPElectric<T >::a;
	using DCPElectric<T >::b;
	using DCPElectric<T >::c;
	using DCPElectric<T >::e_old;
	using DCPElectric<T >::p_now;
	using DCPElectric<T >::p_old;
	using DCPElectric<T >::q_now;
	using DCPElectric<T >::q_old;
};

template <typename T> class DCPEy: public DCPElectric<T>
{
public:
	DCPEy(const int * const idx, int size, double epsilon,
			const double * const a, int a_i_size, int a_j_size,
			const double * const b, int b_i_size, int b_j_size,
			const double * const c, int c_size) :
				DCPElectric<T>(idx, size, epsilon, a, a_i_size, a_j_size, b, b_i_size, b_j_size, c, c_size)
	{
	}

	void update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			double dz, double dx, double dt, double n)
	{
		T dps_sum = static_cast<T >(0.);
		for (unsigned int i = 0; i < a.size(); ++i)
		{
			dps_sum += a[i][0] * q_old[i] + (a[i][1] - 1) * q_now[i];
		}

		T cps_sum = static_cast<T >(0.);
		for (unsigned int i = 0; i < b.size(); ++i)
		{
			cps_sum += b[i][0] * p_old[i] + (b[i][1] - 1) * p_now[i];
		}

		T e_now = ey(i,j,k);
		T e_new = c[0] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx)
				+ c[1] * (dps_sum - cps_sum) - c[2] * e_old + c[3] * e_now;

		std::vector<T > q_new(a.size());
		for (unsigned int i = 0; i < q_new.size(); ++i)
		{
			q_new[i] = a[i][0] * q_old[i] + a[i][1] * q_now[i] - a[i][2] * e_now;
		}
		std::copy(q_now.begin(), q_now.end(), q_old.begin());
		std::copy(q_new.begin(), q_new.end(), q_now.begin());

		std::vector<T > p_new(b.size());
		for (unsigned int i = 0; i < p_new.size(); ++i)
		{
			p_new[i] = b[i][0] * p_old[i] + b[i][1] * p_now[i] + b[i][2] * (e_old - e_new) + b[i][3] * e_now;
		}
		std::copy(p_now.begin(), p_now.end(), p_old.begin());
		std::copy(p_new.begin(), p_new.end(), p_now.begin());

		e_old = e_now;
		ey(i,j,k) = e_new;
	}

protected:
	using DCPElectric<T >::epsilon;
	using DCPElectric<T >::a;
	using DCPElectric<T >::b;
	using DCPElectric<T >::c;
	using DCPElectric<T >::e_old;
	using DCPElectric<T >::p_now;
	using DCPElectric<T >::p_old;
	using DCPElectric<T >::q_now;
	using DCPElectric<T >::q_old;
};

template <typename T> class DCPEz: public DCPElectric<T>
{
public:
	DCPEz(const int * const idx, int size, double epsilon,
			const double * const a, int a_i_size, int a_j_size,
			const double * const b, int b_i_size, int b_j_size,
			const double * const c, int c_size) :
				DCPElectric<T>(idx, size, epsilon, a, a_i_size, a_j_size, b, b_i_size, b_j_size, c, c_size)
	{
	}

	void update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			double dx, double dy, double dt, double n)
	{
		T dps_sum = static_cast<T >(0.);
		for (unsigned int i = 0; i < a.size(); ++i)
		{
			dps_sum += a[i][0] * q_old[i] + (a[i][1] - 1) * q_now[i];
		}

		T cps_sum = static_cast<T >(0.);
		for (unsigned int i = 0; i < b.size(); ++i)
		{
			cps_sum += b[i][0] * p_old[i] + (b[i][1] - 1) * p_now[i];
		}

		T e_now = ez(i,j,k);
		T e_new = c[0] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy)
				+ c[1] * (dps_sum - cps_sum) - c[2] * e_old + c[3] * e_now;

		std::vector<T > q_new(a.size());
		for (unsigned int i = 0; i < a.size(); ++i)
		{
			q_new[i] = a[i][0] * q_old[i] + a[i][1] * q_now[i] - a[i][2] * e_now;
		}
		std::copy(q_now.begin(), q_now.end(), q_old.begin());
		std::copy(q_new.begin(), q_new.end(), q_now.begin());

		std::vector<T > p_new(b.size());
		for (unsigned int i = 0; i < b.size(); ++i)
		{
			p_new[i] = b[i][0] * p_old[i] + b[i][1] * p_now[i] + b[i][2] * (e_old - e_new) + b[i][3] * e_now;
		}
		std::copy(p_now.begin(), p_now.end(), p_old.begin());
		std::copy(p_new.begin(), p_new.end(), p_now.begin());

		e_old = e_now;
		ez(i,j,k) = e_new;
	}

protected:
	using DCPElectric<T >::epsilon;
	using DCPElectric<T >::a;
	using DCPElectric<T >::b;
	using DCPElectric<T >::c;
	using DCPElectric<T >::e_old;
	using DCPElectric<T >::p_now;
	using DCPElectric<T >::p_old;
	using DCPElectric<T >::q_now;
	using DCPElectric<T >::q_old;
};

template <typename T> class DCPHx: public DielectricHx<T>
{
public:
	DCPHx(const int * const idx, int size, double mu = 1) :
		DielectricHx<T>(idx, size, mu)
	{
	}
};

template <typename T> class DCPHy: public DielectricHy<T>
{
public:
	DCPHy(const int * const idx, int size, double mu = 1) :
		DielectricHy<T>(idx, size, mu)
	{
	}
};

template <typename T> class DCPHz: public DielectricHz<T>
{
public:
	DCPHz(const int * const idx, int size, double mu = 1) :
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

#endif /*PW_DCP_HH_*/
