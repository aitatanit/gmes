/* The implementation of Drude-critical points model based on the following article.
 * P. G. Etchegoin, E. C. Le Ru, and M. Meyer, "An analytic model for the
 * optical properties of gold," J. Chem. Phys. 125, 164705, 2001.
 */

#ifndef PW_DCP_HH_
#define PW_DCP_HH_

#include <complex>
#include <numeric>
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
/**************************************************/
/* Auxiliary Differential Equation Implementation */
/**************************************************/

template <typename T> class DCPElectric: public MaterialElectric<T>
{
public:
	DCPElectric(const int * const idx, int size, double epsilon,
			const double * const a, int a_i_size, int a_j_size,
			const double * const b, int b_i_size, int b_j_size,
			const double * const c, int c_size) :
		MaterialElectric<T>(idx, size), epsilon(epsilon),
		c(c, c + c_size),
		e_old(static_cast<T>(0)),
		q_now(a_i_size, static_cast<T>(0)),
		q_old(a_i_size, static_cast<T>(0)),
		p_now(b_i_size, static_cast<T>(0)),
		p_old(b_i_size, static_cast<T>(0))
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

	double get_epsilon() const
	{
		return epsilon;
	}

	void set_epsilon(double epsilon)
	{
		this->epsilon = epsilon;
	}

	T dps_sum(const T& init)
	{
		T sum(init);
		for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i)
		{
			sum += a[i][0] * q_old[i] + (a[i][1] - 1) * q_now[i];
		}

		return sum;
	}

	T cps_sum(const T& init)
	{
		T sum(init);
		for (typename std::vector<T>::size_type i = 0; i < b.size(); ++i)
		{
			sum += b[i][0] * p_old[i] + (b[i][1] - 1) * p_now[i];
		}

		return sum;
	}

	void update_q(const T& e_now)
	{
		std::vector<T> q_new(a.size());
		for (typename std::vector<T>::size_type i = 0; i < q_new.size(); ++i)
		{
			q_new[i] = a[i][0] * q_old[i] + a[i][1] * q_now[i] + a[i][2] * e_now;
		}
		std::copy(q_now.begin(), q_now.end(), q_old.begin());
		std::copy(q_new.begin(), q_new.end(), q_now.begin());
	}

	void update_p(const T& e_old, const T& e_now, const T& e_new)
	{
		std::vector<T> p_new(b.size());
		for (typename std::vector<T>::size_type i = 0; i < p_new.size(); ++i)
		{
			p_new[i] = b[i][0] * p_old[i] + b[i][1] * p_now[i] + b[i][2] * (e_old - e_new) + b[i][3] * e_now;
		}
		std::copy(p_now.begin(), p_now.end(), p_old.begin());
		std::copy(p_new.begin(), p_new.end(), p_now.begin());
	}

protected:
	double epsilon;
	std::vector<std::vector<double> > a;
	std::vector<std::vector<double> > b;
	std::vector<double> c;
	T e_old;
	std::vector<T> q_now;
	std::vector<T> q_old;
	std::vector<T> p_now;
	std::vector<T> p_old;
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
		T e_now = ex(i,j,k);
		T e_new = c[0] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz)
				+ c[1] * (dps_sum(static_cast<T>(0)) + cps_sum(static_cast<T>(0)))
				+ c[2] * e_old + c[3] * e_now;

		update_q(e_now);
		update_p(e_old, e_now, e_new);

		e_old = e_now;
		ex(i,j,k) = e_new;
	}

protected:
	using DCPElectric<T>::epsilon;
	using DCPElectric<T>::a;
	using DCPElectric<T>::b;
	using DCPElectric<T>::c;
	using DCPElectric<T>::e_old;
	using DCPElectric<T>::p_now;
	using DCPElectric<T>::p_old;
	using DCPElectric<T>::q_now;
	using DCPElectric<T>::q_old;
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
		T e_now = ey(i,j,k);
		T e_new = c[0] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx)
				+ c[1] * (dps_sum(static_cast<T>(0)) + cps_sum(static_cast<T>(0)))
				+ c[2] * e_old + c[3] * e_now;

		update_q(e_now);
		update_p(e_old, e_now, e_new);

		e_old = e_now;
		ey(i,j,k) = e_new;
	}

protected:
	using DCPElectric<T>::epsilon;
	using DCPElectric<T>::a;
	using DCPElectric<T>::b;
	using DCPElectric<T>::c;
	using DCPElectric<T>::e_old;
	using DCPElectric<T>::p_now;
	using DCPElectric<T>::p_old;
	using DCPElectric<T>::q_now;
	using DCPElectric<T>::q_old;
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
		T e_now = ez(i,j,k);
		T e_new = c[0] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy)
				+ c[1] * (dps_sum(static_cast<T>(0)) + cps_sum(static_cast<T>(0)))
				+ c[2] * e_old + c[3] * e_now;

		update_q(e_now);
		update_p(e_old, e_now, e_new);

		e_old = e_now;
		ez(i,j,k) = e_new;
	}

protected:
	using DCPElectric<T>::epsilon;
	using DCPElectric<T>::a;
	using DCPElectric<T>::b;
	using DCPElectric<T>::c;
	using DCPElectric<T>::e_old;
	using DCPElectric<T>::p_now;
	using DCPElectric<T>::p_old;
	using DCPElectric<T>::q_now;
	using DCPElectric<T>::q_old;
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

/*********************************************************/
/* Piecewise-Linear Recursive Convolution Implementation */
/*********************************************************/

template <typename T> class DCPPLRCElectric: public MaterialElectric<T>
{
public:
	DCPPLRCElectric(const int * const idx, int size, double epsilon,
			const double * const a, int a_i_size, int a_j_size,
			const std::complex<double> * const b_c, int b_c_i_size, int b_c_j_size,
			const double * const c, int c_size) :
		MaterialElectric<T>(idx, size), epsilon(epsilon),
		c(c, c + c_size),
		psi_dp_re(a_i_size, .0),
		psi_dp_im(a_i_size, .0),
		psi_cp_re(b_c_i_size, std::complex<double>(0)),
		psi_cp_im(b_c_i_size, std::complex<double>(0))
	{
		for (int i = 0; i < a_i_size; i++)
		{
			std::vector<double> tmp(a + i * a_j_size, a + (i + 1) * a_j_size);
			this->a.push_back(tmp);
		}

		for (int i = 0; i < b_c_i_size; i++)
		{
			std::vector<std::complex<double> > tmp(b_c + i * b_c_j_size, b_c + (i + 1) * b_c_j_size);
			this->b.push_back(tmp);
		}
	}

	double get_epsilon() const
	{
		return epsilon;
	}

	void set_epsilon(double epsilon)
	{
		this->epsilon = epsilon;
	}

	void update_psi_dp(const std::complex<double>& e_now, const std::complex<double>& e_new)
	{
		for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i)
		{
			psi_dp_re[i] = a[i][0] * e_new.real() + a[i][1] * e_now.real() + a[i][2] * psi_dp_re[i];
			psi_dp_im[i] = a[i][0] * e_new.imag() + a[i][1] * e_now.imag() + a[i][2] * psi_dp_im[i];
		}
	}

	void update_psi_cp(const std::complex<double>& e_now, const std::complex<double>& e_new)
	{
		for (typename std::vector<std::complex<double> >::size_type i = 0; i < b.size(); ++i)
		{
			psi_cp_re[i] = b[i][0] * e_new.real() + b[i][1] * e_now.real()
					+ b[i][2] * psi_cp_re[i];
			psi_cp_im[i] = b[i][0] * e_new.imag() + b[i][1] * e_now.imag()
					+ b[i][2] * psi_cp_im[i];
		}
	}

	std::complex<double> psi_total()
	{
		double psi_re = std::accumulate(psi_dp_re.begin(), psi_dp_re.end(), double(0))
				+ std::accumulate(psi_cp_re.begin(), psi_cp_re.end(), std::complex<double>(0)).real();
		double psi_im = std::accumulate(psi_dp_im.begin(), psi_dp_im.end(), double(0))
				+ std::accumulate(psi_cp_im.begin(), psi_cp_im.end(), std::complex<double>(0)).real();
		return std::complex<double>(psi_re, psi_im);
	}

protected:
	double epsilon;
	std::vector<std::vector<double> > a;
	std::vector<std::vector<std::complex<double> > > b;
	std::vector<double> c;
	// *_re and *_im are for the real and imaginary part of the e-field, respectively.
	std::vector<double> psi_dp_re, psi_dp_im;
	std::vector<std::complex<double> > psi_cp_re, psi_cp_im;
};

template <typename S, typename T>
static inline T& assign(const std::complex<S>& in, T& out)
{
	return out = static_cast<T>(in.real());
}

template <typename S, typename T>
static inline std::complex<T>& assign(const std::complex<S>& in, std::complex<T>& out)
{
	out.real() = static_cast<T>(in.real());
	out.imag() = static_cast<T>(in.imag());
	return out;
}

template <typename T> class DCPPLRCEx: public DCPPLRCElectric<T>
{
public:
	DCPPLRCEx(const int * const idx, int size, double epsilon,
			const double * const a, int a_i_size, int a_j_size,
			const std::complex<double> * const b_c, int b_c_i_size, int b_c_j_size,
			const double * const c, int c_size) :
				DCPPLRCElectric<T>(idx, size, epsilon, a, a_i_size, a_j_size, b_c, b_c_i_size, b_c_j_size, c, c_size)
	{
	}

	void update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			double dy, double dz, double dt, double n)
	{
		std::complex<double> e_now = ex(i,j,k);
		std::complex<double> e_new = c[0] * e_now
				+ c[1] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz)
				+ c[2] * psi_total();

		update_psi_dp(e_now, e_new);
		update_psi_cp(e_now, e_new);
		assign(e_new, ex(i,j,k));
	}

protected:
	using DCPPLRCElectric<T>::epsilon;
	using DCPPLRCElectric<T>::a;
	using DCPPLRCElectric<T>::b;
	using DCPPLRCElectric<T>::c;
	using DCPPLRCElectric<T>::psi_dp_re;
	using DCPPLRCElectric<T>::psi_dp_im;
	using DCPPLRCElectric<T>::psi_cp_re;
	using DCPPLRCElectric<T>::psi_cp_im;
	using DCPPLRCElectric<T>::update_psi_dp;
	using DCPPLRCElectric<T>::update_psi_cp;
	using DCPPLRCElectric<T>::psi_total;
};

template <typename T> class DCPPLRCEy: public DCPPLRCElectric<T>
{
public:
	DCPPLRCEy(const int * const idx, int size, double epsilon,
			const double * const a, int a_i_size, int a_j_size,
			const std::complex<double> * const b_c, int b_c_i_size, int b_c_j_size,
			const double * const c, int c_size) :
				DCPPLRCElectric<T>(idx, size, epsilon, a, a_i_size, a_j_size, b_c, b_c_i_size, b_c_j_size, c, c_size)
	{
	}

	void update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			double dz, double dx, double dt, double n)
	{
		std::complex<double> e_now = ey(i,j,k);
		std::complex<double> e_new = c[0] * e_now
				+ c[1] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx)
				+ c[2] * psi_total();

		update_psi_dp(e_now, e_new);
		update_psi_cp(e_now, e_new);

		assign(e_new, ey(i,j,k));
	}

protected:
	using DCPPLRCElectric<T>::epsilon;
	using DCPPLRCElectric<T>::a;
	using DCPPLRCElectric<T>::b;
	using DCPPLRCElectric<T>::c;
	using DCPPLRCElectric<T>::psi_dp_re;
	using DCPPLRCElectric<T>::psi_dp_im;
	using DCPPLRCElectric<T>::psi_cp_re;
	using DCPPLRCElectric<T>::psi_cp_im;
	using DCPPLRCElectric<T>::update_psi_dp;
	using DCPPLRCElectric<T>::update_psi_cp;
	using DCPPLRCElectric<T>::psi_total;
};

template <typename T> class DCPPLRCEz: public DCPPLRCElectric<T>
{
public:
	DCPPLRCEz(const int * const idx, int size, double epsilon,
			const double * const a, int a_i_size, int a_j_size,
			const std::complex<double> * const b_c, int b_c_i_size, int b_c_j_size,
			const double * const c, int c_size) :
				DCPPLRCElectric<T>(idx, size, epsilon, a, a_i_size, a_j_size, b_c, b_c_i_size, b_c_j_size, c, c_size)
	{
	}

	void update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			double dx, double dy, double dt, double n)
	{
		std::complex<double> e_now = ez(i,j,k);
		std::complex<double> e_new = c[0] * e_now
				+ c[1] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy)
				+ c[2] * psi_total();

		update_psi_dp(e_now, e_new);
		update_psi_cp(e_now, e_new);

		assign(e_new, ez(i,j,k));
	}

protected:
	using DCPPLRCElectric<T>::epsilon;
	using DCPPLRCElectric<T>::a;
	using DCPPLRCElectric<T>::b;
	using DCPPLRCElectric<T>::c;
	using DCPPLRCElectric<T>::psi_dp_re;
	using DCPPLRCElectric<T>::psi_dp_im;
	using DCPPLRCElectric<T>::psi_cp_re;
	using DCPPLRCElectric<T>::psi_cp_im;
	using DCPPLRCElectric<T>::update_psi_dp;
	using DCPPLRCElectric<T>::update_psi_cp;
	using DCPPLRCElectric<T>::psi_total;
};

template <typename T> class DCPPLRCHx: public DielectricHx<T>
{
public:
	DCPPLRCHx(const int * const idx, int size, double mu = 1) :
		DielectricHx<T>(idx, size, mu)
	{
	}
};

template <typename T> class DCPPLRCHy: public DielectricHy<T>
{
public:
	DCPPLRCHy(const int * const idx, int size, double mu = 1) :
		DielectricHy<T>(idx, size, mu)
	{
	}
};

template <typename T> class DCPPLRCHz: public DielectricHz<T>
{
public:
	DCPPLRCHz(const int * const idx, int size, double mu = 1) :
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
