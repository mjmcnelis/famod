
#include <stdlib.h>
#include <math.h>
#include "anisotropic_integrands.hpp"
#include "rfunctions.hpp"


// thermal energy density
double Eeq_integrand(double pbar, double mbar, int sign)
{
	// gauss laguerre (a = 2)
	return sqrt(pbar*pbar + mbar*mbar) * exp(pbar) / (exp(sqrt(pbar*pbar + mbar*mbar)) + (double)sign);
}

// thermal pressure
double Peq_integrand(double pbar, double mbar, int sign)
{
	// gauss laguerre (a = 2)
	return pbar*pbar / sqrt(pbar*pbar + mbar*mbar) * exp(pbar) / (exp(sqrt(pbar*pbar + mbar*mbar)) + (double)sign);
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                1D ANISOTROPIC FUNCTION INTEGRANDS                ::
//                                                                  ::
//     Momentum bar integrand of the anisotropic functions.         ::
//	   (written in gla form) Integrate with gauss_aniso_1D.         ::
//                                                                  ::
//     Ea_integrand        PTa_integrand        PLa_integrand       ::
//     I402m1_integrand    I421m1_integrand                         ::
//																	::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


double Ea_integrand(double pbar, double ax, double az, double mbar, int sign)
{
	// gauss laguerre (a = 2)
	return pbar * R200(pbar,ax,az,mbar) * exp(pbar) / (exp(sqrt(pbar*pbar + mbar*mbar)) + (double)sign);
}


double PTa_integrand(double pbar, double ax, double az, double mbar, int sign)
{
	// gauss laguerre (a = 2)
	return pbar * R201(pbar,ax,az,mbar) * exp(pbar) / (exp(sqrt(pbar*pbar + mbar*mbar)) + (double)sign);
}


double PLa_integrand(double pbar, double ax, double az, double mbar, int sign)
{
	// gauss laguerre (a = 2)
	return pbar * R220(pbar,ax,az,mbar) * exp(pbar) / (exp(sqrt(pbar*pbar + mbar*mbar)) + (double)sign);
}


double I402m1_integrand(double pbar, double ax, double az, double mbar, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Eabar) + (double)sign;
	// fa*(1-sign*fa) = exp(Eabar)/(exp(Eabar)+sign)^2
	// gauss laguerre (a = 2)
	return (pbar*pbar / Eabar) * R402(pbar,ax,az,mbar) * exp(pbar+Eabar) / (qstat*qstat);
}


double I421m1_integrand(double pbar, double ax, double az, double mbar, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Eabar) + (double)sign;
	// gauss laguerre (a = 3)
	return (pbar*pbar / Eabar) * R421(pbar,ax,az,mbar) * exp(pbar+Eabar) / (qstat*qstat);
}



//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                TEST 3D ANISOTROPIC FUNCTION INTEGRANDS           ::
//                                                                  ::
//      Test functions to show gauss_test_aniso_3D works.           ::
//	    Compare with the corresponding functions integrated with    ::
//      gauss_aniso_1D and gauss_mod_aniso_3D w/ no modification    ::
//                                                                  ::
//         test3D_Ea         test3D_PTa         test3D_PLa          ::
//																	::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


double test3D_Ea_integrand(double xphi, double costheta, double pbar, double ax, double az, double mbar, int sign)
{
	double costheta2 = costheta * costheta;
	// gauss laguerre (a = 2)
	return pbar * sqrt(ax*ax*(1.0-costheta2) + az*az*costheta2 + (mbar*mbar)/(pbar*pbar)) * exp(pbar) / (exp(sqrt(pbar*pbar + mbar*mbar)));
}


double test3D_PTa_integrand(double xphi, double costheta, double pbar, double ax, double az, double mbar, int sign)
{
	double costheta2 = costheta * costheta;
	double sintheta2 = 1.0 - costheta2;
	// gauss laguerre (a = 2)
	return pbar * sintheta2 / sqrt(ax*ax*sintheta2 + az*az*costheta2 + (mbar*mbar)/(pbar*pbar)) * exp(pbar) / (exp(sqrt(pbar*pbar + mbar*mbar)) + (double)sign);
}


double test3D_PLa_integrand(double xphi, double costheta, double pbar, double ax, double az, double mbar, int sign)
{
	double costheta2 = costheta * costheta;
	double sintheta2 = 1.0 - costheta2;
	// gauss laguerre (a = 2)
	return pbar * costheta2 / sqrt(ax*ax*sintheta2 + az*az*costheta2 + (mbar*mbar)/(pbar*pbar)) * exp(pbar) / (exp(sqrt(pbar*pbar + mbar*mbar)) + (double)sign);
}



//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//              MODDED 3D ANISOTROPIC FUNCTION INTEGRANDS             ::
//                                                                    ::
//      Hydrodynamic moments of the modded anisotropic distribution   ::
//      function. Integrate with gauss_mod_aniso_3D. F                ::
//                                                                    ::
//						 (10 components of T^munu)                    ::
//                                                                    ::
//              mod_Ea            mod_PTa           mod_PLa           ::
//              mod_pixx          mod_pixy          mod_Wxz           ::
//              mod_Wyz           mod_Ttx           mod_Tty           ::
//              mod_Ttz                                               ::
// 																	  ::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::




double modEa_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign)
{
	// angles of pstar, pbar -> radial momentum of pstar_bar
	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta * costheta); // sintheta: right domain

	// coordinate transformation p = A * B * pstar
	// A = vahydro momentum matrix (see main)
	// B = diag(ax,ax,az) is the ahydro momentum matrix
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta}; // unit vector of pstar
	// B * unit vector of pstar
	unit_vec[0] = ax * unit_vec[0];
	unit_vec[1] = ax * unit_vec[1];
	unit_vec[2] = az * unit_vec[2];

	// A * (B * unit_vec)
	for(int i = 0; i < n; i++) for(int k = 0; k < n; k++) p_unit[i] += A[i][k] * unit_vec[k];
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	// E(p(pstar)) / |pstar| from coordinate transformation
	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	// gauss laguerre (a = 2)
	return pbar * energy_unit * exp(pbar) / (exp(sqrt(mbar*mbar + pbar*pbar)) + (double)sign);
}


double modPTa_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign)
{
	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta * costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	unit_vec[0] = ax * unit_vec[0];
	unit_vec[1] = ax * unit_vec[1];
	unit_vec[2] = az * unit_vec[2];

	for(int i = 0; i < n; i++) for(int k = 0; k < n; k++) p_unit[i] += A[i][k] * unit_vec[k];

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	// p_i(pstar) / |pstar| vector
	double px_unit = p_unit[0];
	double py_unit = p_unit[1];

	// gauss laguerre (a = 2)
	return pbar * (px_unit * px_unit + py_unit * py_unit) / energy_unit * exp(pbar) / (exp(sqrt(mbar*mbar + pbar*pbar)) + (double)sign);


}


double modPLa_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign)
{
	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta * costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	unit_vec[0] = ax * unit_vec[0];
	unit_vec[1] = ax * unit_vec[1];
	unit_vec[2] = az * unit_vec[2];

	for(int i = 0; i < n; i++) for(int k = 0; k < n; k++) p_unit[i] += A[i][k] * unit_vec[k];

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double pz_unit = p_unit[2];

	// gauss laguerre (a = 2)
	return pbar * pz_unit * pz_unit / energy_unit * exp(pbar) / (exp(sqrt(mbar*mbar + pbar*pbar)) + (double)sign);
}


double modpixx_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign)
{
	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta * costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	unit_vec[0] = ax * unit_vec[0];
	unit_vec[1] = ax * unit_vec[1];
	unit_vec[2] = az * unit_vec[2];

	for(int i = 0; i < n; i++) for(int k = 0; k < n; k++) p_unit[i] += A[i][k] * unit_vec[k];

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double px_unit = p_unit[0];
	double py_unit = p_unit[1];

	// gauss laguerre (a = 2)
	return pbar * (px_unit * px_unit - py_unit * py_unit) / energy_unit * exp(pbar) / (exp(sqrt(mbar*mbar + pbar*pbar)) + (double)sign);
}


double modpixy_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign)
{
	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta * costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	unit_vec[0] = ax * unit_vec[0];
	unit_vec[1] = ax * unit_vec[1];
	unit_vec[2] = az * unit_vec[2];

	for(int i = 0; i < n; i++) for(int k = 0; k < n; k++) p_unit[i] += A[i][k] * unit_vec[k];

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double px_unit = p_unit[0];
	double py_unit = p_unit[1];

	// gauss laguerre (a = 2)
	return pbar * px_unit * py_unit / energy_unit * exp(pbar) / (exp(sqrt(mbar*mbar + pbar*pbar)) + (double)sign);
}


double modWxz_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign)
{
	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta * costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	unit_vec[0] = ax * unit_vec[0];
	unit_vec[1] = ax * unit_vec[1];
	unit_vec[2] = az * unit_vec[2];

	for(int i = 0; i < n; i++) for(int k = 0; k < n; k++) p_unit[i] += A[i][k] * unit_vec[k];

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double px_unit = p_unit[0];
	double pz_unit = p_unit[2];

	// gauss laguerre (a = 2)
	return pbar * px_unit * pz_unit / energy_unit * exp(pbar) / (exp(sqrt(mbar*mbar + pbar*pbar)) + (double)sign);
}


double modWyz_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign)
{
	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta * costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	unit_vec[0] = ax * unit_vec[0];
	unit_vec[1] = ax * unit_vec[1];
	unit_vec[2] = az * unit_vec[2];

	for(int i = 0; i < n; i++) for(int k = 0; k < n; k++) p_unit[i] += A[i][k] * unit_vec[k];

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double py_unit = p_unit[1];
	double pz_unit = p_unit[2];

	// gauss laguerre (a = 2)
	return pbar * py_unit * pz_unit / energy_unit * exp(pbar) / (exp(sqrt(mbar*mbar + pbar*pbar)) + (double)sign);
}


double modTtx_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign)
{
	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta * costheta);

	double px_unit = 0.0;
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	unit_vec[0] = ax * unit_vec[0];
	unit_vec[1] = ax * unit_vec[1];
	unit_vec[2] = az * unit_vec[2];

	for(int i = 0; i < n; i++) for(int k = 0; k < n; k++) px_unit += A[0][k] * unit_vec[k];

	// gauss laguerre (a = 2)
	return pbar * px_unit * exp(pbar) / (exp(sqrt(mbar*mbar + pbar*pbar)) + (double)sign);
}


double modTty_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign)
{
	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta * costheta);

	double py_unit = 0.0;
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	unit_vec[0] = ax * unit_vec[0];
	unit_vec[1] = ax * unit_vec[1];
	unit_vec[2] = az * unit_vec[2];

	for(int i = 0; i < n; i++) for(int k = 0; k < n; k++) py_unit += A[1][k] * unit_vec[k];

	// gauss laguerre (a = 2)
	return pbar * py_unit * exp(pbar) / (exp(sqrt(mbar*mbar + pbar*pbar)) + (double)sign);
}


double modTtz_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign)
{
	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta * costheta);

	double pz_unit = 0.0;
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	unit_vec[0] = ax * unit_vec[0];
	unit_vec[1] = ax * unit_vec[1];
	unit_vec[2] = az * unit_vec[2];

	for(int i = 0; i < n; i++) for(int k = 0; k < n; k++) pz_unit += A[2][k] * unit_vec[k];

	// gauss laguerre (a = 2)
	return pbar * pz_unit * exp(pbar) / (exp(sqrt(mbar*mbar + pbar*pbar)) + (double)sign);
}










// additional anisotropic integrands that appear
// in jacobian of root finding algoritm
/////////////////////////////////////////////////////////////////////////


double I2001_integrand(double pbar, double ax, double az, double mbar, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Eabar) + (double)sign;
	// gauss laguerre (a = 3)
	return Eabar * R200(pbar,ax,az,mbar) * exp(pbar+Eabar) / (qstat*qstat);
}


double I2011_integrand(double pbar, double ax, double az, double mbar, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Eabar) + (double)sign;
	// gauss laguerre (a = 3)
	return Eabar * R201(pbar,ax,az,mbar) * exp(pbar+Eabar) / (qstat*qstat);
}


double I2201_integrand(double pbar, double ax, double az, double mbar, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Eabar) + (double)sign;
	// gauss laguerre (a = 3)
	return Eabar * R220(pbar,ax,az,mbar) * exp(pbar+Eabar) / (qstat*qstat);
}


double I401m1_integrand(double pbar, double ax, double az, double mbar, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Eabar) + (double)sign;
	// gauss laguerre (a = 3)
	return (pbar*pbar / Eabar) * R401(pbar,ax,az,mbar) * exp(pbar+Eabar) / (qstat*qstat);
}


double I420m1_integrand(double pbar, double ax, double az, double mbar, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Eabar) + (double)sign;
	// gauss laguerre (a = 3)
	return (pbar*pbar / Eabar) * R420(pbar,ax,az,mbar) * exp(pbar+Eabar) / (qstat*qstat);
}


double I440m1_integrand(double pbar, double ax, double az, double mbar, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Eabar) + (double)sign;
	// gauss laguerre (a = 3)
	return (pbar*pbar / Eabar) * R440(pbar,ax,az,mbar) * exp(pbar+Eabar) / (qstat*qstat);
}
















































