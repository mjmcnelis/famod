
#include <stdlib.h>
#include <math.h>
#include "anisotropic_integrands.hpp"
#include "rfunctions.hpp"


// equilibrium net baryon density
double nBeq_integrand(double pbar, double mbar, double aB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	// gla (a = 1)
	return bn * pbar * exp(pbar) / (exp(Ebar-bn*aB)+a);
}

// equilibrium energy density
double Eeq_integrand(double pbar, double mbar, double aB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	// gla (a = 2)
	return Ebar * exp(pbar) / (exp(Ebar-bn*aB)+a);
}

// equilibrium pressure
double Peq_integrand(double pbar, double mbar, double aB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	// gla (a = 2)
	return pbar*pbar/Ebar * exp(pbar) / (exp(Ebar-bn*aB)+a);
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

double nBa_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	// gla (a = 1)
	return bn * pbar * exp(pbar) / (exp(Eabar-bn*aBt)+a);
}

double Ea_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	// gla (a = 2)
	return pbar * R200(pbar,ax,az,mbar) * exp(pbar) / (exp(Eabar-bn*aBt)+a);
}


double PTa_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	// gla (a = 2)
	return pbar * R201(pbar,ax,az,mbar) * exp(pbar) / (exp(Eabar-bn*aBt)+a);
}


double PLa_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	// gla (a = 2)
	return pbar * R220(pbar,ax,az,mbar) * exp(pbar) / (exp(Eabar-bn*aBt)+a);
}

double bn_I1001_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// gla (a = 2)
	return bn * Eabar * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
}

double bn2_I1000_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// gla (a = 1)
	return bn * bn * pbar * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
}


double bn_I2000_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// gla (a = 2)
	return bn * pbar * R200(pbar,ax,az,mbar) * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
}

double bn_I2010_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// gla (a = 2)
	return bn * pbar * R201(pbar,ax,az,mbar) * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
}


double bn_I2200_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// gla (a = 2)
	return bn * pbar * R220(pbar,ax,az,mbar) * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
}



double bn_I301m1_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// gla (a = 2)
	return bn * pbar * pbar / Eabar * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
}

double bn_I320m1_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// gla (a = 2)
	return bn * pbar * pbar / Eabar * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
}

double I402m1_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// fa*(1-sign*fa) = exp(Eabar)/(exp(Eabar)+sign)^2
	// gla (a = 2)
	return (pbar*pbar / Eabar) * R402(pbar,ax,az,mbar) * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
}


double I421m1_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// gla (a = 3)
	return (pbar*pbar / Eabar) * R421(pbar,ax,az,mbar) * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
}


double I2001_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// gla (a = 3)
	return Eabar * R200(pbar,ax,az,mbar) * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
}


double I2011_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// gla (a = 3)
	return Eabar * R201(pbar,ax,az,mbar) * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
}


double I2201_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// gla (a = 3)
	return Eabar * R220(pbar,ax,az,mbar) * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
}


double I401m1_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// gla (a = 3)
	return (pbar*pbar / Eabar) * R401(pbar,ax,az,mbar) * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
}


double I420m1_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// gla (a = 3)
	return (pbar*pbar / Eabar) * R420(pbar,ax,az,mbar) * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
}


double I440m1_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Eabar-bn*aBt)+a;
	// gla (a = 3)
	return (pbar*pbar / Eabar) * R440(pbar,ax,az,mbar) * exp(pbar+Eabar-bn*aBt)/(qstat*qstat);
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




double modEa_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
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

	// gla (a = 2)
	return pbar * energy_unit * exp(pbar) / (exp(Eabar-bn*aBt)+a);
}


double modPTa_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

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

	// gla (a = 2)
	return pbar * (px_unit * px_unit + py_unit * py_unit) / energy_unit * exp(pbar) / (exp(Eabar-bn*aBt)+a);


}


double modPLa_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

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

	// gla (a = 2)
	return pbar * pz_unit * pz_unit / energy_unit * exp(pbar) / (exp(Eabar-bn*aBt)+a);
}


double modpixx_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

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

	// gla (a = 2)
	return pbar * (px_unit * px_unit - py_unit * py_unit) / energy_unit * exp(pbar) / (exp(Eabar-bn*aBt)+a);
}


double modpixy_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

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

	// gla (a = 2)
	return pbar * px_unit * py_unit / energy_unit * exp(pbar) / (exp(Eabar-bn*aBt)+a);
}


double modWxz_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

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

	// gla (a = 2)
	return pbar * px_unit * pz_unit / energy_unit * exp(pbar) / (exp(Eabar-bn*aBt)+a);
}


double modWyz_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

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

	// gla (a = 2)
	return pbar * py_unit * pz_unit / energy_unit * exp(pbar) / (exp(Eabar-bn*aBt)+a);
}


double modTtx_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta * costheta);

	double px_unit = 0.0;
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	unit_vec[0] = ax * unit_vec[0];
	unit_vec[1] = ax * unit_vec[1];
	unit_vec[2] = az * unit_vec[2];

	for(int k = 0; k < n; k++) px_unit += A[0][k] * unit_vec[k];

	// gla (a = 2)
	return pbar * px_unit * exp(pbar) / (exp(Eabar-bn*aBt)+a);
}


double modTty_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta * costheta);

	double py_unit = 0.0;
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	unit_vec[0] = ax * unit_vec[0];
	unit_vec[1] = ax * unit_vec[1];
	unit_vec[2] = az * unit_vec[2];

	for(int k = 0; k < n; k++) py_unit += A[1][k] * unit_vec[k];

	// gla (a = 2)
	return pbar * py_unit * exp(pbar) / (exp(Eabar-bn*aBt)+a);
}


double modTtz_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign)
{
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta * costheta);

	double pz_unit = 0.0;
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	unit_vec[0] = ax * unit_vec[0];
	unit_vec[1] = ax * unit_vec[1];
	unit_vec[2] = az * unit_vec[2];

	for(int k = 0; k < n; k++) pz_unit += A[2][k] * unit_vec[k];

	// gla (a = 2)
	return pbar * pz_unit * exp(pbar) / (exp(Eabar-bn*aBt)+a);
}

















































