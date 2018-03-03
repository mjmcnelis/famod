
#include <stdlib.h>

#ifndef ANISOTROPIC_INTEGRANDS_H

#define ANISOTROPIC_INTEGRANDS_H

double Eeq_integrand(double pbar, double mbar, double aB, int baryon, int sign);
double Peq_integrand(double pbar, double mbar, double aB, int baryon, int sign);
double nBeq_integrand(double pbar, double mbar, double aB, int baryon, int sign);

double Ea_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);
double PTa_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);
double PLa_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);
double nBa_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);

// for the anisotropic chapman-enskog approximation and jacobian in root finder
double I402m1_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);
double I421m1_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);

// remaining anisotropic integrands for jacobian in root finder
double I2001_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);
double I2011_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);
double I2201_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);
double I401m1_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);
double I420m1_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);
double I440m1_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);

double bn_I1001_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);
double bn2_I1000_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);
double bn_I2000_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);
double bn_I2010_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);
double bn_I2200_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);

double bn_I301m1_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);
double bn_I320m1_integrand(double pbar, double ax, double az, double mbar, double aBt, int baryon, int sign);

double modEa_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign);
double modPTa_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign);
double modPLa_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign);
double modpixx_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign);
double modpixy_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign);
double modWxz_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign);
double modWyz_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign);
double modTtx_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign);
double modTty_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign);
double modTtz_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, double aBt, int baryon, int sign);


#endif