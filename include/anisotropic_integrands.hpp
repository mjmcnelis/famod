
#include <stdlib.h>

#ifndef ANISOTROPIC_INTEGRANDS_H

#define ANISOTROPIC_INTEGRANDS_H

double Eeq_integrand(double pbar, double mbar, int sign);
double Peq_integrand(double pbar, double mbar, int sign);

double Ea_integrand(double pbar, double ax, double az, double mbar, int sign);
double PTa_integrand(double pbar, double ax, double az, double mbar, int sign);
double PLa_integrand(double pbar, double ax, double az, double mbar, int sign);

// for the anisotropic chapman-enskog approximation and jacobian in root finder
double I402m1_integrand(double pbar, double ax, double az, double mbar, int sign);
double I421m1_integrand(double pbar, double ax, double az, double mbar, int sign);

// remaining anisotropic integrands for jacobian in root finder
double I2001_integrand(double pbar, double ax, double az, double mbar, int sign);
double I2011_integrand(double pbar, double ax, double az, double mbar, int sign);
double I2201_integrand(double pbar, double ax, double az, double mbar, int sign);
double I401m1_integrand(double pbar, double ax, double az, double mbar, int sign);
double I420m1_integrand(double pbar, double ax, double az, double mbar, int sign);
double I440m1_integrand(double pbar, double ax, double az, double mbar, int sign);


double test3D_Ea_integrand(double xphi, double costheta, double pbar, double ax, double az, double mbar, int sign);
double test3D_PTa_integrand(double xphi, double costheta, double pbar, double ax, double az, double mbar, int sign);
double test3D_PLa_integrand(double xphi, double costheta, double pbar, double ax, double az, double mbar, int sign);


double modEa_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign);
double modPTa_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign);
double modPLa_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign);
double modpixx_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign);
double modpixy_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign);
double modWxz_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign);
double modWyz_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign);
double modTtx_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign);
double modTty_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign);
double modTtz_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign);


#endif