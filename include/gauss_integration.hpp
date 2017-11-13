
#include <stdlib.h>

#ifndef GAUSS_INTEGRATION_H

#define GAUSS_INTEGRATION_H

double Gauss_Thermo_1D(double thermal_1D_integrand(double pbar, double mbar, int sign), double * pbar_root, double * pbar_weight, int pbar_pts, double mbar, int sign);

double Gauss_Aniso_1D(double aniso_1D_integrand(double pbar, double ax, double az, double mbar, int sign), double * pbar_root, double * pbar_weight, int pbar_pts, double ax, double az, double mbar, int sign);

double gauss_test_aniso_3D(double test_aniso_3D_integrand(double xphi, double costheta, double pbar, double ax, double az, double mbar, int sign), double * xphi_root, double * xphi_weight, double * costheta_root, double * costheta_weight, double * pbar_root, double * pbar_weight, int xphi_pts, int costheta_pts, int pbar_pts, double ax, double az, double mbar, int sign);

double Gauss_Mod_Aniso_3D(double mod_aniso_3D_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign), double * xphi_root, double * xphi_weight, double * costheta_root, double * costheta_weight, double * pbar_root, double * pbar_weight, int xphi_pts, int costheta_pts, int pbar_pts, double ax, double az, double ** A, int n, double mbar, int sign);


#endif
