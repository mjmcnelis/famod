
#include <stdlib.h>
#include "gauss_integration.hpp"

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                    1D and 3D GAUSS INTEGRATION                     ::
//                                                                    ::
//     Compute 1D anisotropic integrals over radial momentum          ::
//     bar using Gauss Laguerre quadrature. Compute 3D test           ::
//	   and modded anisotropic integrals. The two angular integrals    ::
//     are computed with Clenshaw Curtis quadrature and the radial    ::
//     momentum bar with Gauss Laguerre quadrature.                   ::
//                                                                    ::
//     gauss_aniso_1D    gauss_test_aniso_3D    gauss_mod_aniso_3D    ::
//                                                                    ::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



// for thermal equilibrium functions
double Gauss_Thermo_1D(double thermal_1D_integrand(double pbar, double mbar, int sign), double * pbar_root, double * pbar_weight, int pbar_pts, double mbar, int sign)
{
	double sum = 0.0;
	for(int k = 0; k < pbar_pts; k++) sum += pbar_weight[k] * thermal_1D_integrand(pbar_root[k], mbar, sign);
	return sum;
}



double Gauss_Aniso_1D(double aniso_1D_integrand(double pbar, double ax, double az, double mbar, int sign), double * pbar_root, double * pbar_weight, int pbar_pts, double ax, double az, double mbar, int sign)
{
	double sum = 0.0;
	for(int k = 0; k < pbar_pts; k++) sum += pbar_weight[k] * aniso_1D_integrand(pbar_root[k], ax, az, mbar, sign);
	return sum;
}



// only for testing 3D integration: precursor for the mod version
double gauss_test_aniso_3D(double test_aniso_3D_integrand(double xphi, double costheta, double pbar, double ax, double az, double mbar, int sign), double * xphi_root, double * xphi_weight, double * costheta_root, double * costheta_weight, double * pbar_root, double * pbar_weight, int xphi_pts, int costheta_pts, int pbar_pts, double ax, double az, double mbar, int sign)
{
	double sum = 0.0;
	for(int i = 0; i < xphi_pts; i++)
	{
		for(int j = 0; j < costheta_pts; j++)
		{
			for(int k = 0; k < pbar_pts; k++)
			{
				sum += xphi_weight[i] * costheta_weight[j] * pbar_weight[k] * test_aniso_3D_integrand(xphi_root[i], costheta_root[j], pbar_root[k], ax, az, mbar, sign);
			}
		}
	}
	return sum;
}


double Gauss_Mod_Aniso_3D(double mod_aniso_3D_integrand(double xphi, double costheta, double pbar, double ax, double az, double ** A, int n, double mbar, int sign), double * xphi_root, double * xphi_weight, double * costheta_root, double * costheta_weight, double * pbar_root, double * pbar_weight, int xphi_pts, int costheta_pts, int pbar_pts, double ax, double az, double ** A, int n, double mbar, int sign)
{
	double sum = 0.0;
	for(int i = 0; i < xphi_pts; i++)
	{
		for(int j = 0; j < costheta_pts; j++)
		{
			for(int k = 0; k < pbar_pts; k++)
			{
				sum += xphi_weight[i] * costheta_weight[j] * pbar_weight[k] * mod_aniso_3D_integrand(xphi_root[i], costheta_root[j], pbar_root[k], ax, az, A, n, mbar, sign);
			}
		}
	}
	return sum;
}





















