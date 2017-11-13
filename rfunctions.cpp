
#include <stdlib.h>
#include <math.h>
#include "rfunctions.hpp"
#include "hypergeometric.hpp"

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                            R FUNCTIONS                           ::
//                                                                  ::
//      R_nlq functions that appear in the anisotropic integrands   ::
//      Ea, PTa, PLa, I402-1, I421-1, respectively                  ::
//                                                                  ::
//             R200     R201     R220     R402      R421            ::
//                                                                  ::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


//Ea, I2001
double R200(double pbar, double ax, double az, double mbar)
{
	double mop2 = (mbar*mbar) / (pbar*pbar);
	double x = (ax*ax - az*az) / (az*az + mop2);
	return tfuncE(x) * sqrt(az*az + mop2);
}


//PTa, I2011
double R201(double pbar, double ax, double az, double mbar)
{
	double mop2 = (mbar*mbar) / (pbar*pbar);
	double x = (ax*ax - az*az) / (az*az + mop2);
	return tfuncPT(x) / sqrt(az*az + mop2);
}


//PLa, I2201
double R220(double pbar, double ax, double az, double mbar)
{
	double mop2 = (mbar*mbar) / (pbar*pbar);
	double x = (ax*ax - az*az) / (az*az + mop2);
	return tfuncPL(x) / sqrt(az*az + mop2);
}


//I402-1
double R402(double pbar, double ax, double az, double mbar)
{
	double mop2 = (mbar*mbar) / (pbar*pbar);
	double x = (ax*ax - az*az) / (az*az + mop2);
	return tfunc402(x) / sqrt(az*az + mop2);
}


//I421-1
double R421(double pbar, double ax, double az, double mbar)
{
	double mop2 = (mbar*mbar) / (pbar*pbar);
	double x = (ax*ax - az*az) / (az*az + mop2);
	return tfunc421(x) / sqrt(az*az + mop2);
}


//I401-1
double R401(double pbar, double ax, double az, double mbar)
{
	double mop2 = (mbar*mbar) / (pbar*pbar);
	double x = (ax*ax - az*az) / (az*az + mop2);
	return tfunc401(x) * sqrt(az*az + mop2);
}


//I420-1
double R420(double pbar, double ax, double az, double mbar)
{
	double mop2 = (mbar*mbar) / (pbar*pbar);
	double x = (ax*ax - az*az) / (az*az + mop2);
	return tfunc420(x) * sqrt(az*az + mop2);
}


//I440-1
double R440(double pbar, double ax, double az, double mbar)
{
	double mop2 = (mbar*mbar) / (pbar*pbar);
	double x = (ax*ax - az*az) / (az*az + mop2);
	return tfunc440(x) / sqrt(az*az + mop2);
}








