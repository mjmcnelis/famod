
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "anisotropicvariables.hpp"
#include "anisotropic_integrands.hpp"
#include "gauss_integration.hpp"
#define EPS_MIN 1.0e-16


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                    LUP LINEAR EQUATION SOLVER                    ::
//                                                                  ::
//     Solves linear equation Ax = b using LU decomposition         ::
//	   with implicit partial pivoting (LUP). To directly solve      ::
//     for Ax = b, run these two functions concurrently:            ::
//                                                                  ::
//            LUP_decomposition              LUP_solve              ::
//                                                                  ::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


void LUP_decomposition(double ** A, int n, int * pvector)
{
	int i;
	int j;
	int k;
	int imax;
	double big;
	double sum;
	double temp;
	double implicit_scale[n];
	for(i = 0; i < n; i++)
	{
		pvector[i] = i;
	}
	for(i = 0; i < n; i++)
	{
		big = 0.0;
		for(j = 0; j < n; j++)
		{
			temp = fabs(A[i][j]);
			if(temp > big)
			{
				big = temp;
			}
		}
		if(big == 0.0)
		{
			printf("Singular matrix in the routine");
			break;
		}
		implicit_scale[i] = 1.0 / big;
	}
	// LU Decomposition
	for(j = 0; j < n; j++)
	{
		for(i = 0; i < j; i++)
		{
			sum = A[i][j];
			for(k = 0; k < i; k++)
			{
				sum -= A[i][k] * A[k][j];
			}
			A[i][j] = sum;
		}
		big = 0.0;
		for(i = j; i < n; i++)
		{
			sum = A[i][j];
			for(k = 0; k < j; k++)
			{
				sum -= A[i][k] * A[k][j];
			}
			A[i][j] = sum;

			temp = implicit_scale[i] * fabs(sum);
			if(temp >= big)
			{
				big = temp;
				imax = i;
			}
		}
		if(j != imax)
		{
			for(k = 0; k < n; k++)
				{
					temp = A[imax][k];
					A[imax][k] = A[j][k];
					A[j][k] = temp;
				}
			implicit_scale[imax] = implicit_scale[j];
		}

		pvector[j] = imax;

		if(A[j][j] == 0.0)
		{
			A[j][j] = EPS_MIN;
		}
		if(j != n-1)
		{
			temp = 1.0 / A[j][j];
			for(i = j+1; i < n; i++)
			{
				A[i][j] *= temp;
			}
		}
	}
}

void LUP_solve(double ** PA, int n, int * pvector, double b[])
{
	int i;
	int j;
	int m = -1;
	int ip;
	double sum;
	// Forward substitution routine for Ly = b
	for(i = 0; i < n; i++)
	{
		ip = pvector[i];
		sum = b[ip];
		b[ip] = b[i];
		if(m != -1)
		{
			for(j = m; j <= i-1; j++)
			{
				sum -= PA[i][j] * b[j];
			}
		}
		else if(sum)
		{
			m = i;
		}
		b[i] = sum;
	}
	// Backward substitution routine for Ux = y
	for(i = n-1; i >= 0; i--)
	{
		sum = b[i];
		for(j = i+1; j < n; j++)
		{
			sum -= PA[i][j] * b[j];
		}
		b[i] = sum / PA[i][i];
	}
}





//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                       ANISOTROPIC VARIABLES                      ::
//                                                                  ::
//     Compute the anisotropic variables (lambda,ax,az) from        ::
//	   the quantities (e_kinetic,pl,pt)       		                ::
//                                                                  ::
//                     find_anisotropic_variables                   ::
//																	::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


void find_anisotropic_variables(double e, double pl, double pt, double nB, double * mass, int * degeneracy, int * baryon, int * sign, int number_resonances, double * pbar_rootN, double * pbar_weightN, double * pbar_rootT, double * pbar_weightT, double * pbar_rootJ, double * pbar_weightJ, int pbar_pts, double * lambda, double * ax, double * az, double * aBt)
{
	const double Ea = e;
	const double PTa = pt;
	const double PLa = pl;
	const double nBa = nB;  // make it number 4

	double dof;

	double lambdai = *lambda;			                     // initial guess for anisotropic variables
	double axi = *ax;
	double azi = *az;
	double aBti = *aBt;

	const int n = 4; 										 // dimension space
	double X[n] = {lambdai, axi, azi, aBti};			     // initial guess vector
  	double F[n];  											 // F vector (root equation: F[X] = 0)
	// J = Jacobian of F
	double **J = (double **) malloc(n * sizeof(double *));
	for(int i = 0; i < n; i++) J[i] = (double *) malloc(n * sizeof(double));
 	int pvector[n];									  		 // permutation vector





	// prefactors and factors of F/J elements: to be evaluated at ith iteration of X
 	double commonfactori;

  	double lambdai2;
  	double lambdai3;
  	double axi2;
  	double azi2;
  	double lambdaiaxi3;
  	double lambdaiazi3;

  	double factorEai;
  	double factorPTai;
  	double factorPLai;

  	double factorI2001;
  	double factorI2011;
  	double factorI2201;

  	double factorI401m1;
  	double factorI420m1;
  	double factorI402m1; 
  	double factorI421m1;
  	double factorI440m1;

  	double factornBai;


	// anisotropic functions evaluated at ith iteration of X
	double Eai;
	double PTai;
	double PLai;

	double I2001;
	double I2011;
	double I2201;

	double I401m1;
	double I420m1;
	double I402m1;
	double I421m1;
	double I440m1;

	double nBai;

	double bn_I2000;
	double bn_I2010;
	double bn_I2200;

	double bn_I1001;
	double bn2_I1000; 

	



	int i = 0;				  // starting ith iteration
	int Nmax = 200;			  // max number of iterations
	double dXnorm2;           // L2-norm of dX iteration
	double Fnorm2;		      // L2-norm of F
	double tolX = 1.0e-8;    // tolerance for X
	double tolF = 1.0e-12;    // what's the scale I should use?

	// Find anisotropic variables using 3D Newton Method
	do
	{
		// Evaluate factors and prefactors of F/J
		lambdai2 = lambdai * lambdai;
	    lambdai3 = lambdai2 * lambdai;
	    axi2 = axi * axi;
	    azi2 = azi * azi;
	    lambdaiaxi3 = lambdai * axi * axi2;
	    lambdaiazi3 = lambdai * azi * azi2;

	    commonfactori = axi2 * azi * lambdai2 / (4.0*M_PI*M_PI);

	    factorEai = commonfactori * lambdai2;
	    factorPTai = commonfactori * axi2 * lambdai2 / 2.0;
	    factorPLai = commonfactori * azi2 * lambdai2;

	    factorI2001 = commonfactori * lambdai3;
	    factorI2011 = commonfactori * axi2 * lambdai3 / 2.0;
	    factorI2201 = commonfactori * azi2 * lambdai3;

	    factorI401m1 = factorI2011;
	    factorI420m1 = factorI2201;

	    factorI402m1 = commonfactori * axi2 * axi2 * lambdai3 / 8.0;
	    factorI421m1 = commonfactori * axi2 * azi2 * lambdai3 / 2.0;
	    factorI440m1 = commonfactori * azi2 * azi2 * lambdai3;



	    factornBai = axi2 * azi * lambdai3 / (2.0*M_PI*M_PI);



	    // reset F elements to 0
	    Eai = 0.0;
	    PTai = 0.0;
	    PLai = 0.0;

	    nBai = 0.0;

	    // Evaluate anisotropic functions for F  (1D Gauss Laguerre integrals)
	    for(int k = 0; k < number_resonances; k++)
	    {
	    	dof = (double)degeneracy[k];
	    	
	    	Eai += dof * factorEai * Gauss_Aniso_1D(Ea_integrand, pbar_rootT, pbar_weightT, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);
	    	PTai += dof * factorPTai * Gauss_Aniso_1D(PTa_integrand, pbar_rootT, pbar_weightT, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);
	    	PLai += dof * factorPLai * Gauss_Aniso_1D(PLa_integrand, pbar_rootT, pbar_weightT, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);
	    	if(baryon[k] != 0)
	    	{
	    		nBai += dof * factornBai * Gauss_Aniso_1D(nBa_integrand, pbar_rootN, pbar_weightN, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);
	    	}
	    }


	    ////////////////////////////
	    //                        //
	    //   F  =  Eai - Ea       //
	    //         PTai - PTa     //
	    //         PLai - PLa     //
	    //         nBai - nBa     //
	    //                        //
	    ////////////////////////////


		F[0] = Eai - Ea;
    	F[1] = PTai - PTa;
    	F[2] = PLai - PLa;
    	F[3] = nBai - nBa;

		// Calculate L2-norm of F
    	Fnorm2 = sqrt(fabs(F[0]*F[0] + F[1]*F[1] + F[2]*F[2] + F[3]*F[3]));


    	// reset J elements to 0

    	I2001 = 0.0;
    	I2011 = 0.0;
    	I2201 = 0.0;

    	I401m1 = 0.0;
	    I420m1 = 0.0;
	    I402m1 = 0.0;
	    I421m1 = 0.0;
	    I440m1 = 0.0;

	    bn_I2000 = 0.0;
	    bn_I2010 = 0.0;
	    bn_I2200 = 0.0;

	    bn_I1001 = 0.0; 
	    bn2_I1000 = 0.0;

    	// sum over resonances 

	    for(int k = 0; k < number_resonances; k++)
	    {
	    	dof = (double)degeneracy[k];

	    	I2001 += dof * factorI2001 * Gauss_Aniso_1D(I2001_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);
		    I2011 += dof * factorI2011 * Gauss_Aniso_1D(I2011_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);
		    I2201 += dof * factorI2201 * Gauss_Aniso_1D(I2201_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);
		    //I401m1 += dof * factorI401m1 * Gauss_Aniso_1D(I401m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);
		    //I420m1 += dof * factorI420m1 * Gauss_Aniso_1D(I420m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);
		    I402m1 += dof * factorI402m1 * Gauss_Aniso_1D(I402m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);
		    I421m1 += dof * factorI421m1 * Gauss_Aniso_1D(I421m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);
		    I440m1 += dof * factorI440m1 * Gauss_Aniso_1D(I440m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);

		    if(baryon[k] != 0)
		    {
		    	bn_I2000 += dof * factorEai * Gauss_Aniso_1D(bn_I2000_integrand, pbar_rootT, pbar_weightT, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]); 
		    	bn_I2010 += dof * factorPTai * Gauss_Aniso_1D(bn_I2010_integrand, pbar_rootT, pbar_weightT, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]); 
		    	bn_I2200 += dof * factorPLai * Gauss_Aniso_1D(bn_I2200_integrand, pbar_rootT, pbar_weightT, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);
		    	 
		    	bn_I1001 += dof * factornBai * Gauss_Aniso_1D(bn_I1001_integrand, pbar_rootT, pbar_weightT, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);
		    	bn2_I1000 += dof * factornBai * Gauss_Aniso_1D(bn2_I1000_integrand, pbar_rootN, pbar_weightN, pbar_pts, axi, azi, mass[k]/lambdai, aBti, baryon[k], sign[k]);
		    }
	    }
    	// Evaluate anisotropic functions for J (more 1D Gauss-Laguerre integrals)
	    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	    // need to add the 4th row and columns associated with nB and aB

	    ////////////////////////////////////////////////////////////////////////////////////////////
	    //                                                                                        //
	    //    J  =  I2001/lambda2       2*I401m1/lambda/ax3    I420m1/lambda/az3     bn_I2000     //
	    //                                                                                        //
	    //          I2011/lambda2       4*I402m1/lambda/ax3    I421m1/lambda/az3     bn_I2010     //
	    //                                                                                        //
	    //          I2201/lambda2       2*I421m1/lambda/ax3    I440m1/lambda/az3     bn_I2200     //
	    //         																				  //
	    //          bn_I1001/lambda2          2nBa/ax                nBa/az          bn2_I1000    //
	    //                                                                                        //
	    ////////////////////////////////////////////////////////////////////////////////////////////

	    // row 1
	    J[0][0] = I2001/lambdai2;
	    //J[0][1] = 2.0 * I401m1 / lambdaiaxi3;  // use identity J401m1 = lambda*ax2*(Ea+PTa)
	    J[0][1] = 2.0*(Eai+PTai)/axi;
	    //J[0][2] = I420m1 / lambdaiazi3;        // used identity J420m1 = lambda*az2*(Ea+PLa)
	    J[0][2] = (Eai+PLai)/azi;
	    J[0][3] = bn_I2000;   

	    // row 2
	    J[1][0] = I2011/lambdai2;
	    J[1][1] = 4.0*I402m1/lambdaiaxi3;
	    J[1][2] = I421m1/lambdaiazi3;       // identity J421m1 = ax2*az2*lambda*(PTa-PLa)/(ax2-az2) isn't numerically stable right now
	    //J[1][2] = (PTai - PLai)*axi2/(axi2-azi2)/azi;
	    J[1][3] = bn_I2010;   


	    // row 3
	    J[2][0] = I2201/lambdai2;
	    J[2][1] = 2.0 * I421m1 / lambdaiaxi3;
	    J[2][2] = I440m1 / lambdaiazi3;
	    J[2][3] = bn_I2200;   


	    // row 4
	    //J[3][0] = 3.0*nBai/lambdai;
	    J[3][0] = bn_I1001/lambdai2; 
	    J[3][1] = 2.0*nBai/axi;
	    J[3][2] = nBai/azi;
	    J[3][3] = bn2_I1000;  



	    // Solve matrix equation: J * dX = - F
	    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	    for(int k = 0; k < n; k++) F[k] = - F[k];  // change sign of F first
	    // LU solver routine
	    LUP_decomposition(J, n, pvector);          // LUP of J now stored in J (pvector also updated)
	    LUP_solve(J, n, pvector, F);               // dX is now stored in F
	    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	    // update X_i
	    for(int k = 0; k < n; k++) X[k] += F[k];


	    // calculate L2-norm of dX iteration
	    dXnorm2 = sqrt(fabs(F[0]*F[0] + F[1]*F[1] + F[2]*F[2] + F[3]*F[3]));


		// update individual variables
		lambdai = X[0];
   	 	axi = X[1];
    	azi = X[2];
    	aBti = X[3];

    	//printf("%f\n",lambdai);

		i++;

	} while((dXnorm2 > tolX) && (i < Nmax));

	cout << "Iterations: " << i << endl;

	// final answer
	*lambda = lambdai;
	*ax = axi;
	*az = azi;
	*aBt = aBti;

	// free allocated memory in Jacobian matrix
	for (int i = 0; i < n; i++) free(J[i]);
    free(J);
}










