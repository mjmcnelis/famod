#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string>
#include <string.h>
#include <iostream>
#include <iomanip>
using namespace std;
#include <sstream>
#include <fstream>
#include "anisotropicvariables.hpp"
#include "gauss_integration.hpp"
#include "anisotropic_integrands.hpp"

//#include <gsl/gsl_sf.h>

#define GEV_TO_INVERSE_FM 5.067731


// temporary
const int a = 21;
const int gla_pts = 64;
double root_gla[a][gla_pts];
double weight_gla[a][gla_pts];

int load_gauss_laguerre_data()
{
  FILE *fp;

  stringstream laguerre_roots_weights;
  laguerre_roots_weights << "gla_roots_weights_" << gla_pts << "_points.txt";

  if((fp = fopen(laguerre_roots_weights.str().c_str(), "r")) == NULL)
  {
     return 1;
  }
  for(int i = 0; i < a; i++)
  {
   for(int j = 0; j < gla_pts; j++)
   {
      if(fscanf(fp, "%i %lf %lf", &i, &root_gla[i][j], &weight_gla[i][j])!= 3)
      	{
        	printf("error loading roots/weights of Gauss-Laguerre Quadradture at %d %d\n", i, j);
    		return 1;
    	}
   }
  }
  fclose(fp);
  return 0;
}

void free_2D(double ** M, int n)
{
	for (int i = 0; i < n; i++) free(M[i]);
    free(M);
}



int main()
{
	// set up hadron resonance gas particles' mass, degeneracy, baryon, sign
	FILE *HRG;
	stringstream resonances;
	resonances << "pdg.dat";
	HRG = fopen(resonances.str().c_str(),"r");

	// pdg.dat contains (anti)mesons and baryons, but not antibaryons
	// so had to add antibaryons manually
	int N_mesons, N_baryons;

	fscanf(HRG, "%d", &N_mesons);	// read 1st line: number of mesons
	fscanf(HRG, "%d", &N_baryons);	// read 2nd line: number of baryons

	int N_resonances = N_mesons + 2*N_baryons;	// total number of resonances

	int particle_id;
	char name[20];
	double mass[N_resonances]; // [GeV] units in file
	double width;
	int degeneracy[N_resonances];
	int baryon[N_resonances], strange, charm, bottom, isospin;
	double charge;
	int decays;

	int m = 0; // antibaryon data marker

	// load data of mesons+baryons
	for(int k = 0; k < N_mesons+N_baryons; k++)
	{
		fscanf(HRG, "%d %s %lf %lf %d %d %d %d %d %d %lf %d", &particle_id, name, &mass[k], &width, &degeneracy[k], &baryon[k], &strange, &charm, &bottom, &isospin, &charge, &decays);

		if(baryon[k] == 1)	// manually add antibaryon data at end of array
		{
			mass[m+N_mesons+N_baryons] = mass[k];
			degeneracy[m+N_mesons+N_baryons] = degeneracy[k];
			baryon[m+N_mesons+N_baryons] = -1;
			m++;
		}
	}


	int sign[N_resonances];				   // sign array for bose/fermi distributions
	for(int k = 0; k < N_resonances; k++)  // degeneracy = 2*spin + 1
	{
		if(degeneracy[k] % 2 == 0)
			sign[k] = 1;  				   // fermion
		else if(degeneracy[k] % 2 == 1)
			sign[k] = -1; 				   // boson

		mass[k] *= GEV_TO_INVERSE_FM;      // convert resonance masses to fm^-1
	}

	fclose(HRG);



	//  Set up the pbar roots/weights for anisotropic integrands

	const int aN = 1; // gla (a = 1)
	const int aT = 2; // gla (a = 2)
	const int aJ = 3; // gla (a = 3)

	double * pbar_rootN = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_weightN = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_rootT = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_weightT = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_rootJ = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_weightJ = (double *)malloc(gla_pts * sizeof(double));


	printf("Start loading gauss data...");
	int num_error;
	if((num_error = load_gauss_laguerre_data()) != 0)
	{
		fprintf(stderr, "Error loading gauss data (%d)!\n", num_error);
		return 1;
	}
	printf("done\n\n");


	// Set momentum bar roots-weights
	for(int i = 0; i < gla_pts; i++)
	{
		pbar_rootN[i] = root_gla[aN][i];
		pbar_weightN[i] = weight_gla[aN][i];
		pbar_rootT[i] = root_gla[aT][i];
		pbar_weightT[i] = weight_gla[aT][i];
		pbar_rootJ[i] = root_gla[aJ][i];
		pbar_weightJ[i] = weight_gla[aJ][i];
	}


	// Chebyshev root-weight generator (temporary)

	// phi = M_PI * (1 + xphi)        (variable substitution)

	const int Ncheby = 32;

	double * cheby_root = (double *)malloc((Ncheby+1) * sizeof(double));
	double * cheby_weight = (double *)malloc((Ncheby+1) * sizeof(double));

	double cheby_interval = M_PI / (double)Ncheby; // spacings of the cosine argument

	cheby_weight[0] = 1.0 / ((double)Ncheby * (double)Ncheby - 1.0);
	cheby_weight[Ncheby] = cheby_weight[0];

	for(int i = 0; i < Ncheby+1; i++) cheby_root[i] = cos((double)i * cheby_interval);

	double subw = 0.0;

	if(Ncheby % 2 == 0)
	{
		for(int i = 1; i < Ncheby; i++)
		{
			for(int j = 1; j < Ncheby/2; j++)
			{
				subw += 2.0 / (1.0 - 4.0 * (double)j * (double)j) * cos(2.0 * (double)i * (double)j * cheby_interval);
			}
			cheby_weight[i] = 2.0 / (double)Ncheby * (1.0 + subw + cos((double)i * M_PI) / (1.0 - (double)Ncheby * (double)Ncheby));
			subw = 0.0;
		}
	}
	else printf("Ncheby is not even!\n");



	// set angular roots-weights
	const int angle_pts = Ncheby+1; // number of angular evaluation points
	double * xphi_root = (double *)malloc(angle_pts * sizeof(double));
	double * xphi_weight = (double *)malloc(angle_pts * sizeof(double));
	double * costheta_root = (double *)malloc(angle_pts * sizeof(double));
	double * costheta_weight = (double *)malloc(angle_pts * sizeof(double));

	for(int i = 0; i < angle_pts; i++)
	{
		xphi_root[i] = cheby_root[i];
		xphi_weight[i] = cheby_weight[i];
		costheta_root[i] = cheby_root[i];
		costheta_weight[i] = cheby_weight[i];
	}



	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//                                                       ::
	//                  TESTING GENERAL CASE                 ::
	//                                                       ::
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	const double T = 0.135 * GEV_TO_INVERSE_FM;       // temperature in fm^-1
	const double aB = 4.0;							  // chemical potential over temperature
	double ax = 0.9;
	double az = 0.9;    // it's because of the rounding errors for ax ~ az
	double lambda = 0.17 * GEV_TO_INVERSE_FM;
	double aBt = 2.8;
	// double ax = 1.11875604453775;                           // alpha_perp
	// double az = 0.752040817012744;                          // alpha_L
	// double lambda = 0.156200599527327 * GEV_TO_INVERSE_FM;  // lambda in fm^-1


    // Calculate equilibrium energy density, pressure and net baryon density of hadron resonance gas

	double factEeq = pow(T,4) / (2.0*M_PI*M_PI);
	double factPeq = pow(T,4) / (6.0*M_PI*M_PI);
	double factnBeq = pow(T,3) / (2.0*M_PI*M_PI);

	double Eeq = 0.0;
	double Peq = 0.0;
	double nBeq = 0.0;

	double dof;

	for(int k = 0; k < N_resonances; k++)
	{
		dof = (double)degeneracy[k];
		Eeq += dof * factEeq * Gauss_Thermo_1D(Eeq_integrand, pbar_rootT, pbar_weightT, gla_pts, mass[k]/T, aB, baryon[k], sign[k]);
		Peq += dof * factPeq * Gauss_Thermo_1D(Peq_integrand, pbar_rootT, pbar_weightT, gla_pts, mass[k]/T, aB, baryon[k], sign[k]);
		if(baryon[k] != 0)
		{
			nBeq += dof * factnBeq * Gauss_Thermo_1D(nBeq_integrand, pbar_rootN, pbar_weightN, gla_pts, mass[k]/T, aB, baryon[k], sign[k]);
		}
	}


	// printf("Calculate proton energy density contribution\n");
	// double mass_p = 0.938 * GEV_TO_INVERSE_FM;
	// double degeneracy_p = 2.0;

	// double Ep_numerical = degeneracy_p * factEeq * Gauss_Thermo_1D(Eeq_integrand, pbar_rootT, pbar_weightT, gla_pts, mass_p/T, aB, 1, 1);

	// double Ep_series = 0.0;
	// double prefactor;
	// int imax = 3;
	// printf("Proton energy series:\n");
	// for(int i = 1; i < imax; i++)
	// {
	// 	prefactor = degeneracy_p * pow(-1.0,i+1) * exp((double)i*aB) * factEeq / pow((double)i,4);

	// 	Ep_series += prefactor * Gauss_Thermo_1D(Eeq_integrand, pbar_rootT, pbar_weightT, gla_pts, (double)i*mass_p/T, aB, 0, 0);
	// 	printf("%f\n", Ep_series);
	// }
	// cout << "Proton energy series =    " << setprecision(15) <<  Ep_series << endl;
	// cout << "Proton energy numerical = " << setprecision(15) << Ep_numerical << "\n" << endl;


	printf("Hadron resonance gas:\n");
	cout << "e = " << setprecision(15) << Eeq / GEV_TO_INVERSE_FM << " GeV/fm^3" << endl;
	cout << "p = " << setprecision(15) << Peq / GEV_TO_INVERSE_FM << " GeV/fm^3" << endl;
	cout << "nB = " << setprecision(15) << nBeq << " fm^-3" << endl;


	// choose ahydro quantities
	double e = Eeq;
	double pt = 0.8 * Peq;
	double pl = 0.8 * Peq;
	double nB = nBeq;

	printf("\nAnisotropic hydro input:\n");
	cout << "e_a = " << setprecision(15) << Eeq << endl;
	cout << "pt_a = " << setprecision(15) << pt << endl;
	cout << "pl_a = " << setprecision(15) << pl << endl;
	cout << "nB_a = " << setprecision(15) << nB << endl;
	printf("\n");

	// bug if isotropic pressures: pt = pl
	find_anisotropic_variables(e, pl, pt, nB, mass, degeneracy, baryon, sign, N_resonances, pbar_rootN, pbar_weightN, pbar_rootT, pbar_weightT, pbar_rootJ, pbar_weightJ, gla_pts, &lambda, &ax, &az, &aBt);

	printf("\nThermodynamic variables:\n");
	cout << "T = " << setprecision(5) << T / GEV_TO_INVERSE_FM << " GeV" << endl;
	cout << "aB = " << setprecision(5) << aB << endl;
	printf("\nAnisotropic variables:\n");
	cout << "lambda = " << lambda / GEV_TO_INVERSE_FM << " GeV" << endl;
	cout << "ax = " << setprecision(8) << ax << endl;
	cout << "az = " << setprecision(8) << az << endl;
	cout << "aBt = " << setprecision(5) << aBt << endl;

	// mbar = mass / lambda
	double mbar[N_resonances];
	for(int k = 0; k < N_resonances; k++) mbar[k] = mass[k] / lambda;


	// set up prefactor_1D * 1D gauss integral
	double factEa = pow(ax,2) * pow(az,1) * pow(lambda,4) / (4.0*M_PI*M_PI);
	double factPTa = pow(ax,4) * pow(az,1) * pow(lambda,4) / (8.0*M_PI*M_PI);
	double factPLa = pow(ax,2) * pow(az,3) * pow(lambda,4) / (4.0*M_PI*M_PI);
	double factnBa = pow(ax,2) * pow(az,1) * pow(lambda,3) / (2.0*M_PI*M_PI);
	double factI402m1 = pow(ax,6) * pow(az,1) * pow(lambda,5) / (32.0*M_PI*M_PI);
	double factI421m1 = pow(ax,4) * pow(az,3) * pow(lambda,5) / (8.0*M_PI*M_PI);


	double Ea = 0.0;
	double PTa = 0.0;
	double PLa = 0.0;
	double nBa = 0.0;
	double I402m1 = 0.0;
	double I421m1 = 0.0;

	// sum over all resonances
	for(int k = 0; k < N_resonances; k++)
	{
		dof = (double)degeneracy[k];
		Ea += dof * factEa * Gauss_Aniso_1D(Ea_integrand, pbar_rootT, pbar_weightT, gla_pts, ax, az, mbar[k], aBt, baryon[k], sign[k]);
		PTa += dof * factPTa * Gauss_Aniso_1D(PTa_integrand, pbar_rootT, pbar_weightT, gla_pts, ax, az, mbar[k], aBt, baryon[k], sign[k]);
		PLa += dof * factPLa * Gauss_Aniso_1D(PLa_integrand, pbar_rootT, pbar_weightT, gla_pts, ax, az, mbar[k], aBt, baryon[k], sign[k]);
		I402m1 += dof * factI402m1 * Gauss_Aniso_1D(I402m1_integrand, pbar_rootJ, pbar_weightJ, gla_pts, ax, az, mbar[k], aBt, baryon[k], sign[k]);
		I421m1 += dof * factI421m1 * Gauss_Aniso_1D(I421m1_integrand, pbar_rootJ, pbar_weightJ, gla_pts, ax, az, mbar[k], aBt, baryon[k], sign[k]);
		if(baryon[k] != 0)
		{
			nBa += dof * factnBa * Gauss_Aniso_1D(nBa_integrand, pbar_rootN, pbar_weightN, gla_pts, ax, az, mbar[k], aBt, baryon[k], sign[k]);
		}
	}

	printf("\nAnisotropic hydro output:\n");
	cout << "Ea = " << setprecision(15) << Ea << endl;
	cout << "PTa = " << setprecision(15) << PTa << endl;
	cout << "PLa = " << setprecision(15) << PLa << endl;
	cout << "nBa = " << setprecision(15) << nBa << endl;

	// order of magnitude relations
	//double pi_order = I402m1 / (ax * ax * lambda * PTa);
	//double W_order = (ax+az) / sqrt(2.0*ax*az) * I421m1 / (ax * az * lambda * sqrt(PLa*PTa));
	// piperp << pi_order * PTa
	// Wperp << W_order * sqrt(PLa * PTa)
	//cout << "pi_order: " << setprecision(15) << pi_order << endl;
	//cout << "W_order:  " << W_order << endl;


	// viscous inputs
	double pixx = 0.5 * PTa;
	double pixy = 0.0 * PTa;
	double Wxz = 0.0 * sqrt(PLa*PTa);
	double Wyz = 0.0 * sqrt(PLa*PTa);
	double Ttx = 0.0;
	double Tty = 0.0;
	double Ttz = 0.0;


	// momentum rescaling matrix A for viscous anisotropic hydrodynamics

	// sample from anisotropic p_prime distribution
	// p = A * p_prime;

	// in the integral calculation combine the vahydro rescaling matrix A with the ahydro deformation matrix B
	// p = (A * B) p_prime is the variable transformation (the exponential is then isotropic in p_prime space)



	double factorT = ax * ax * az / (ax + az);
	double factorL = az * ax * az / (ax + az);


	const int n = 3;
  	double **A = (double **) malloc(n * sizeof(double *));
  	for(int i = 0; i < n; i++) A[i] = (double *) malloc(n * sizeof(double));



    A[0][0] = 1.0 + 0.5*pixx*ax*ax*lambda/I402m1;	A[0][1] = 0.5*pixy*ax*ax*lambda/I402m1;			A[0][2] = Wxz*factorT*lambda/I421m1;

    A[1][0] = 0.5*pixy*ax*ax*lambda/I402m1;			A[1][1] = 1.0 - 0.5*pixx*ax*ax*lambda/I402m1;	A[1][2] = Wyz*factorT*lambda/I421m1;

  	A[2][0] = Wxz*factorL*lambda/I421m1;			A[2][1] = Wyz*factorL*lambda/I421m1;			A[2][2] = 1.0;



	// renormalize particle density: f'' = particle_renormalization * f'
	// detA calculated by hand
	double detA = A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]) + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
	double particle_renormalization = 1.0 / detA;
	//cout << "Anisotropic particle density renormalization:  " << setprecision(5) << particle_renormalization << "\n" << endl;

	// detC = detA * detB: modded anisotropic particle density scale


	// particle renormalization scaling
	// double particle_renormalization = (ax * ax * az) / (ax_mod * ay_mod * az_mod);
	// cout << "Particle density ratio: " << setprecision(15) << particle_renormalization << endl;

	// detB: anisotropic particle density scale
	double detB = ax * ax * az;
	double factor = detB * pow(lambda,4) / (8.0*M_PI*M_PI);

	// set up prefactor_mod * mod_3D integral (already renormalized by 1/detA)
	double factmodEa = factor;
	double factmodPTa = 0.5 * factor;
	double factmodPLa =  factor;
	double factmodpixx = 0.5 * factor;
	double factmodpixy = factor;
	double factmodWxz =  factor;
	double factmodWyz =  factor;
	double factmodTtx =  factor;
	double factmodTty =  factor;
	double factmodTtz =  factor;

	double modEa = 0.0;
	double modPTa = 0.0;
	double modPLa = 0.0;
	double modnBa = 0.0;
	double modpixx = 0.0;
	double modpixy = 0.0;
	double modWxz = 0.0;
	double modWyz = 0.0;
	double modTtx = 0.0;
	double modTty = 0.0;
	double modTtz = 0.0;


	// Compute modification outputs (sum over resonances)
	for(int k = 0; k < N_resonances; k++)
	{
		dof = (double)degeneracy[k];

		if(baryon[k] != 0)
		{
			modnBa += dof * factnBa * Gauss_Aniso_1D(nBa_integrand, pbar_rootN, pbar_weightN, gla_pts, ax, az, mbar[k], aBt, baryon[k], sign[k]);
		}

		modEa += dof * factmodEa * Gauss_Mod_Aniso_3D(modEa_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], aBt, baryon[k], sign[k]);

		modPTa += dof * factmodPTa * Gauss_Mod_Aniso_3D(modPTa_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], aBt, baryon[k], sign[k]);

		modPLa += dof * factmodPLa * Gauss_Mod_Aniso_3D(modPLa_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], aBt, baryon[k], sign[k]);

		modpixx += dof * factmodpixx * Gauss_Mod_Aniso_3D(modpixx_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], aBt, baryon[k], sign[k]);

		modpixy += dof * factmodpixy * Gauss_Mod_Aniso_3D(modpixy_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], aBt, baryon[k], sign[k]);

		modWxz += dof * factmodWxz * Gauss_Mod_Aniso_3D(modWxz_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], aBt, baryon[k], sign[k]);

		modWyz += dof * factmodWyz * Gauss_Mod_Aniso_3D(modWyz_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], aBt, baryon[k], sign[k]);

		modTtx += dof * factmodTtx * Gauss_Mod_Aniso_3D(modTtx_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], aBt, baryon[k], sign[k]);

		modTty += dof * factmodTty * Gauss_Mod_Aniso_3D(modTty_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], aBt, baryon[k], sign[k]);

		modTtz += dof * factmodTtz * Gauss_Mod_Aniso_3D(modTtz_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], aBt, baryon[k], sign[k]);
	}


	// Txx and Tyy inputs/outputs
	double Txx = PTa + pixx;
	double Tyy = PTa - pixx;
	double modTxx = modPTa + modpixx;
	double modTyy = modPTa - modpixx;


	// print input/output results for comparison

	printf("\n");

	cout << setprecision(4) << "nB:        " << nBa << "         " << setprecision(3) << (modnBa / nBa - 1.0) * 100 << " % error" << "\n" << "modnBa:     " << setprecision(4) << modnBa << endl;

	printf("\n");

	cout << setprecision(4) << "E:        " << Ea << "         " << setprecision(3) << (modEa / Ea - 1.0) * 100 << " % error" << "\n" << "modE:     " << setprecision(4) << modEa << endl;

	printf("\n");

	cout << setprecision(4) << "PL:       " << PLa << "        " << setprecision(3) << (modPLa / PLa - 1.0) * 100 << " % error" << "\n" << "modPL:    " << setprecision(4) << modPLa << endl;

	printf("\n");

	cout << setprecision(4) << "Txx:      " << Txx << "        " << setprecision(3) << (modTxx / Txx - 1.0) * 100 << " % error" << "\n" << "modTxx:   " << setprecision(4) << modTxx << endl;

	printf("\n");

	cout << setprecision(4) << "Tyy:      " << Tyy << "        " << setprecision(3) << (modTyy / Tyy - 1.0) * 100 << " % error" << "\n" << "modTyy:   " << setprecision(4) << modTyy << endl;

	printf("\n");

	cout << setprecision(4) << "Txy:      " << pixy << "       "  << setprecision(3) << (modpixy / pixy - 1.0) * 100 << " % error" << "\n" << "modTxy:   " << setprecision(4) << modpixy << endl;

	printf("\n");

	cout << setprecision(4) << "Txz:      " << Wxz << "       "  << setprecision(3) << (modWxz / Wxz - 1.0) * 100 << " % error" << "\n" << "modTxz:   " << setprecision(4) << modWxz << endl;

	printf("\n");

	cout << setprecision(4) << "Tyz:      " << Wyz << "       "  << setprecision(3) << (modWyz / Wyz - 1.0) * 100 << " % error" << "\n" << "modTyz:   " << setprecision(4) << modWyz << endl;

	printf("\n");

	cout << setprecision(4) << "Ttx:      " << Ttx << "       "  << "\n" << "modTtx:   " << setprecision(4) << modTtx << endl;

	printf("\n");

	cout << setprecision(4) << "Tty:      " << Tty << "       "  << "\n" << "modTty:   " << setprecision(4) << modTty << endl;

	printf("\n");

	cout << setprecision(4) << "Ttz:      " << Ttz << "       "  << "\n" << "modTtz:   " << setprecision(4) << modTtz << endl;


	printf("\n");

	printf("Plots:\n");

	cout << setprecision(4) << "dE/E:         " << (modEa / Ea - 1.0) << endl;

	cout << setprecision(4) << "dPL/PL:       " << (modPLa / PLa - 1.0) << endl;

	cout << setprecision(4) << "dPT/PT:       " << (modPTa / PTa - 1.0) << endl;

	printf("\n");

	cout << setprecision(4) << "pixx/PT:      " << pixx / PTa << endl;
	cout << setprecision(4) << "modpixx/PT:   " << modpixx / PTa << endl;
	cout << setprecision(4) << "dpixx/PT:     " << (modpixx - pixx) / PTa << endl;

	printf("\n");

	cout << setprecision(4) << "Txx/PT:       " << Txx / PTa << endl;
	cout << setprecision(4) << "modTxx/PT:    " << modTxx / PTa << endl;
	cout << setprecision(4) << "dTxx/PT:      " << (modTxx - Txx) / PTa << endl;

	printf("\n");

	cout << setprecision(4) << "Tyy/PT:       " << Tyy / PTa << endl;
	cout << setprecision(4) << "modTyy/PT:    " << modTyy / PTa << endl;
	cout << setprecision(4) << "dTyy/PT:      " << (modTyy - Tyy) / PTa << endl;

	printf("\n");

	cout << setprecision(4) << "Txy/PT:       " << pixy / PTa << endl;
	cout << setprecision(4) << "modTxy/PT:    " << modpixy / PTa << endl;
	cout << setprecision(4) << "dTxy/PT:      " << (modpixy - pixy) / PTa << endl;

	printf("\n");

	cout << setprecision(4) << "Txz / sqrt(PLPT):       " << Wxz / sqrt(PLa*PTa) << endl;
	cout << setprecision(4) << "modTxz / sqrt(PLPT):    " << modWxz / sqrt(PLa*PTa)<< endl;
	cout << setprecision(4) << "dTxz / sqrt(PLPT):      " << (modWxz - Wxz) / sqrt(PLa*PTa)<< endl;

	printf("\n");

	cout << setprecision(4) << "Tyz / sqrt(PLPT):       " << Wyz / sqrt(PLa*PTa) << endl;
	cout << setprecision(4) << "modTyz / sqrt(PLPT):    " << modWyz / sqrt(PLa*PTa) << endl;
	cout << setprecision(4) << "dTyz / sqrt(PLPT):      " << (modWyz - Wyz) / sqrt(PLa*PTa) << endl;

	printf("\n\n\n");

	printf("Freeing memory...");

	free(pbar_rootN);
	free(pbar_weightN);
	free(pbar_rootT);
	free(pbar_weightT);
	free(pbar_rootJ);
	free(pbar_weightJ);
	free(cheby_root);
	free(cheby_weight);
	free(xphi_root);
	free(xphi_weight);
	free(costheta_root);
	free(costheta_weight);

	// free(particle_id);
	// free_2D(name,20);
	// free(mass); // file units = [GeV]
	// free(width);
	// free(degeneracy);
	// free(baryon);
	// free(strange);
	// free(charm);
	// free(bottom);
	// free(isospin);
	// free(charge);
	// free(decays);
	// free(sign);
	// free(mbar);


	free_2D(A,n);

	printf("done\n\n");

	return 0;
}






