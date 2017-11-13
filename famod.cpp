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
const int gla_pts = 32;
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
	//load hadron resonance gas particles' mass, degeneracy, baryon, sign
	FILE *HRG;
	stringstream resonances;
	resonances << "pdg.dat";
	HRG = fopen(resonances.str().c_str(),"r");

	// pdg.dat contains (anti)mesons and baryons, but not antibaryons
	// so add antibaryons manually
	int N_mesons, N_baryons, N_antibaryons;
	// read 1st line: number of mesons
	fscanf(HRG, "%d", &N_mesons);
	// read 2nd line: number of baryons
	fscanf(HRG, "%d", &N_baryons);

	N_antibaryons = N_baryons;

	// total number of resonances
	int N_resonances = N_mesons + N_baryons + N_antibaryons;

	int particle_id;
	char name[20];
	double mass[N_resonances]; // file units = [GeV]
	double width;
	int degeneracy[N_resonances];
	int baryon[N_resonances], strange, charm, bottom, isospin;
	double charge;
	int decays;

	int m = 0; // antibaryon marker

	// load data of mesons+baryons
	for(int k = 0; k < N_mesons+N_baryons; k++)
	{
		fscanf(HRG, "%d %s %lf %lf %d %d %d %d %d %d %lf %d", &particle_id, name, &mass[k], &width, &degeneracy[k], &baryon[k], &strange, &charm, &bottom, &isospin, &charge, &decays);

		if(baryon[k] == 1)
		{
			// manually add antibaryon data at end of array
			mass[m+N_mesons+N_baryons] = mass[k];
			degeneracy[m+N_mesons+N_baryons] = degeneracy[k];
			baryon[m+N_mesons+N_baryons] = -1;
			m++;
		}
	}

	// sign array for bose/fermi distributions
	int sign[N_resonances];
	for(int k = 0; k < N_resonances; k++)
	{
		// degeneracy = 2*spin + 1
		if(degeneracy[k] % 2 == 0)
			sign[k] = 1;  // fermions
		else if(degeneracy[k] % 2 == 1)
			sign[k] = -1; // bosons
		// convert resonance masses to fm^-1
		mass[k] *= GEV_TO_INVERSE_FM;
	}






	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	//  Set up the pbar roots/weights for:   Ea, PTa, PLa, pixx, pixy, Wxz, Wyz, Ttx, Tty, Ttyz (aT = 2)
	//                                       I402m1, I421m1 (aR = 3)
	//										 Na (aN = 1)
	const int aT = 2; // associated Laguerre polynomials a = 2 for Ea, PTa, PLa and other T^munu components
	const int aJ = 3; // associated Laguerre polynomials a = 3 for I402m1, I421m1
	double * pbar_rootT = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_weightT = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_rootJ = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_weightJ = (double *)malloc(gla_pts * sizeof(double));



	// Load gauss laguerre roots-weights
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

	const double T = 0.155 * GEV_TO_INVERSE_FM;       // temperature in fm^-1
	double ax = 1.0;
	double az = 1.0;    // it's because of the rounding errors for ax ~ az
	double lambda = T;
	// double ax = 1.11875604453775;                           // alpha_perp
	// double az = 0.752040817012744;                            // alpha_L
	// double lambda = 0.156200599527327 * GEV_TO_INVERSE_FM;  // lambda in fm^-1
	// convert resonance masses to fm^-1
	//for(int k = 0; k < number_resonances; k++) mass[k] *= GEV_TO_INVERSE_FM;




    // calculate hadron resonance gas equilibrium energy density and pressure

	double factEeq = pow(T,4) / (2.0*M_PI*M_PI);
	double factPeq = pow(T,4) / (6.0*M_PI*M_PI);

	double Eeq = 0.0, Peq = 0.0;

	double dof;

	for(int k = 0; k < N_resonances; k++)
	{
		dof = (double)degeneracy[k];
		Eeq += dof * factEeq * Gauss_Thermo_1D(Eeq_integrand, pbar_rootT, pbar_weightT, gla_pts, mass[k]/T, sign[k]);
		Peq += dof * factPeq * Gauss_Thermo_1D(Peq_integrand, pbar_rootT, pbar_weightT, gla_pts, mass[k]/T, sign[k]);
	}
	// are these reasonable values
	printf("Hadron resonance gas:\n");
	cout << "e = " << setprecision(15) << Eeq / GEV_TO_INVERSE_FM << " GeV/fm^3" << endl;
	cout << "p = " << setprecision(15) << Peq / GEV_TO_INVERSE_FM << " GeV/fm^3" << endl;


	// select the bulk pressure
	double e = Eeq;
	double pt = 1.2 * Peq;
	double pl = 0.6 * Peq;

	cout << "e = " << setprecision(15) << Eeq << endl;
	cout << "pt = " << setprecision(15) << pt << endl;
	cout << "pl = " << setprecision(15) << pl << endl;

	// there's some bug in the code at this step.

	find_anisotropic_variables(e, pl, pt, mass, degeneracy, sign, N_resonances, pbar_rootT, pbar_weightT, pbar_rootJ, pbar_weightJ, gla_pts, &lambda, &ax, &az);

	cout << T / GEV_TO_INVERSE_FM << endl;
	cout << lambda / GEV_TO_INVERSE_FM << endl;
	cout << ax << endl;
	cout << az << endl;

	// mbar = mass / lambda
	double mbar[N_resonances];
	for(int k = 0; k < N_resonances; k++) mbar[k] = mass[k] / lambda;


	// set up prefactor_1D * 1D gauss integral
	double factEa = pow(ax,2) * pow(az,1) * pow(lambda,4) / (4.0*M_PI*M_PI);
	double factPTa = pow(ax,4) * pow(az,1) * pow(lambda,4) / (8.0*M_PI*M_PI);
	double factPLa = pow(ax,2) * pow(az,3) * pow(lambda,4) / (4.0*M_PI*M_PI);
	double factI402m1 = pow(ax,6) * pow(az,1) * pow(lambda,5) / (32.0*M_PI*M_PI);
	double factI421m1 = pow(ax,4) * pow(az,3) * pow(lambda,5) / (8.0*M_PI*M_PI);


	double Ea = 0.0, PTa = 0.0, PLa = 0.0, I402m1 = 0.0, I421m1 = 0.0;

	// sum over all resonances
	for(int k = 0; k < N_resonances; k++)
	{
		dof = (double)degeneracy[k];
		Ea += dof * factEa * Gauss_Aniso_1D(Ea_integrand, pbar_rootT, pbar_weightT, gla_pts, ax, az, mbar[k], sign[k]);
		PTa += dof * factPTa * Gauss_Aniso_1D(PTa_integrand, pbar_rootT, pbar_weightT, gla_pts, ax, az, mbar[k], sign[k]);
		PLa += dof * factPLa * Gauss_Aniso_1D(PLa_integrand, pbar_rootT, pbar_weightT, gla_pts, ax, az, mbar[k], sign[k]);
		I402m1 += dof * factI402m1 * Gauss_Aniso_1D(I402m1_integrand, pbar_rootJ, pbar_weightJ, gla_pts, ax, az, mbar[k], sign[k]);
		I421m1 += dof * factI421m1 * Gauss_Aniso_1D(I421m1_integrand, pbar_rootJ, pbar_weightJ, gla_pts, ax, az, mbar[k], sign[k]);
	}

	cout << "e = " << Ea << endl;
	cout << "pt = " << PTa << endl;
	cout << "pl = " << PLa << endl;

	// order of magnitude relations
	double pi_order = I402m1 / (ax * ax * lambda * PTa);
	double W_order = (ax+az) / sqrt(2.0*ax*az) * I421m1 / (ax * az * lambda * sqrt(PLa*PTa));
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

	double modEa = 0.0, modPTa = 0.0, modPLa = 0.0, modpixx = 0.0, modpixy = 0.0,
	modWxz = 0.0, modWyz = 0.0, modTtx = 0.0, modTty = 0.0, modTtz = 0.0;


	// Compute modification outputs (sum over resonances)
	for(int k = 0; k < N_resonances; k++)
	{
		dof = (double)degeneracy[k];

		modEa += dof * factmodEa * Gauss_Mod_Aniso_3D(modEa_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], sign[k]);

		modPTa += dof * factmodPTa * Gauss_Mod_Aniso_3D(modPTa_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], sign[k]);

		modPLa += dof * factmodPLa * Gauss_Mod_Aniso_3D(modPLa_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], sign[k]);

		modpixx += dof * factmodpixx * Gauss_Mod_Aniso_3D(modpixx_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], sign[k]);

		modpixy += dof * factmodpixy * Gauss_Mod_Aniso_3D(modpixy_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], sign[k]);

		modWxz += dof * factmodWxz * Gauss_Mod_Aniso_3D(modWxz_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], sign[k]);

		modWyz += dof * factmodWyz * Gauss_Mod_Aniso_3D(modWyz_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], sign[k]);

		modTtx += dof * factmodTtx * Gauss_Mod_Aniso_3D(modTtx_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], sign[k]);

		modTty += dof * factmodTty * Gauss_Mod_Aniso_3D(modTty_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], sign[k]);

		modTtz += dof * factmodTtz * Gauss_Mod_Aniso_3D(modTtz_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, ax, az, A, n, mbar[k], sign[k]);
	}


	// Txx and Tyy inputs/outputs
	double Txx = PTa + pixx;
	double Tyy = PTa - pixx;
	double modTxx = modPTa + modpixx;
	double modTyy = modPTa - modpixx;


	// print input/output results for comparison

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






