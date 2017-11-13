
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hypergeometric.hpp"

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                      HYPERGEOMETRIC FUNCTIONS                    ::
//                                                                  ::
//      Hypergeometric functions that appear in the R-functions     ::
//      of the anisotropic integrands Ea, PTa, PLa, I402-1, and     ::
//		I421-1 respectively                                         ::
//                                                                  ::
//      tfuncE     tfuncPT     tfuncPL     tfunc402     tfunc421    ::
//                                                                  ::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


double tfuncE(double x)
{
	double result = 0.0;
	// if(x > -1.0)
	// {
	// 	double x1 = (double)x;
	// 	double x2 = x1 * x1;
	// 	double x3 = x2 * x1;
	// 	double x4 = x3 * x1;
	// 	double x5 = x4 * x1;
	// 	double x6 = x5 * x1;
	// 	double x7 = x6 * x1;
	// 	double x8 = x7 * x1;
	// 	double x9 = x8 * x1;
	// 	double x10 = x9 * x1;
	// 	double x11 = x10 * x1;
	// 	double x12 = x11 * x1;
	// 	double x13 = x12 * x1;
	// 	double x14 = x13 * x1;
	// 	double x15 = x14 * x1;
	// 	double x16 = x15 * x1;
	// 	result = (63.21975933390082 + 380.08563149296145*x1 + 810.2307212043552*x2 + 585.0129686787305*x3 - 406.6714931011936*x4 - 1003.2879972442737*x5 -
	// 	633.4140610182495*x6 - 63.549331542631045*x7 + 114.54636074466283*x8 + 62.011233002217395*x9 + 13.750682148451842*x10 +
	// 	1.4604810282159304*x11 + 0.07252235160524895*x12 + 0.001567857330443711*x13 + 0.000013272370337623199*x14 +
	// 	3.660522335594211e-8*x15 + 2.1251579740535422e-11*x16)/
	// 	(31.609879640022022 + 179.5061893290641*x + 347.38729112508554*x2 + 187.77466277293303*x3 - 247.39515928702926*x4 -
	// 	414.05598505483044*x5 - 196.62441712412095*x6 + 13.796301325027898*x7 + 47.20534483149715*x8 + 17.872426954407864*x9 +
	// 	2.9248857075582086*x10 + 0.22164902064018988*x11 + 0.007474229419816371*x12 + 0.00010364823859082527*x13 +
	// 	5.187322323753321e-7*x14 + 7.174936385850154e-10*x15 + 1.175247329467607e-13*x16);
	// }
	if(x > 0.0)
		result = 1.0 + (1.0 + x) * atan(sqrt(x))/sqrt(x);
	else if(x < 0.0 && x > -1.0)
		result = 1.0 + (1.0 + x) * atanh(sqrt(-x))/sqrt(-x);
	else if(x == 0.0)
		result = 2.0;
	else if(x <= -1.0)
		throw "tfuncE outside domain";
	return result;
}



double tfuncPT(double x)
{
	double result = 0.0;
	// if(x > -1.0)
	// {
	// 	double x1 = (double)x;
	// 	double x2 = x1 * x1;
	// 	double x3 = x2 * x1;
	// 	double x4 = x3 * x1;
	// 	double x5 = x4 * x1;
	// 	double x6 = x5 * x1;
	// 	double x7 = x6 * x1;
	// 	double x8 = x7 * x1;
	// 	double x9 = x8 * x1;
	// 	double x10 = x9 * x1;
	// 	double x11 = x10 * x1;
	// 	double x12 = x11 * x1;
	// 	double x13 = x12 * x1;
	// 	double x14 = x13 * x1;
	// 	double x15 = x14 * x1;
	// 	double x16 = x15 * x1;
	// 	result = (-5727.780276245886 - 16656.917423117204*x1 - 12506.70489796712*x2 + 3533.6365352963835*x3 + 4046.8720268359707*x4 - 3633.102683844837*x5 +
	//      396.3241468432578*x6 + 10198.392446287147*x7 + 14219.720976200637*x8 + 10469.868083059493*x9 + 4309.651396082802*x10 +
	//      908.7020129548255*x11 + 85.53600401445127*x12 + 3.0784088039656656*x13 + 0.035776352295545606*x14 +
	//      0.00010211822898423643*x15 + 3.095818533435384e-8*x16)/
	//    	(-4295.835228827275 - 14211.022782730835*x - 13959.795028128196*x2 - 97.67311371562538*x3 + 4529.779558944654*x4 -
	//      1934.3935160008843*x5 - 871.0058673904322*x6 + 8044.423723562404*x7 + 13564.43982013535*x8 + 11468.467869290651*x9 +
	//      5647.292037835587*x10 + 1541.8604492584877*x11 + 208.53309737071402*x12 + 11.992399966996727*x13 + 0.24786621578918133*x14 +
	//      0.0014919746956868216*x15 + 1.689074902336641e-6*x16);
 //   	}
	if(x > 0.0)
		// can I get rounding error propagation when |x| is very small?
		// yes, accuracy to 4/3 turns around when |x| = 1e-9 and lower
		result = (1.0 + (x - 1.0) * atan(sqrt(x))/sqrt(x)) / x;
	else if(x < 0.0 && x > -1.0)
		result = (1.0 + (x - 1.0) * atanh(sqrt(-x))/sqrt(-x)) / x;
	else if(x == 0.0)
		result = 4.0 / 3.0;
	else if(x <= -1.0)
		throw "tfuncPT outside domain";
	return result;
}



double tfuncPL(double x)
{
	double result = 0.0;
	// if(x > -1.0)
	// {
	// 	double x1 = (double)x;
	// 	double x2 = x1 * x1;
	// 	double x3 = x2 * x1;
	// 	double x4 = x3 * x1;
	// 	double x5 = x4 * x1;
	// 	double x6 = x5 * x1;
	// 	double x7 = x6 * x1;
	// 	double x8 = x7 * x1;
	// 	double x9 = x8 * x1;
	// 	double x10 = x9 * x1;
	// 	double x11 = x10 * x1;
	// 	double x12 = x11 * x1;
	// 	double x13 = x12 * x1;
	// 	double x14 = x13 * x1;
	// 	double x15 = x14 * x1;
	// 	double x16 = x15 * x1;
	// 	result = (456.0717685013072 + 6138.497799670426*x1 + 18415.553046770605*x2 + 19691.7172421894*x3 + 2898.2127968410837*x4 - 3614.6523331460667*x5 +
	//      11389.611697391874*x6 + 21474.826796208286*x7 + 14871.90337130695*x8 + 5219.601958438456*x9 + 964.6137551732895*x10 +
	//      89.7736374260045*x11 + 3.866194043736338*x12 + 0.06912592028121776*x13 + 0.00044196551466947056*x14 +
	//      7.664395594410113e-7*x15 + 1.539335288358236e-10*x16)/
	//    	(684.1076528076261 + 9344.568232901041*x + 29433.605444732875*x2 + 34655.91053674811*x3 + 9179.8694508414*x4 - 5423.728248733478*x5 +
	//      16156.69594073634*x6 + 35777.20320424328*x7 + 27919.78865897606*x8 + 11204.510344427128*x9 + 2454.1159153285016*x10 +
	//      285.7976201827534*x11 + 16.501831546099268*x12 + 0.42876486959188886*x13 + 0.0044186084191126206*x14 +
	//      0.000014810477376582228*x15 + 1.0341668485008954e-8*x16);
	// }
	if(x > 0.0)
		// rounding error propagation when |x| is very small
		// accuracy to 2/3 turns around when |x| = 1e-9 and lower
		result = (-1.0 + (x + 1.0) * atan(sqrt(x))/sqrt(x)) / x;
	else if(x < 0.0 && x > -1.0)
		result = (-1.0 + (x + 1.0) * atanh(sqrt(-x))/sqrt(-x)) / x;
	else if(x == 0.0)
		result = 2.0 / 3.0;
	else if(x <= -1.0)
		throw "tfuncPL is outside domain";
	return result;
}



double tfunc402(double x)
{
	double result = 0.0;
	// if(x > -1.0)
	// {
	// 	double x1 = (double)x;
	// 	double x2 = x1 * x1;
	// 	double x3 = x2 * x1;
	// 	double x4 = x3 * x1;
	// 	double x5 = x4 * x1;
	// 	double x6 = x5 * x1;
	// 	double x7 = x6 * x1;
	// 	double x8 = x7 * x1;
	// 	double x9 = x8 * x1;
	// 	double x10 = x9 * x1;
	// 	double x11 = x10 * x1;
	// 	double x12 = x11 * x1;
	// 	double x13 = x12 * x1;
	// 	double x14 = x13 * x1;
	// 	double x15 = x14 * x1;
	// 	double x16 = x15 * x1;
	// 	result = (-782.6241458110353 + 2877.4937140681823*x1 + 11666.945973644832*x2 + 8209.405189330819*x3 - 1491.7416191099312*x4 + 2351.8535113214166*x5 +
	//      4676.460585239748*x6 + 233.1658473312586*x7 + 3847.444607895792*x8 + 8218.087686921586*x9 + 5185.441182809942*x10 +
	//      1350.9463388535519*x11 + 144.6813790307555*x12 + 5.730779869275215*x13 + 0.07199605279826411*x14 + 0.0002191395134614306*x15 +
	//      7.00328863503514e-8*x16)/
	//    (-733.7101279175023 + 2383.202992648816*x + 12168.766064060474*x2 + 12071.780843634386*x3 + 942.438393895715*x4 + 1269.190777932744*x5 +
	//      5587.3197958329*x6 + 1747.4662503295506*x7 + 3439.1378626879673*x8 + 9328.567205945315*x9 + 7754.063179602607*x10 +
	//      2770.050106683396*x11 + 439.9397079946307*x12 + 28.332872673857835*x13 + 0.639908201951777*x14 + 0.004139070893590553*x15 +
	//      4.968575689502532e-6*x16);
	// }
	if(x > 0.0)
		// more prone to rounding error when x is small due to x^2 / x^2
		// accuracy to 16/15 turns around when |x| = 1e-6 and lower  (can be numerically unstable)
		result = (3.0 * (x - 1.0) + (3.0 + (x * (3.0*x - 2.0))) * atan(sqrt(x))/sqrt(x)) / (4.0*x*x);
	else if(x < 0.0 && x > -1.0)
		result = (3.0 * (x - 1.0) + (3.0 + (x * (3.0*x - 2.0))) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x*x);
	else if(x == 0.0)
		result = 16.0 / 15.0;
	else if(x <= -1.0)
		throw "tfunc402 outside domain";
	return result;
}



double tfunc421(double x)
{
	double result = 0.0;
	// if(x > -1.0)
	// {
	// 	double x1 = (double)x;
	// 	double x2 = x1 * x1;
	// 	double x3 = x2 * x1;
	// 	double x4 = x3 * x1;
	// 	double x5 = x4 * x1;
	// 	double x6 = x5 * x1;
	// 	double x7 = x6 * x1;
	// 	double x8 = x7 * x1;
	// 	double x9 = x8 * x1;
	// 	double x10 = x9 * x1;
	// 	double x11 = x10 * x1;
	// 	double x12 = x11 * x1;
	// 	double x13 = x12 * x1;
	// 	double x14 = x13 * x1;
	// 	double x15 = x14 * x1;
	// 	double x16 = x15 * x1;
	// 	result = (-420.65125779019604 + 464.28029537808857*x1 + 3246.0671906277103*x2 + 2098.237776364318*x3 - 1061.9549124835105*x4 + 532.7983558552136*x5 +
	//      914.6141784226661*x6 - 2418.6530038230967*x7 - 2652.892648012059*x8 - 416.8478662195916*x9 + 363.4944854153988*x10 +
	//      145.17547168536962*x11 + 16.972162791609925*x12 + 0.6693300232722256*x13 + 0.008148267235988629*x14 +
	//      0.000023934980557864453*x15 + 7.4068454059322425e-9*x16)/
	//    (-1577.442188404056 + 1290.35348266619*x + 12766.771913831271*x2 + 11195.127927858319*x3 - 2403.8898621960525*x4 + 675.83675733029*x5 +
	//      4295.662325468036*x6 - 8343.96307877952*x7 - 12618.627251409462*x8 - 3834.7232116058217*x9 + 1267.6644115040936*x10 +
	//      897.4182513912277*x11 + 159.92103119712186*x12 + 10.237581260406852*x13 + 0.22263552961633887*x14 + 0.001382377475563293*x15 +
	//      1.6010885095334335e-6*x16);
	// }
	if(x > 0.0)
		// more prone to rounding error when x is small due to x^2 / x^2
		// accuracy to 4/15 turns around when |x| = 1e-6 and lower  (can be numerically unstable)
		result = (3.0 + x + (x + 1.0) * (x - 3.0) * atan(sqrt(x))/sqrt(x)) / (4.0*x*x);
	else if(x < 0.0 && x > -1.0)
		result = (3.0 + x + (x + 1.0) * (x - 3.0) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x*x);
	else if(x == 0.0)
		result = 4.0 / 15.0;
	else if(x <= -1.0)
		throw "tfunc421 outside domain";
	return result;
}


double tfunc401(double x)
{
	double result = 0.0;
	if(x > 0.0)
		// accuracy to 4/3 turns around when |x| = 1e-10 and lower (rounding error)
		result = (1.0 + 3.0*x + (3.0*x - 1.0) * (x + 1.0) * atan(sqrt(x))/sqrt(x)) / (4.0*x);
	else if(x < 0.0 && x > -1.0)
		result = (1.0 + 3.0*x + (3.0*x - 1.0) * (x + 1.0) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x);
	else if(x == 0.0)
		result = 4.0 / 3.0;
	else if(x <= -1.0)
		throw "tfunc401 outside domain";
	return result;
}


double tfunc420(double x)
{
	double result = 0.0;
	if(x > 0.0)
		result = (x - 1.0 + (x + 1.0) * (x + 1.0) * atan(sqrt(x))/sqrt(x)) / (4.0*x);
	else if(x < 0.0 && x > -1.0)
		result = (x - 1.0 + (x + 1.0) * (x + 1.0) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x);
	else if(x == 0.0)
		result = 2.0 / 3.0;
	else if(x <= -1.0)
		throw "tfunc420 outside domain";
	return result;
}


double tfunc440(double x)
{
	double result = 0.0;
	if(x > 0.0)
		// more prone to rounding error when x is small due to x^2 / x^2
		// accuracy to 2/5 turns around when |x| = 1e-6 and lower  (can be numerically unstable)
		result = (-(3.0 + 5.0*x) + 3.0 * (x + 1.0) * (x + 1.0) * atan(sqrt(x))/sqrt(x)) / (4.0*x*x);
	else if(x < 0.0 && x > -1.0)
		result = (-(3.0 + 5.0*x) + 3.0 * (x + 1.0) * (x + 1.0) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x*x);
	else if(x == 0.0)
		result = 0.4;
	else if(x <= -1.0)
		throw "tfunc440 outside domain";
	return result;
}




