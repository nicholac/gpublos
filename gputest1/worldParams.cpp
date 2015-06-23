/*
 * worldParams.cpp
 *
 *  Created on: 7 Apr 2015
 *      Author: dusted-ipro
 */


#include <worldParams.h>


double gravity = 9.82;
//TODO: Alter this based on the altitude
double airDens = 1.225;
//TODO: Allow setting of this based on weapon vars
//Using Sreamlined body coefficient
double dragCoef = 0.04; //Coef for sphere
double pi = 3.1415926; // pi

double worldParams::sigma(const double& area) {
	/** Pre calculate drag equation part
	 *
	 */
	double out = dragCoef*area*0.5*airDens;
	return out;
}






