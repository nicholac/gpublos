/*
 * motionUtils.cpp
 *
 *  Created on: 18 Mar 2015
 *      Author: dusted-ipro
 */

#include <motionUtils.h>
#include <worldParams.h>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

/**
 * Calculates new projectile position
 * \param currPos The new position of the projectile (output)
 * \param prevPos The previous position of the projectile
 * \param vel Previous velocity (at previous timestep)
 * \param deltaV Change in velocity already calculated as a result of drag
 * \param timeStep Time step thats being used in the modelling
 * \param grav Gravity acceleration constant
 *
 * @author Chris Nicholas
 */
void motionUtils::motion(std::shared_ptr<position>& currPos, const std::shared_ptr<position>& prevPos,
							const velocity& vel, const acceleration& acc, const double& timeStep) {

	//Calc resultant movement and new position
	currPos->x = prevPos->x + ((vel.x * timeStep) + (0.5*acc.x*(pow(timeStep, 2.0))));
	currPos->y = prevPos->y + ((vel.y * timeStep) + (0.5*acc.y*(pow(timeStep, 2.0))));
	currPos->z = prevPos->z + ((vel.z * timeStep) + (0.5*acc.z*(pow(timeStep, 2.0))));
}

void motionUtils::drag(double& dragForce, const double& vel, const double& mortSigma) {
	/**
	 * Calculates the drag force for a particular velocity vector
	 */
	dragForce = mortSigma*pow(vel, 2.0);
	if (vel > 0.0){
		dragForce*=-1.0;
	}
}

double motionUtils::calibreToarea(const double& calibre) {
	/**Calculates frontal area based on calibre (estimation)
	 * Calibre in metres
	 */
	double out = pi*(pow((0.5*calibre), 2.0));
	return out;
}











