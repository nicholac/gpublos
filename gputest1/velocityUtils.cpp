/*
 * VelocityUtils.cpp
 *
 *  Created on: 18 Mar 2015
 *      Author: dusted-ipro
 */

#include <velocityUtils.h>
#include <velocity.h>
#include <cmath>



velocity velocityUtils::calcInitVel(const double& muzzleVel, const double& az, const double& el) {
	//Az, el and combined velocity to velocity vectors
	///param az x-y azimuth in radians, measured from positive x axis anti-clockwise through to 360 (0-2pi)
	//param el elevation in radians (0=straight up) 0 to pi (but realistically we only use to 1/2 pi
	//TODO: Make this a referenced shared poitner to avoid moving memory about / copying
	velocity outVel;
	outVel.x = muzzleVel * sin(el) * cos(az);
	outVel.y = muzzleVel * sin(el) * sin(az);
	outVel.z = muzzleVel * cos(el);
	return outVel;
}



/*
VelocityUtils::~VelocityUtils() {
	// TODO Auto-generated destructor stub
}
*/
