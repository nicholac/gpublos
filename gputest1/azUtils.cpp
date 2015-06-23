/*
 * AzUtils.cpp
 *
 *  Created on: 18 Mar 2015
 *      Author: dusted-ipro
 */

#include <azUtils.h>
#include <position.h>
#include <worldParams.h>
#include <cmath>
#include <memory>
#include <vector>


double azUtils::calcAz(const double& x0, const double& y0, const double& x1, const double& y1) {
	//Azimuth between two points - outputs in radians, converts output to 0-2pi for spherical conversion
	//Range is 0-2pi from x axis anti-clockwise round
	//Outputs in radians
	double dX = x1-x0;
	double dY = y1-y0;
	double theta = atan2(dY, dX);
	if (theta < 0.0) {
		theta+=(2.0*pi);
	}
	return theta;
}

double azUtils::deg2Rad(const double& degrees) {
	double radian;
	radian = degrees * (pi/180.0);
	return radian;
}

float azUtils::rad2Deg(const float& radians)  {
	float degrees;
	degrees = radians * (180.0/pi);
	return degrees;
}

void azUtils::normPositions(const double& launchUTMx, const double& launchUTMy,
		const double& tgtUTMx, const double& tgtUTMy,
		double& outLaunchX, double& outLaunchY, double& outTgtX, double& outTgtY) {

	//Input coords to nearest metre with top left as coord, so for processing we move to pixel centre (+0.5)
	outTgtX = tgtUTMx - launchUTMx;
	outTgtY = tgtUTMy - launchUTMy;
	outLaunchX = 0.0;
	outLaunchY = 0.0;
}

void azUtils::unNormPositions(const double& launchUTMx, const double& launchUTMy, std::vector<std::shared_ptr<position>>& traj) {

	for (size_t i=0; i<traj.size(); i++) {
		traj[i]->x+=launchUTMx;
		traj[i]->y+=launchUTMy;
	}
}


