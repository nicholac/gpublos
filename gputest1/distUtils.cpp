/*
 * distUtils.cpp
 *
 *  Created on: 18 Mar 2015
 *      Author: dusted-ipro
 */

#include <distUtils.h>
#include <cmath>
#include <position.h>
#include <azUtils.h>
#include <velocityUtils.h>
#include <velocity.h>
#include <motionUtils.h>
#include <recShot.h>
#include <iostream>


double distUtils::vecDist(const std::shared_ptr<position>& pos1,
								const std::shared_ptr<position>& pos2) {
	//Euclid distance between two 3d points
	double dX = pos2->x-pos1->x;
	double dY = pos2->y-pos1->y;
	double dZ = pos2->z-pos1->z;
	double distance = sqrt(pow(dX, 2.0)+pow(dY, 2.0)+pow(dZ, 2.0));
	return distance;
}

double distUtils::vecDist2d(const std::shared_ptr<position>& pos1,
								const std::shared_ptr<position>& pos2) {
	//2d (x,y) Euclid distance between two 3d points
	double dX = pos2->x-pos1->x;
	double dY = pos2->y-pos1->y;
	double distance = sqrt(pow(dX, 2.0)+pow(dY, 2.0));
	return distance;
}

double distUtils::dotProd(const std::shared_ptr<position>& launchPos,
						const std::shared_ptr<position>& tgtPos,
						const std::shared_ptr<position>& missPos) {
	//Dot product calc for launch to landing etc
	double vec1x;
	double vec1y;
	double vec2x;
	double vec2y;
	double out;
	//Get differences
	vec1x = launchPos->x - tgtPos->x;
	vec1y = launchPos->y - tgtPos->y;
	vec2x = missPos->x - tgtPos->x;
	vec2y = missPos->y - tgtPos->y;
	//Do dot product
	out = (vec1x*vec2x)+(vec1y*vec2y);
	return out;
}

double distUtils::mortorRange(const double& initVel, const double& timeStep, const double& el, const double& mortSigma,
								const double& mortMass) {
	/*
	*  Calculates the max range for a mortor
	*  Input: Launch position, initial velocity
	*  TODO: Change this to use the same recording shot as everything else
	*/
	velocity vel;
	auto launchPos = std::make_shared<position>( position {0,0,0});
	auto tgtPos = std::make_shared<position>( position {500,500,0});
	double dist = 0.0;
	//Run this shot
	auto trj = recShot::launchAir(launchPos, tgtPos, initVel, el, timeStep, mortSigma, mortMass);
	//cout << trj.back()->x << " // " << trj.back()->y << " // " << trj.back()->z << endl;
	//Work out the range the shot covered
	dist = distUtils::vecDist2d(launchPos, trj.back());
	return dist;
}






