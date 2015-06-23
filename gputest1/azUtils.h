/*
 * AzUtils.h
 *
 *  Created on: 18 Mar 2015
 *      Author: dusted-ipro
 */

#ifndef SRC_AZUTILS_H_
#define SRC_AZUTILS_H_

#include <memory>
#include <position.h>
#include <vector>

class azUtils {
public:
	//! Calculate Azimuth between two points
	/*!
	\param x0, y0 double xPosition
	\param x1, y1 double yPosition
	\return double Azimuth
	\sa deg2rad, rad2deg
	*/
	static double calcAz(const double& x0, const double& y0, const double& x1, const double& y1);

	//! Degrees to Radians
	/*!
	\param double degrees
	\return double radians
	\sa rad2deg
	*/
	static double deg2Rad(const double& degrees);

	//! Radians to Degrees
	/*!
	\param double radians
	\return double degrees
	\sa deg2rad
	*/
	static float rad2Deg(const float& radians);


	//! Normalise Tgt position coords with launch position as 0.0, 0.0 (origin)
	/*! UTM coords are to nearest whole metre
	\param double inUTMx Launch UTM X position
	\param double inUTMy Launch UTM X position
	\param double inUTMx Target UTM X position - normed
	\param double inUTMy Target UTM X position - normed
	\param double outUnitX Target output X position
	\param double outUnitY Target output Y position
	\return void altered UTM X and Y positions
	\sa calcAz
	*/
	static void normPositions(const double& launchUTMx, const double& launchUTMy,
								const double& tgtUTMx, const double& tgtUTMy,
								double& outLaunchX, double& outLaunchY, double& outTgtX, double& outTgtY);

	//! Un-Normalise a trajectory of positions back to launch position as UTM Input
	/*!
	\param double inUTMx Launch UTM X position
	\param double inUTMy Launch UTM X position
	\return void Altered Vector of the trajectory
	\sa normPositions
	*/
	static void unNormPositions(const double& launchUTMx, const double& launchUTMy,
									std::vector<std::shared_ptr<position>>& traj);

};

#endif /* SRC_AZUTILS_H_ */
