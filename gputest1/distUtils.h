/*
 * distUtils.h
 *
 *  Created on: 18 Mar 2015
 *      Author: dusted-ipro
 */

#ifndef SRC_DISTUTILS_H_
#define SRC_DISTUTILS_H_

#include <vector>
#include <memory> // this is for the safe pointer wrapper
#include <position.h>
using namespace std;


class distUtils {
	public:

		//! 3D Distance Between two positions
		/*!
		\param pos1: 3D Position contained within a shared_ptr
		\param pos2: 3D Position contained within a shared_ptr
		\return distance: distance between the two points
		\sa vecDist2d
		*/
		static double vecDist(const std::shared_ptr<position>& pos1,
								const std::shared_ptr<position>& pos2);

		//! 2D Distance Between two positions
		/*!
		\param pos1: 3D Position contained within a shared_ptr
		\param pos2: 3D Position contained within a shared_ptr
		\return distance: 2d (x, y) distance between the two points
		\sa vecDist
		*/
		static double vecDist2d(const std::shared_ptr<position>& pos1,
										const std::shared_ptr<position>& pos2);

		//! Dot Product of the two vectors made from three points: launch-->tgt, miss-->tgt.
		/*!
		\param launchPos: 3D Position contained within a shared_ptr
		\param tgtPos: 3D Position contained within a shared_ptr
		\param missPos: 3D Position contained within a shared_ptr
		\return dotProduct
		\sa optUtils::optimise
		*/
		static double dotProd(const std::shared_ptr<position>& launchPos,
								const std::shared_ptr<position>& tgtPos,
								const std::shared_ptr<position>& missPos);

		//! Maximum range for a mortor (2d)
		/*!
		\param initVel: Muzzle Velocity
		\param timeStep: Accuracy of the modelling - in seconds.
		\param el: Barrel elevation in degrees from horizontal (0 to 90)
		\return distance: 2d (x, y) distance from launch to impact
		\sa predictUtils
		*/
		static double mortorRange(const double& initVel, const double& timeStep, const double& el, const double& mortSigma,
									const double& mortMass);

};


#endif /* SRC_DISTUTILS_H_ */
