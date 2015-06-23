/*
 * motionUtils.h
 *
 *  Created on: 18 Mar 2015
 *      Author: dusted-ipro
 */

#ifndef SRC_MOTIONUTILS_H_
#define SRC_MOTIONUTILS_H_

#include <memory>
#include <position.h>
#include <velocity.h>
#include <acceleration.h>

class motionUtils {
	public:

		/*! Alters position based on velocity and acceleration at given timeStep
		\param position currPos New position after movement (output)
		\param position prevPs Previous position (before movement)
		\param velocity vel x,y,z velocity (m/s) at previous position
		\param acceleration acc x,y,z acceleration (m/s^2) for this move
		\param double timeStep modelling step time (secs)
		\return void altered currPos
		\sa drag, recShot::launchAir
		*/
		static void motion(std::shared_ptr<position>& currPos, const std::shared_ptr<position>& prevPos,
									const velocity& vel, const acceleration& acc, const double& timeStep);

		/*! Drag function
		\param double vel velocity in one dimension (m/s)
		\param double mortSigma pre-computed partial drag equation (drag coef, area, airdens)
		\param double dragForce (output)
		\return void altered dragForce
		\sa drag, recShot::launchAir, worldParams::sigma, main::mortMuzVel
		*/
		static void drag(double& dragForce, const double& vel, const double& mortSigma);

		/**Frontal area estimation from calibre (area of circle from calibre)
		 \param calibre Calibre of weapon (metres)
	 	 \return double area Frontal area estimation based on calibre
	 	 \sa drag, motion, worldParams::sigma, main::calibre
		 */
		static double calibreToarea(const double& calibre);

};

#endif /* SRC_MOTIONUTILS_H_ */
