/*
 * worldparams.h
 *
 *  Created on: 7 Apr 2015
 *      Author: dusted-ipro
 */

#ifndef SRC_WORLDPARAMS_H_
#define SRC_WORLDPARAMS_H_


//3d velocity vector
extern double gravity;
extern double airDens;
//Coef for Sphere
extern double dragCoef;
extern double pi;

/**Drag Equation pre-calculation
 * \param dragCoef Coefficient of drag for the object (e.g. 0.42 sphere, 0.04 streamlined body)
 * \param airDens Density of air (1.25 at ground level)
 * \para area Reference frontal area of body
 * \sa motionUtils
 */

class worldParams {
public:
		static double sigma(const double& area);
};


#endif /* SRC_WORLDPARAMS_H_ */
