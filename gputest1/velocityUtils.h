/*
 * VelocityUtils.h
 *
 *  Created on: 18 Mar 2015
 *      Author: dusted-ipro
 */

#ifndef SRC_VELOCITYUTILS_H_
#define SRC_VELOCITYUTILS_H_
/*
struct velocity {
	//3d velocity vector
	double x;
	double y;
	double z;
};
*/

#include <velocity.h>

class velocityUtils {
public:
	static velocity calcInitVel(const double& muzzleVel, const double& az, const double& el);
};

#endif /* SRC_VELOCITYUTILS_H_ */
