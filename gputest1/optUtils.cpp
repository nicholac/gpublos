/*
 * optimise.cpp
 *
 *  Created on: 18 Mar 2015
 *      Author: dusted-ipro
 */

#include <optUtils.h>
#include <recShot.h>
#include <vector>
#include <memory> // this is for the safe pointer wrapper
#include <position.h>
#include <velocity.h>
#include <velocityUtils.h>
#include <azUtils.h>
#include <motionUtils.h>
#include <distUtils.h>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <outputUtils.h>
using namespace std;

double optUtils::optimise(const std::shared_ptr<position>& launchPos, const std::shared_ptr<position>& tgtPos,
							const double& mortMuzVel, const double& tgtSize, const double& timeStep,
							const double& mortSigma, const double& mortMass, const double& seedEl) {

	/*
	 * Shot optimisation class - Finds Optimum Elevation using dot products between landing positions
	 * Assumes basic trajectory checks before running
	 * Can be seeded with close elevation (e.g. from firing table) to lower iteration count
	 * No checks in here for shots that wont reach zenith - these should have been removed pre-running
	 * Output: Optimum shot elevation (degrees)
	 * -555 - recShot trajectory fail (zenith not reached most likely)
	 */

	//Finds best shot / first angle for hitting target
	//Seed jump in degrees
	double seedJump = 0.5;
	bool contin = true;
	//This is the flag for reversing the shot optimiser
	int i = 0;
	//Decreasing from max elevation start
	double el = seedEl;
	//float el = 85.0;
	double dotP = 0.0;
	double prevDot = 0.0;
	bool zenithChk = false;
	//Hit outputs:
	//hit - std::vector with trajectory
	//Optimiser not reaching zenith = -777.0
	//Optimiser taking too long = -888.0
	//Other fail: -999.0
	auto missPos = std::make_shared<position>(position {0,0,0});
	//Out elevation
	double optEl = 0.0;

	//Optimise start
	while (contin == true) {
		//Launch - recording
		auto outTrjPtr = recShot::launchAir(launchPos, tgtPos, mortMuzVel, el, timeStep, mortSigma, mortMass);
		missPos = outTrjPtr.back();
		if (missPos->x == -99.0){
			//The Shot failed (zenith etc) - we shouldnt get to this position
			optEl = -555.0;
			contin = false;

		}
		//cout << "Dist: " << distUtils::vecDist(tgtPos, missPos) << endl;
		//Check Distance from last point to tgt - 2d
		if ( distUtils::vecDist(tgtPos, missPos) < tgtSize ) {
			//Output the optimal elevation for shot
			optEl = el;
			break;
		}
		else {
			//Missed - optimise
			//if first run - launch again with arbitrary seed
			if (i == 0) {
				el=el-seedJump;
				//Work out dot products using last recorded position
				chkDotP(tgtPos, launchPos, missPos, dotP);
				}

			else {
				//Record previous dotP and calc current
				prevDot = dotP;
				chkDotP(tgtPos, launchPos, missPos, dotP);
				//Optimise the elevation
				optElev(dotP, prevDot, el, seedJump, contin, el);
			}
		}
		//Iter Counter
		i++;

		//Abort if the elevation is not reaching tgt
		if (el < 0.0){
			optEl = -666.0;
			contin = false;
		}
		//Check if we are achieving tgt height if optimiser running too long
		if (std::abs(seedJump) < 0.00005) {
			if ( zenithChk == true ) {
				//We've already checked zenith and its good - its failing to find solution - break
				optEl = -888.0;
				contin = false;
			}
			else {
				//Check we are optimising and reaching the tgt height
				if ( maxZenith(outTrjPtr) > tgtPos->z ) {
					//Remember we've checked it but dont break
					zenithChk = true;
				}
				else {
					//Not achieving tgt height - fail
					optEl = -777.0;
					contin = false;
				}

			}
		}

	} //Opt while fin
	//Cleanup
	//Output
	return optEl;
}


double optUtils::maxZenith(const std::vector<std::shared_ptr<position>>& inTrj) {
	//Gets maximum height (at zenith) of input trajectory
	double zChk = 0.0;
	for (std::shared_ptr<position> pos : inTrj) {
		if ( pos->z > zChk ) {
			zChk = pos->z;
		}
	}
	return zChk;
}

void optUtils::chkDotP(const std::shared_ptr<position>& tgtPos, const std::shared_ptr<position>& launchPos,
					const std::shared_ptr<position>& missPos, double& dotP) {
	//Work out dot products using last recorded position
	dotP = distUtils::dotProd(launchPos, tgtPos, missPos);
}


void optUtils::optElev(const double& dotP, const double& prevDot,
				double& el, double& seedJump, bool& contin, double& optEl) {
	//Work out direction we are out
	if (dotP < 0.0) {
		if (prevDot < 0.0) {
			//Too far and prev was too far (negative dot prod between the vectors)
			//Increase elevation, no change in seed
			el = el+seedJump;
		}
		if (prevDot > 0.0) {
			//Too far and prev too close:
			//alter seedjump and increase el
			seedJump=seedJump/4.0;
			el=el+seedJump;
		}
	}
	else if (dotP > 0.0) {
		if (prevDot > 0.0) {
			//Too close and prev too close
			//decrease elevation, no change seed
			el = el-seedJump;
		}
		if (prevDot < 0.0) {
			//Too close and prev too far
			//Change seed, decrease el
			seedJump = seedJump/4.0;
			el = el-seedJump;
		}
	}
	//Something has gone wrong
	else {
		//Other Optimiser failure
		cout << prevDot << " prevdot // dotp: " << dotP << endl;
		optEl = -999.0;
		contin = false;
	}

}


