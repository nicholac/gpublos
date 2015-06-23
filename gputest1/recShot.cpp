//
//  recShot.cpp
//  gputest1
//
//  Created by Chris Nicholas on 11/06/2015.
//  Copyright (c) 2015 Chris Nicholas. All rights reserved.
//

#include "recShot.h"
#include <vector>
#include <memory> // this is for the safe pointer wrapper
#include <position.h>
#include <velocity.h>
#include <acceleration.h>
#include <velocityUtils.h>
#include <azUtils.h>
#include <motionUtils.h>
#include <distUtils.h>
#include <worldParams.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;


std::vector<std::shared_ptr<position>> recShot::launchAir(const std::shared_ptr<position>& launchPos,
                                                          const std::shared_ptr<position>& tgtPos,
                                                          const double& v0, const double& elDeg,
                                                          const double& timeStep, const double& mortSigma, const double& mortMass) {
    /* *
     * Air Resistance Single Shot
     * param: elDeg - elevation in degrees from horizontal
     * Returns: vector containing shared ptr<positions>
     */
    
    //int zenith = 0;
    double floorZ = -999999.0;
    bool floorMove = false;
    bool zenFail = false;
    
    //Azimuth of shot - radians
    double az = azUtils::calcAz(launchPos->x, launchPos->y, tgtPos->x, tgtPos->y);
    
    //elevation of shot - input as degrees from horizontal, change to degrees from vertical
    double el = azUtils::deg2Rad((90.0-elDeg));
    //TODO:Make velocity a shared ptr
    //Get velocity components at Launch
    velocity vel = velocityUtils::calcInitVel(v0, az, el);
    
    //Initial time
    double t = 0.0;
    
    //Initial Accelerations
    acceleration acc = acceleration();
    acceleration avgAcc = acceleration();
    acc.x = 0.0;
    acc.y = 0.0;
    acc.z = 0.0;
    
    //Init Gravity Force
    float gForce = mortMass*gravity;
    
    //Forces
    double forceX;
    double forceY;
    double forceZ;
    double combForce;
    velocity combDrag();
    
    //Init starting Position and check position for loop
    auto currPos = std::make_shared<position>(position {0.0,0.0,0.0});
    currPos->x = launchPos->x;
    currPos->y = launchPos->y;
    currPos->z = launchPos->z;
    
    auto prevPos = std::make_shared<position>(position {0.0,0.0,0.0});
    
    //Out positions
    auto outPosVec = std::vector<std::shared_ptr<position>>();
    
    //Push start position
    recPos(outPosVec, launchPos);
    
    //Start trajectory loop
    while (currPos->z > floorZ) { //While not below floor
        //Record Previous Position
        prevPos->x = currPos->x;
        prevPos->y = currPos->y;
        prevPos->z = currPos->z;
        //Total Time
        t = t+timeStep;
        //Move Position
        motionUtils::motion(currPos, prevPos, vel, acc, timeStep);
        //Check if its passed zenith
        if (floorMove == false) {
            zenCheck(currPos, tgtPos, prevPos, floorMove, floorZ, zenFail);
            //Check if the zenith failed to reach tgt height
            if (zenFail == true){
                //Return null result
                recShot::recPos(outPosVec, std::make_shared<position>(position {-99,-99,-99}));
                return outPosVec;
            }
        }
        //TODO:  Check if its outside the bounds of DEM
        //Record the position
        recPos(outPosVec, currPos);
        //Drag Forces
        //Combined xyz velocity
        double combV = sqrt(pow(vel.x, 2.0)+pow(vel.y, 2.0)+pow(vel.z, 2.0));
        //Drag Force xyz
        motionUtils::drag(combForce, combV, mortSigma);
        //Todo: improve memory copying here in velocity
        //This actually gives the combined drag - not a velocity but reusing the class
        velocity combDrag = velocityUtils::calcInitVel(combForce, az, el);
        forceX = combDrag.x;
        forceY = combDrag.y;
        forceZ = combDrag.z;
        
        //Gravity Force
        forceZ-=gForce;
        //Verlet Integration
        avgAcc.x = 0.5*((forceX/mortMass)+acc.x);
        avgAcc.y = 0.5*((forceY/mortMass)+acc.y);
        avgAcc.z = 0.5*((forceZ/mortMass)+acc.z);
        //Alter Velocity
        vel.x+=avgAcc.x*timeStep;
        vel.y+=avgAcc.y*timeStep;
        vel.z+=avgAcc.z*timeStep;
        //Reset Verlet Accelerations
        acc.x = avgAcc.x;
        acc.y = avgAcc.y;
        acc.z = avgAcc.z;
        
    } //End while not below floor loop
    //Return the trajectory
    return outPosVec;
}

void recShot::recPos(std::vector<std::shared_ptr<position>>& vec, const std::shared_ptr<position>& pos) {
    /*
     //Make a safe shared position - copied because x,y,z will be re-used
     //Push it to vector at same time - referenced
     //TODO: THis is likely not the best way to do this
     */
    vec.push_back(std::make_shared<position>( position {pos->x, pos->y, pos->z}));
}

void recShot::zenCheck(const std::shared_ptr<position>& currPos, const std::shared_ptr<position>& tgtPos,
                       std::shared_ptr<position>& prevPos, bool& floorMove, double& floorZ, bool& zenFail) {
    //Check for zenith passing
    if (prevPos->z > currPos->z) {
        //Passed Zenith - Check if we reached tgt height
        if (prevPos->z < tgtPos->z){
            //Didnt reach the tgt height
            zenFail = true;
        }
        floorZ = tgtPos->z;
        floorMove = true;
    }
}