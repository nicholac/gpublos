//
//  simplekern.c
//  gputest1
//
//  Created by Chris Nicholas on 05/06/2015.
//  Copyright (c) 2015 Chris Nicholas. All rights reserved.
//


//Utility Functions
//
//void rk4Motion(position& outPos, const position& inPos, velocity& outVel, const velocity& inVel,
//                            acceleration& outAcc, const acceleration& inAcc, const double& timeStep, const double& mortSigma,
//                            const double& mortMass, const double& az, const double& el) {
//    /**
//     * RK4 Integrator
//     */
//    auto pos1 = position {0.0,0.0,0.0};
//    auto pos2 = position {0.0,0.0,0.0};
//    auto pos3 = position {0.0,0.0,0.0};
//    auto pos4 = position {0.0,0.0,0.0};
//    auto tayPos = position {0.0,0.0,0.0};
//    
//    auto vel1 = velocity {0.0,0.0,0.0};
//    auto vel2 = velocity {0.0,0.0,0.0};
//    auto vel3 = velocity {0.0,0.0,0.0};
//    auto vel4 = velocity {0.0,0.0,0.0};
//    auto tayVel = velocity {0.0,0.0,0.0};
//    
//    auto acc1 = acceleration {0.0,0.0,0.0};
//    auto acc2 = acceleration {0.0,0.0,0.0};
//    auto acc3 = acceleration {0.0,0.0,0.0};
//    auto acc4 = acceleration {0.0,0.0,0.0};
//    //auto tayAcc = acceleration {0.0,0.0,0.0};
//    
//    rk4Evaluate(pos1, inPos, vel1, inVel, acc1, inAcc, 0.0, mortSigma, mortMass, az, el);
//    rk4Evaluate(pos2, pos1, vel2, vel1, acc2, acc1, timeStep*0.5, mortSigma, mortMass, az, el);
//    rk4Evaluate(pos3, pos2, vel3, vel2, acc3, acc2, timeStep*0.5, mortSigma, mortMass, az, el);
//    rk4Evaluate(pos4, pos3, vel4, vel3, acc4, acc3, timeStep, mortSigma, mortMass, az, el);
//    
//    //Taylor Series Expansion
//    //outVel = (vel1+((vel2+vel3)*2.0)+vel4) * (1.0/6.0);
//    //outAcc = (acc1+((acc2+acc3)*2.0)+acc4) * (1.0/6.0);
//    //outPos = (pos1+((pos2+pos3)*2.0)+pos4) * (1.0/6.0);
//    
//    motionUtils::taylorExPos( pos1,  pos2,  pos3,  pos4,  tayPos);
//    motionUtils::taylorExVel( vel1,  vel2,  vel3,  vel4,  tayVel);
//    motionUtils::taylorExAcc( acc1,  acc2,  acc3,  acc4,  outAcc);
//    
//    //Calc final state
//    outPos = inPos + (tayVel * timeStep);
//    outVel = inVel + (outAcc * timeStep);
//    
//    
//    return;
//    
//}
//
//
//void taylorExVel( const float4 *vel1,  const float4 *vel2,  const float4 *vel3, const float4 *vel4, float4 *outVel) {
//    
//    outVel = (vel1+((vel2+vel3)*2.0)+vel4) * (1.0/6.0);
//    
//}
//
//void taylorExAcc( const float4 *acc1, const float4 *acc2, const float4 *acc3, const float4 *acc4, float4 *outAcc) {
//    
//    outAcc = (acc1+((acc2+acc3)*2.0)+acc4) * (1.0/6.0);
//    
//}
//
//void taylorExPos( const float4 *pos1, const float4 *pos2, const float4 *pos3, const float4 *pos4, float4 *outPos) {
//    
//    outPos = (pos1+((pos2+pos3)*2.0)+pos4) * (1.0/6.0);
//}
//
//

//void rk4Evaluate(float4 *outPos, float4 *inPos, float4 *outVel, float4 *inVel,
//                              float4 *outAcc, float4 *inAcc,  float *timeStep,  float *mortSigma,
//                               float *mortMass,  float *az,  float *el) {
//    /**
//     * RK4 Single Step Evaluate Function
//     float4 - s0 = x,
//              s1 = y,
//              s2 = z,
//              s3 = SPARE
//     */
//    
//    //Forces
//    __private double combVel;
//    __private double netForce;
//    __private double dragForce;
//    //Normalised Velocity
//    __private float4 normVel;
//    
//    //Drag Components
//    __private float4 drag;
//    
//    //Euler Velocity
//    outVel = inVel + (inAcc * timeStep);
//    //Euler Position
//    outPos = inPos + (outVel * timeStep) + ((inAcc * pow(timeStep, 2.0f))*0.5f);
//    
//    //Drag and accels
//    combVel = sqrt(pow(*outVel.x, 2.0f)+pow(outVel.y, 2.0f)+pow(outVel.z, 2.0f));
//    //motionUtils::drag(netForce, combVel, mortSigma, outPos.z);
//    dragForce = mortSigma*1.225f*pow(combVel, 2.0f);
//    //Normalise vector
//    normVel = outVel / combVel;
//    //Drag Components
//    drag = (normVel * netForce)*-1.0f;
//    //Add Gravity force
//    drag.z+=((mortMass*9.801f)*-1.0f);
//    //Acceleration components
//    outAcc = drag/mortMass;
//
//    return;
//}




__kernel void rangeChk(__global float8* demArr, float4 chkPos, float minDist, float maxDist)
    {
        //Reduces by Range
        int gid = get_global_id (0);

        //Work out the geo coords for this location
        __private float diffX, diffY, chkDist;
        __private float8 launchPos;
        __private int cnt;
        launchPos = demArr[gid];
        //Check Ranges
        diffX = fabs(chkPos.x-launchPos.x);
        diffY = fabs(chkPos.y-launchPos.y);
        chkDist = sqrt(pow(diffX, 2)+pow(diffY, 2));
        if (chkDist < minDist || chkDist > maxDist) {
            //Set the no run flag
            demArr[gid].s5 = 0.0f;
        }
        
    }

__kernel void testBasic(__global float* worldArr, __global float* outArr)
{
    //Mega Simple Addition - to test loop running
    int gid = get_global_id (0);
    __private float data;
    data = worldArr[gid];
    outArr[gid] = 2.0f+1.0f;
}


__kernel void thetaDistReduce(__global float* thetaArr, __global float* distArr, __global float8* demArr, int demIdx)
{
    //Calculate and record the mimumim distance from all shots - but on the device in the demArr
    int gid = get_global_id (0);
    __private float dist = distArr[gid];
    __private float theta = thetaArr[gid];
    if (demArr[demIdx].s5 == 1.0f) {
        //Good to process
        if (demArr[demIdx].s4 > dist){
            //Write dist and theta
            demArr[demIdx].s4 = dist;
            demArr[demIdx].s3 = theta;
            //printf("%f, %f, %f, %f\n", dist, theta, demArr[demIdx].s4, demArr[demIdx].s3);
        }
    }
    
}


__kernel void shotOpt(__global float* thetaArr, __global float* distArr, float4 tgtPos, float8 laPos, float d_timeStep)
{
    //Shot Kernel for Optimisation / Brute Forcing the theta shots (parrellelised based on all theta possibilities)
    
    int gid = get_global_id (0);
    
    //Check for run flag
    if (laPos.s5 == 1.0f){
    
        __private float4 outPos;
        outPos.s0 = 0.0f,
        outPos.s1 = 0.0f,
        outPos.s2 = 0.0f,
        outPos.s3 = 0.0f;
        __private float4 inPos;
        inPos.s0 = 0.0f,
        inPos.s1 = 0.0f,
        inPos.s2 = 0.0f,
        inPos.s3 = 0.0f;
        __private float4 initPos;
        initPos.s0 = 0.0f,
        initPos.s1 = 0.0f,
        initPos.s2 = 0.0f,
        initPos.s3 = 0.0f;
        __private float4 outVel;
        outVel.s0 = 0.0f,
        outVel.s1 = 0.0f,
        outVel.s2 = 0.0f,
        outVel.s3 = 0.0f;
        __private float4 inVel;
        inVel.s0 = 0.0f,
        inVel.s1 = 0.0f,
        inVel.s2 = 0.0f,
        inVel.s3 = 0.0f;
        __private float4 initVel;
        initVel.s0 = 0.0f,
        initVel.s1 = 0.0f,
        initVel.s2 = 0.0f,
        initVel.s3 = 0.0f;
        __private float4 outAcc;
        outAcc.s0 = 0.0f,
        outAcc.s1 = 0.0f,
        outAcc.s2 = 0.0f,
        outAcc.s3 = 0.0f;
        __private float4 inAcc;
        inAcc.s0 = 0.0f,
        inAcc.s1 = 0.0f,
        inAcc.s2 = 0.0f,
        inAcc.s3 = 0.0f;
        __private float4 initAcc;
        initAcc.s0 = 0.0f,
        initAcc.s1 = 0.0f,
        initAcc.s2 = 0.0f,
        initAcc.s3 = 0.0f;
        __private float4 vel1, vel2, vel3, vel4;
        __private float4 pos1, pos2, pos3, pos4;
        __private float4 acc1, acc2, acc3, acc4;
        __private float4 tayVel;
        tayVel.s0 = 0.0f,
        tayVel.s1 = 0.0f,
        tayVel.s2 = 0.0f,
        tayVel.s3 = 0.0f;
        __private float4 tayPos;
        tayPos.s0 = 0.0f,
        tayPos.s1 = 0.0f,
        tayPos.s2 = 0.0f,
        tayPos.s3 = 0.0f;
        
        __private float timeStep = d_timeStep;
        __private float mortSigma;
        __private float mortMass = 3.2f;
        __private float muzVel = 225.1f;
        __private float calibre = 0.081;
        __private float area;
        __private float dragCoef = 0.15;
        __private float pi = 3.14159265359;
        //Calc Mortsigma
        area = pi*(pow((0.5*calibre), 2));
        //Sigma
        mortSigma = dragCoef*area*0.5;
        
        //Get Theta from GID
        __private float el = thetaArr[gid];
        __private float elRad;
        //Convert to radians
        elRad = ((90.0-el)*0.01745329252);
        
        //Calc Az
        __private float dX = tgtPos.x-laPos.x;
        __private float dY = tgtPos.y-laPos.y;
        __private float az = atan2(dY, dX);
        if (az < 0.0f) {
            az+=(2.0f*pi);
        }
        
        //Forces
        __private float combVel = 0.0f;
        __private float netForce = 0.0f;
        __private float dragForce = 0.0f;
        //Normalised Velocity
        __private float4 normVel;
        normVel.s0 = 0.0f;
        normVel.s1 = 0.0f;
        normVel.s2 = 0.0f;
        normVel.s3 = 0.0f;
        
        //Drag Components
        __private float4 drag;
        drag.s0 = 0.0f;
        drag.s1 = 0.0f;
        drag.s2 = 0.0f;
        drag.s3 = 0.0f;
        
        //Distance Utils
        __private float diffX, diffY, dist;
        
        //Initialise Velocity
        initVel.x = muzVel * sin(elRad) * cos(az);
        initVel.y = muzVel * sin(elRad) * sin(az);
        initVel.z = muzVel * cos(elRad);
        
        //Setup the launch Position
        //inPos is used in the loop as currPos
        //Specifying here to convert from float8 to float4
        inPos.x = laPos.x;
        inPos.y = laPos.y;
        inPos.z = laPos.z;
        inVel = initVel;
        
        __private int cnt = 0;
        while (inPos.z > -1.0f){
            cnt++;
            //Eval 1
            
            //Euler Velocity
            vel1 = inVel + (inAcc * 0.0f);
            //Euler Position
            pos1 = inPos + (vel1 * 0.0f) + ((inAcc * 0.0f)*0.5f);
            
            //Drag and accels
            combVel = sqrt(pow(vel1.x, 2)+pow(vel1.y, 2)+pow(vel1.z, 2));
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.z);
            dragForce = mortSigma*1.225f*powr(combVel, 2);
            //Normalise vector
            normVel = vel1 / combVel;
            //Drag Components
            drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            drag.z+=((mortMass*9.801f)*-1.0f);
            //Acceleration components
            acc1 = drag/mortMass;

            //Eval 2
            //Euler Velocity
            vel2 = vel1 + (acc1 * (timeStep*0.5f));
            //Euler Position
            pos2 = pos1 + (vel2 * (timeStep*0.5f)) + ((acc1 * pow((timeStep*0.5f), 2))*0.5f);

            //Drag and accels
            combVel = sqrt(pow(vel1.x, 2)+pow(vel1.y, 2)+pow(vel1.z, 2));
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.z);
            dragForce = mortSigma*1.225f*pow(combVel, 2.0f);
            //Normalise vector
            normVel = vel2 / combVel;
            //Drag Components
            drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            drag.z+=((mortMass*9.801f)*-1.0f);
            //Acceleration components
            acc2 = drag/mortMass;

            //Eval 3
            //Euler Velocity
            vel3 = vel2 + (acc2 * (timeStep*0.5f));
            //Euler Position
            pos3 = pos2 + (vel3 * (timeStep*0.5f)) + ((acc2 * pow((timeStep*0.5f), 2))*0.5f);
            
            //Drag and accels
            combVel = sqrt(pow(vel3.x, 2)+pow(vel3.y, 2)+pow(vel3.z, 2));
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.z);
            dragForce = mortSigma*1.225f*powr(combVel, 2);
            //Normalise vector
            normVel = vel3 / combVel;
            //Drag Components
            drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            drag.z+=((mortMass*9.801f)*-1.0f);
            //Acceleration components
            acc3 = drag/mortMass;

            //Eval 4
            //Euler Velocity
            vel4 = vel3 + (acc3 * timeStep);
            //Euler Position
            pos4 = pos3 + (vel4 * timeStep) + ((acc3 * pow(timeStep, 2))*0.5f);
            
            //Drag and accels
            combVel = sqrt(pow(vel4.x, 2)+pow(vel4.y, 2)+pow(vel4.z, 2));
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.z);
            dragForce = mortSigma*1.225f*pow(combVel, 2);
            //Normalise vector
            normVel = vel4 / combVel;
            //Drag Components
            drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            drag.z+=((mortMass*9.801f)*-1.0f);
            //Acceleration components
            acc4 = drag/mortMass;
            
            //Taylor Expansion
            tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (1.0f/6.0f);
            inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (1.0f/6.0f);
            tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (1.0f/6.0f);
            
            //Swap ready for next iteration
            inPos = inPos + (tayVel * timeStep);
            inVel = inVel + (inAcc * timeStep);
        
        }
        //Calculate the distance to tgt
        diffX = fabs(inPos.x-tgtPos.x);
        diffY = fabs(inPos.y-tgtPos.y);
        dist = sqrt(powr(diffX, 2)+powr(diffY, 2));
        distArr[gid] = dist;
        //printf("%f, %f, %f, %i\n", inPos.x, inPos.y, inPos.z, cnt);
    } else {
        //no run flag found
    }
}


__kernel void shotISect(__global float* thetaArr, __global float* distArr, float4 tgtPos, float8 laPos, float d_timeStep)
{
    //Shot Kernel for intersection (parrellelised based on DEM)
    
    int gid = get_global_id (0);
    
    __private float4 outPos;
    outPos.s0 = 0.0f,
    outPos.s1 = 0.0f,
    outPos.s2 = 0.0f,
    outPos.s3 = 0.0f;
    __private float4 inPos;
    inPos.s0 = 0.0f,
    inPos.s1 = 0.0f,
    inPos.s2 = 0.0f,
    inPos.s3 = 0.0f;
    __private float4 initPos;
    initPos.s0 = 0.0f,
    initPos.s1 = 0.0f,
    initPos.s2 = 0.0f,
    initPos.s3 = 0.0f;
    __private float4 outVel;
    outVel.s0 = 0.0f,
    outVel.s1 = 0.0f,
    outVel.s2 = 0.0f,
    outVel.s3 = 0.0f;
    __private float4 inVel;
    inVel.s0 = 0.0f,
    inVel.s1 = 0.0f,
    inVel.s2 = 0.0f,
    inVel.s3 = 0.0f;
    __private float4 initVel;
    initVel.s0 = 0.0f,
    initVel.s1 = 0.0f,
    initVel.s2 = 0.0f,
    initVel.s3 = 0.0f;
    __private float4 outAcc;
    outAcc.s0 = 0.0f,
    outAcc.s1 = 0.0f,
    outAcc.s2 = 0.0f,
    outAcc.s3 = 0.0f;
    __private float4 inAcc;
    inAcc.s0 = 0.0f,
    inAcc.s1 = 0.0f,
    inAcc.s2 = 0.0f,
    inAcc.s3 = 0.0f;
    __private float4 initAcc;
    initAcc.s0 = 0.0f,
    initAcc.s1 = 0.0f,
    initAcc.s2 = 0.0f,
    initAcc.s3 = 0.0f;
    __private float4 vel1, vel2, vel3, vel4;
    __private float4 pos1, pos2, pos3, pos4;
    __private float4 acc1, acc2, acc3, acc4;
    __private float4 tayVel;
    tayVel.s0 = 0.0f,
    tayVel.s1 = 0.0f,
    tayVel.s2 = 0.0f,
    tayVel.s3 = 0.0f;
    __private float4 tayPos;
    tayPos.s0 = 0.0f,
    tayPos.s1 = 0.0f,
    tayPos.s2 = 0.0f,
    tayPos.s3 = 0.0f;
    
    __private float timeStep = d_timeStep;
    __private float mortSigma;
    __private float mortMass = 3.2f;
    __private float muzVel = 225.1f;
    __private float calibre = 0.081;
    __private float area;
    __private float dragCoef = 0.15;
    __private float pi = 3.14159265359;
    //Calc Mortsigma
    area = pi*(pow((0.5*calibre), 2));
    //Sigma
    mortSigma = dragCoef*area*0.5;
    
    //Get Theta from GID
    __private float el = thetaArr[gid];
    __private float elRad;
    //Convert to radians
    elRad = ((90.0-el)*0.01745329252);
    
    //Calc Az
    __private float dX = tgtPos.x-laPos.x;
    __private float dY = tgtPos.y-laPos.y;
    __private float az = atan2(dY, dX);
    if (az < 0.0f) {
        az+=(2.0f*pi);
    }
    
    //Forces
    __private float combVel = 0.0f;
    __private float netForce = 0.0f;
    __private float dragForce = 0.0f;
    //Normalised Velocity
    __private float4 normVel;
    normVel.s0 = 0.0f;
    normVel.s1 = 0.0f;
    normVel.s2 = 0.0f;
    normVel.s3 = 0.0f;
    
    //Drag Components
    __private float4 drag;
    drag.s0 = 0.0f;
    drag.s1 = 0.0f;
    drag.s2 = 0.0f;
    drag.s3 = 0.0f;
    
    //Distance Utils
    __private float diffX, diffY, dist;
    
    //Initialise Velocity
    initVel.x = muzVel * sin(elRad) * cos(az);
    initVel.y = muzVel * sin(elRad) * sin(az);
    initVel.z = muzVel * cos(elRad);
    
    //Setup the launch Position
    //inPos is used in the loop as currPos
    inPos.x = laPos.x;
    inPos.y = laPos.y;
    inPos.z = laPos.z;
    inVel = initVel;
    
    __private int cnt = 0;
    while (inPos.z > -1.0f){
        cnt++;
        //Eval 1
        
        //Euler Velocity
        vel1 = inVel + (inAcc * 0.0f);
        //Euler Position
        pos1 = inPos + (vel1 * 0.0f) + ((inAcc * 0.0f)*0.5f);
        
        //Drag and accels
        combVel = sqrt(pow(vel1.x, 2)+pow(vel1.y, 2)+pow(vel1.z, 2));
        //motionUtils::drag(netForce, combVel, mortSigma, outPos.z);
        dragForce = mortSigma*1.225f*powr(combVel, 2);
        //Normalise vector
        normVel = vel1 / combVel;
        //Drag Components
        drag = (normVel * dragForce)*-1.0f;
        //Add Gravity force
        drag.z+=((mortMass*9.801f)*-1.0f);
        //Acceleration components
        acc1 = drag/mortMass;
        
        //Eval 2
        //Euler Velocity
        vel2 = vel1 + (acc1 * (timeStep*0.5f));
        //Euler Position
        pos2 = pos1 + (vel2 * (timeStep*0.5f)) + ((acc1 * pow((timeStep*0.5f), 2))*0.5f);
        
        //Drag and accels
        combVel = sqrt(pow(vel1.x, 2)+pow(vel1.y, 2)+pow(vel1.z, 2));
        //motionUtils::drag(netForce, combVel, mortSigma, outPos.z);
        dragForce = mortSigma*1.225f*pow(combVel, 2.0f);
        //Normalise vector
        normVel = vel2 / combVel;
        //Drag Components
        drag = (normVel * dragForce)*-1.0f;
        //Add Gravity force
        drag.z+=((mortMass*9.801f)*-1.0f);
        //Acceleration components
        acc2 = drag/mortMass;
        
        //Eval 3
        //Euler Velocity
        vel3 = vel2 + (acc2 * (timeStep*0.5f));
        //Euler Position
        pos3 = pos2 + (vel3 * (timeStep*0.5f)) + ((acc2 * pow((timeStep*0.5f), 2))*0.5f);
        
        //Drag and accels
        combVel = sqrt(pow(vel3.x, 2)+pow(vel3.y, 2)+pow(vel3.z, 2));
        //motionUtils::drag(netForce, combVel, mortSigma, outPos.z);
        dragForce = mortSigma*1.225f*powr(combVel, 2);
        //Normalise vector
        normVel = vel3 / combVel;
        //Drag Components
        drag = (normVel * dragForce)*-1.0f;
        //Add Gravity force
        drag.z+=((mortMass*9.801f)*-1.0f);
        //Acceleration components
        acc3 = drag/mortMass;
        
        //Eval 4
        //Euler Velocity
        vel4 = vel3 + (acc3 * timeStep);
        //Euler Position
        pos4 = pos3 + (vel4 * timeStep) + ((acc3 * pow(timeStep, 2))*0.5f);
        
        //Drag and accels
        combVel = sqrt(pow(vel4.x, 2)+pow(vel4.y, 2)+pow(vel4.z, 2));
        //motionUtils::drag(netForce, combVel, mortSigma, outPos.z);
        dragForce = mortSigma*1.225f*pow(combVel, 2);
        //Normalise vector
        normVel = vel4 / combVel;
        //Drag Components
        drag = (normVel * dragForce)*-1.0f;
        //Add Gravity force
        drag.z+=((mortMass*9.801f)*-1.0f);
        //Acceleration components
        acc4 = drag/mortMass;
        
        //Taylor Expansion
        tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (1.0f/6.0f);
        inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (1.0f/6.0f);
        tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (1.0f/6.0f);
        
        //Swap ready for next iteration
        inPos = inPos + (tayVel * timeStep);
        inVel = inVel + (inAcc * timeStep);
        
    }
    //Calculate the distance to tgt
    diffX = fabs(inPos.x-tgtPos.x);
    diffY = fabs(inPos.y-tgtPos.y);
    dist = sqrt(powr(diffX, 2)+powr(diffY, 2));
    distArr[gid] = dist;
    //printf("%f, %f, %f, %i\n", inPos.x, inPos.y, inPos.z, cnt);
    
}





//__kernel void test2d(__global float2* demArr)
//{ //Just sets up the optEl array
//    int gid = get_global_id(0);
//    optElArr[i]=85.0f;
//}
//
//__kernel void setupOptEl(__global float* optElArr)
//    { //Just sets up the optEl array
//        int i = get_global_id(0);
//        optElArr[i]=85.0f;
//    }
//
//__kernel void setupSeedJump(__global float* seedJumpArr)
//    { //Just sets up the optEl array
//        int i = get_global_id(0);
//        seedJumpArr[i]=5.0f;
//    }
//
//__kernel void missDist(__global float2* missPosArr,
//                       float tgtX, float tgtY, float tgtSize)
//
//    {//Distance between miss position and tgt
//        int gid = get_global_id (0);
//        __private float dist, diffX, diffY;
//        __private int proc;
//        diffX = fabs(tgtX-missPosArr[gid].x);
//        diffY = fabs(tgtY-missPosArr[gid].y);
//        dist = sqrt(powr(diffX, 2)+powr(diffY, 2));
//    }
//
//__kernel void checkOpt(__global float* optElArr, __global float* seedJumpArr,
//                       __global float* procIdx)
//
//    {//Check array val and update
//        int gid = get_global_id(0);
//        __private int proc;
//        if (optElArr[gid] < 0.0f){
//            procIdx[gid]=0;
//        }
//        if (seedJumpArr[gid] < 0.0005f){
//            procIdx[gid]=0;
//        }
//    }

//__kernel void optimise(float rasterXSize, float rasterYSize, float gtrans0,
//                      float gtrans1, float gtrans2, float gtrans3, float gtrans4,
//                      float gtrans5, float tgtX, float tgtY,
//                      __global float2* missPosArr,
//                      __global float* dotPArr,
//                      __global float* seedJumpArr,
//                      __global float* optElArr)
//
//    { //Works out dot product and saves related seedjump - basically the optimise kernel
//        int gid = get_global_id (0);
//        __private float vec1x, vec1y, vec2x, vec2y, out, dotP, el, seedJump, tmp, tmp2, pxX, pxY;
//        __private float launchX, launchY;
//        __private float2 v1, v2;
//        if (procIdx[gid] > 0){
//            //Work out pixel location and launch from 1d index
//            tmp = modf(float(gid)/rasterXSize, &pxY);
//            tmp2 = modf(tmp * rasterXSize, &pxX);
//            launchX = gtrans0+(pxX*gtrans1)+pxY*gtrans2;
//            launchY = gtrans3+(pxX*gtrans4)+pxY*gtrans5;
//            //Do Dot Prod
//            vec1x = launchX - tgtX;
//            vec1y = launchY - tgtY;
//            vec2x = missPosArr[gid].x - tgtX;
//            vec2y = missPosArr[gid].y - tgtY;
//            v1 = (vec1x, vec1y);
//            v2 = (vec2x, vec2y);
//            //opencl dot - this is the new value - old is in arr
//            dotP = dot(v1, v2);
//            //Check and Move
//            if (dotP < 0.0f) {
//                if (dotPArr[gid] < 0.0f) {
//                    //Too far and prev was too far (negative dot prod between the vectors)
//                    //Increase elevation, no change in seed
//                    el = optElArr[gid]+seedJumpArr[gid];
//                }
//                if (dotPArr[gid] > 0.0f) {
//                    //Too far and prev too close:
//                    //alter seedjump and increase el
//                    seedJump=seedJumpArr[gid]/4.0f;
//                    el = optElArr[gid]+seedJump;
//                }
//            }
//            else if (dotP > 0.0f) {
//                if (dotPArr[gid] > 0.0f) {
//                    //Too close and prev too close
//                    //decrease elevation, no change seed
//                    el = optElArr[gid]-seedJumpArr[gid];
//                }
//                if (dotPArr[gid] < 0.0f) {
//                    //Too close and prev too far
//                    //Change seed, decrease el
//                    seedJump = seedJumpArr[gid]/4.0f;
//                    el = optElArr[gid]-seedJump;
//                }
//            }
//            //Set Vals
//            seedJumpArr[gid]=seedJump;
//            optElArr[gid] = el;
//        }
//    }


//__kernel void shotStep(__global float8* trjArr, float4 totTime, float4 timestep, float4 tgtPos, float mortMass,
//                   float muzVel, float mortSigma, float gForce,
//                   __global float2* missPosArr)
//
//{ //Single step of the recshot func
//    int gid = get_global_id (0);
//    //Get curr position from world Arr
//    __private float8 currPos;
//    __private float8 prevPos;
//    //OutPos is the new position
//    __private float8 outPos;
//    currPos = worldArr[gid];
//    
//    //Record Previous Position
//    prevPos = currPos;
//    //Total Time
//    totTime+=timeStep;
//    //Move Position
//    currPos.x = prevPos.x + ((vel.x * timeStep) + (0.5f*acc.x*(powr(timeStep, 2))));
//    currPos.y = prevPos.y + ((vel.y * timeStep) + (0.5f*acc.y*(powr(timeStep, 2))));
//    currPos.z = prevPos.z + ((vel.z * timeStep) + (0.5f*acc.z*(powr(timeStep, 2))));
//    //Check if its passed zenith
//    if (floorMove == false) {
//        if (prevPos.z > currPos.z) {
//            floorZ = ptgtPos.z;
//            floorMove = true;
//            //Record zenith
//            launchPos.s7 = currPos.z;
//            //printf("%f, %f, %f, %f, %f, %i, %f, %f\n", launchZ, floorZ, tgtZ, currPosZ, combV, lpCnt, prevPosZ, t);
//        }
//    }
//    lpCnt+=1;
//    //Drag Forces
//    //Combined xyz velocity
//    combV = sqrt(powr(vel.x, 2)+powr(vel.y, 2)+powr(vel.z, 2));
//    //Drag Force xyz
//    dragForce = mortSigma*powr(combV, 2);
//    if (combV > 0.0f){
//        dragForce*=-1.0f;
//    }
//    //This actually gives the combined drag - not a velocity but reusing the class
//    force.x = dragForce * sin(elRad) * cos(az);
//    force.y = dragForce * sin(elRad) * sin(az);
//    force.z = dragForce * cos(elRad);
//    //Gravity Force
//    force.z-=gForce;
//    //Verlet Integration - TODO REMOVE THIS
//    avgAcc.x = 0.5*((force.x/mortMass)+acc.x);
//    avgAcc.y = 0.5*((force.y/mortMass)+acc.y);
//    avgAcc.z = 0.5*((force.z/mortMass)+acc.z);
//    //Alter Velocity - TODO USE FMA HERE
//    vel.x+=avgAcc.x*timeStep;
//    vel.y+=avgAcc.y*timeStep;
//    vel.z+=avgAcc.z*timeStep;
//    //Reset Vertlet
//    acc.x = avgAcc.x;
//    acc.y = avgAcc.y;
//    acc.z = avgAcc.z;
//
//    
//}
//
//__kernel void shot(__global float8* worldArr, float4 tgtPos, float mortMass,
//                        float muzVel, float mortSigma, float gForce)
//
//{
//    //Runs single shot - returns last position
//    //float8(x, y ,z ,optEl, dotp, prevdot, dist, SPARE)
//    int gid = get_global_id (0);
//    __private float diffX, diffY;
//    __private float az;
//    __private float radian;
//    __private float4 vel;
//    __private float combV;
//    __private float4 acc;
//    __private float4 avgAcc;
//    __private float4 force;
//    __private float dragForce;
//    __private float4 currPos;
//    __private float4 prevPos;
//    __private float floorZ;
//    __private float elRad;
//    __private float pi;
//    __private float timeStep;
//    __private float t;
//    __private bool floorMove;
//    __private bool zenFail;
//    __private int lpCnt;
//    __private int resVal;
//    //Landing Point
//    __private float2 missPos;
//    missPos = (float2)(0.0f, 0.0f);
//    lpCnt=0;
//    floorMove = false;
//    t = 0.0f;
//    zenFail = false;
//    //Get float8 from input
//    //float8(x, y ,z ,optEl, dotp, prevdot, dist, SPARE)
//    //pos.x, pos.y, pos.z
//    //pos.s3 - optel
//    //pos.s4 - dotp
//    //pos.s5 - prevDot
//    //pos.s6 - dist land 2 tgt
//    __private float8 launchPos;
//    launchPos = worldArr[gid];
//    __private float4 ptgtPos;
//    ptgtPos = tgtPos;
//
//    //Adding obs height to launch Z and tgtZ
//    launchPos.z+=20.0f;
//    ptgtPos.z+=20.0f;
//    timeStep = 0.01f;
//    pi = 3.1415926f;
//    //Calc Az
//    diffX = fabs(ptgtPos.x-launchPos.x);
//    diffY = fabs(ptgtPos.y-launchPos.y);
//    az = atan2(diffY, diffX);
//    if (az < 0.0f) {
//        az+=(2.0f*pi);
//    }
//    //Deg2rad
//    elRad = (90.0f-launchPos.s3) * (pi/180.0f);
//    //Vel components
//    vel.x = muzVel * sin(elRad) * cos(az);
//    vel.y = muzVel * sin(elRad) * sin(az);
//    vel.z = muzVel * cos(elRad);
//    //Acc Components
//    acc.x = 0.0f;
//    acc.y = 0.0f;
//    acc.z = 0.0f;
//    //Setup Floor
//    floorZ = -99999.0f;
//    //Reset Forces
//    force = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
//    dragForce=0.0f;
//    //Setup currpos
//    currPos.x = launchPos.x;
//    currPos.y = launchPos.y;
//    currPos.z = launchPos.z;
//    //Run Shot
//    while (currPos.z > floorZ) { //While not below floor
//        //Record Previous Position
//        prevPos = currPos;
//        //Total Time
//        t+=timeStep;
//        //Move Position
//        currPos.x = prevPos.x + ((vel.x * timeStep) + (0.5f*acc.x*(powr(timeStep, 2))));
//        currPos.y = prevPos.y + ((vel.y * timeStep) + (0.5f*acc.y*(powr(timeStep, 2))));
//        currPos.z = prevPos.z + ((vel.z * timeStep) + (0.5f*acc.z*(powr(timeStep, 2))));
//        //Check if its passed zenith
//        if (floorMove == false) {
//            if (prevPos.z > currPos.z) {
//                floorZ = ptgtPos.z;
//                floorMove = true;
//                //Record zenith
//                launchPos.s7 = currPos.z;
//                //printf("%f, %f, %f, %f, %f, %i, %f, %f\n", launchZ, floorZ, tgtZ, currPosZ, combV, lpCnt, prevPosZ, t);
//            }
//        }
//        lpCnt+=1;
//        //Drag Forces
//        //Combined xyz velocity
//        combV = sqrt(powr(vel.x, 2)+powr(vel.y, 2)+powr(vel.z, 2));
//        //Drag Force xyz
//        dragForce = mortSigma*powr(combV, 2);
//        if (combV > 0.0f){
//            dragForce*=-1.0f;
//        }
//        //This actually gives the combined drag - not a velocity but reusing the class
//        force.x = dragForce * sin(elRad) * cos(az);
//        force.y = dragForce * sin(elRad) * sin(az);
//        force.z = dragForce * cos(elRad);
//        //Gravity Force
//        force.z-=gForce;
//        //Verlet Integration - TODO REMOVE THIS
//        avgAcc.x = 0.5*((force.x/mortMass)+acc.x);
//        avgAcc.y = 0.5*((force.y/mortMass)+acc.y);
//        avgAcc.z = 0.5*((force.z/mortMass)+acc.z);
//        //Alter Velocity - TODO USE FMA HERE
//        vel.x+=avgAcc.x*timeStep;
//        vel.y+=avgAcc.y*timeStep;
//        vel.z+=avgAcc.z*timeStep;
//        //Reset Vertlet
//        acc.x = avgAcc.x;
//        acc.y = avgAcc.y;
//        acc.z = avgAcc.z;
//    } //End while not below floor loop
//    //Record Miss point in 2float
//    //Record
//    missPosArr[gid] = (float2)(currPos.x, currPos.y);
//}
//
//
//__kernel void intersect(__global float* worldArr,  __global int* result)
//Run trajectories & profiles
//Input - demArr, procIdx, optEls, weapon Params
//Process - intersection between lines
//Output - int results for blos and tlos
//{
//}

//
//__kernel void intersect_segment(__global float* worldArr, __global float* trjArr, float linLen, __global int* result)
//
//{
//    //This function presumes we use the kernel to process the intersect between one lines segment
//    //So trj - line - line segment -< kernels do this bit
//    //Function below (intersect_line) uses one kernel per LINE (so one level back
//    /*
//    int i = get_global_id (0);
//    for (int i=0; i < 100000; i++){
//        int p = i*1000.0;
//    }
//    group_result [i] = input[i]+10.0;
//     */
//    __local float delta;
//    __local float iSectX;
//    __local float iSectY;
//    __local float x1;
//    x1 = 0.0;
//    //TODO: pull lin len in or calculate
//    //float x2 = linLen;
//    __local float x3;
//    x3 = 0.0;
//    //float linLen = linLen;
//    //y1 = inWorld.at(i)->z;
//    //y2 = inWorld.at(i+1)->z;
//    //y3 = inTrj.at(i)->z;
//    //y4 = inTrj.at(i+1)->z;
//    
//    
//    //Get the array index for this kernel from the global memory
//    size_t i = get_global_id (0);
//    //Get the data associated with this index
//    __local float y1;
//    __local float y2;
//    __local float y3;
//    __local float y4;
//    y1 = worldArr[i];
//    y2 = worldArr[i+1];
//    y3 = trjArr[i];
//    y4 = trjArr[i+1];
//    
//    
//    //Check Parallel
//    delta = ((x1-linLen)*(y3-y4))-((y1-y2)*(x3-linLen));
//    if (delta == 0.0) {
//        result[i] = 0;
//    }
//    else {
//        //Check intersect
//        iSectX = ((((x1*y2)-(y1*linLen))*(x3-linLen))-((x1-linLen)*((x3*y4)-(y3*linLen))))/delta;
//        iSectY = ((((x1*y2)-(y1*linLen))*(y3-y4))-((y1-y2)*((x3*y4)-(y3*linLen))))/delta;
//        if (iSectX >= x1) {
//            if (iSectX <= linLen) {
//                if (iSectY >= min(y3, y4)) {
//                    if (iSectY <= max(y3, y4)) {
//                        //Intersect
//                        result[i] = 1;
//                    }
//                }
//                
//            }
//        }
//        //Otherwise not
//        result[i] = 0;
//    }
//    
//    
//}
//

