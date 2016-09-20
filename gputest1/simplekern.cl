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
//    combVel = sqrt(pow(*outVel.s0, 2.0f)+pow(outVel.s1, 2.0f)+pow(outVel.s2, 2.0f));
//    //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
//    dragForce = mortSigma*1.225f*pow(combVel, 2.0f);
//    //Normalise vector
//    normVel = outVel / combVel;
//    //Drag Components
//    drag = (normVel * netForce)*-1.0f;
//    //Add Gravity force
//    drag.s2+=((mortMass*9.801f)*-1.0f);
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
        launchPos = demArr[gid];
        //Check Ranges
        diffX = fabs(chkPos.s0-launchPos.s0);
        diffY = fabs(chkPos.s1-launchPos.s1);
        chkDist = half_sqrt(diffX*diffX+diffY*diffY);
        if (chkDist < minDist || chkDist > maxDist) {
            //Set the no run flag
            demArr[gid].s5 = 0.0f;
            //printf("%f, %f, %f\n", chkDist, minDist, maxDist);
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
    
        float4 outPos;
        outPos.s0 = 0.0f,
        outPos.s1 = 0.0f,
        outPos.s2 = 0.0f,
        outPos.s3 = 0.0f;
         float4 inPos;
        inPos.s0 = 0.0f,
        inPos.s1 = 0.0f,
        inPos.s2 = 0.0f,
        inPos.s3 = 0.0f;
         float4 initPos;
        initPos.s0 = 0.0f,
        initPos.s1 = 0.0f,
        initPos.s2 = 0.0f,
        initPos.s3 = 0.0f;
         float4 outVel;
        outVel.s0 = 0.0f,
        outVel.s1 = 0.0f,
        outVel.s2 = 0.0f,
        outVel.s3 = 0.0f;
         float4 inVel;
        inVel.s0 = 0.0f,
        inVel.s1 = 0.0f,
        inVel.s2 = 0.0f,
        inVel.s3 = 0.0f;
         float4 initVel;
        initVel.s0 = 0.0f,
        initVel.s1 = 0.0f,
        initVel.s2 = 0.0f,
        initVel.s3 = 0.0f;
         float4 outAcc;
        outAcc.s0 = 0.0f,
        outAcc.s1 = 0.0f,
        outAcc.s2 = 0.0f,
        outAcc.s3 = 0.0f;
         float4 inAcc;
        inAcc.s0 = 0.0f,
        inAcc.s1 = 0.0f,
        inAcc.s2 = 0.0f,
        inAcc.s3 = 0.0f;
         float4 initAcc;
        initAcc.s0 = 0.0f,
        initAcc.s1 = 0.0f,
        initAcc.s2 = 0.0f,
        initAcc.s3 = 0.0f;
         float4 vel1, vel2, vel3, vel4;
         float4 pos1, pos2, pos3, pos4;
         float4 acc1, acc2, acc3, acc4;
         float4 tayVel;
        tayVel.s0 = 0.0f,
        tayVel.s1 = 0.0f,
        tayVel.s2 = 0.0f,
        tayVel.s3 = 0.0f;
         float4 tayPos;
        tayPos.s0 = 0.0f,
        tayPos.s1 = 0.0f,
        tayPos.s2 = 0.0f,
        tayPos.s3 = 0.0f;
        
         float timeStep = d_timeStep;
         float mortSigma;
         float mortMass = 3.2f;
         float muzVel = 225.1f;
         float calibre = 0.081f;
         float area;
         float dragCoef = 0.15f;
         float pi = 3.14159265359f;
        //Calc Mortsigma
        area = pi*(pow((0.5f*calibre), 2));
        //Sigma
        mortSigma = dragCoef*area*0.5f;
        
        //Get Theta from GID
         float el = thetaArr[gid];
         float elRad;
        //Convert to radians
        elRad = ((90.0f-el)*0.01745329252f);
        
        //Calc Az
         float dX = tgtPos.s0-laPos.s0;
         float dY = tgtPos.s1-laPos.s1;
         float az = atan2(dY, dX);
        if (az < 0.0f) {
            az+=(2.0f*pi);
        }
        
        //Forces
         float combVel = 0.0f;
         float dragForce = 0.0f;
        //Normalised Velocity
         float4 normVel;
        normVel.s0 = 0.0f;
        normVel.s1 = 0.0f;
        normVel.s2 = 0.0f;
        normVel.s3 = 0.0f;
        
        //Drag Components
         float4 drag;
        drag.s0 = 0.0f;
        drag.s1 = 0.0f;
        drag.s2 = 0.0f;
        drag.s3 = 0.0f;
        
        //Distance Utils
         float diffX, diffY, dist;
        
        //Initialise Velocity
        initVel.s0 = muzVel * sin(elRad) * cos(az);
        initVel.s1 = muzVel * sin(elRad) * sin(az);
        initVel.s2 = muzVel * cos(elRad);
        
        //Setup the launch Position
        //inPos is used in the loop as currPos
        //Specifying here to convert from float8 to float4
        inPos.s0 = laPos.s0;
        inPos.s1 = laPos.s1;
        inPos.s2 = laPos.s2;
        inVel = initVel;
        
        __private int cnt = 0;
        while (inPos.s2 > -1.0f){
            cnt++;
            //Eval 1
            
            //Euler Velocity
            vel1 = inVel + (inAcc * 0.0f);
            //Euler Position
            pos1 = inPos + (vel1 * 0.0f) + ((inAcc * 0.0f)*0.5f);
            
            //Drag and accels
            combVel = sqrt(pow(vel1.s0, 2)+pow(vel1.s1, 2)+pow(vel1.s2, 2));
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            dragForce = mortSigma*1.225f*pow(combVel, 2);
            //Normalise vector
            normVel = vel1 / combVel;
            //Drag Components
            drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            drag.s2+=((mortMass*9.801f)*-1.0f);
            //Acceleration components
            acc1 = drag/mortMass;

            //Eval 2
            //Euler Velocity
            vel2 = vel1 + (acc1 * (timeStep*0.5f));
            //Euler Position
            pos2 = pos1 + (vel2 * (timeStep*0.5f)) + ((acc1 * pow((timeStep*0.5f), 2))*0.5f);

            //Drag and accels
            combVel = sqrt(pow(vel1.s0, 2)+pow(vel1.s1, 2)+pow(vel1.s2, 2));
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            dragForce = mortSigma*1.225f*pow(combVel, 2.0f);
            //Normalise vector
            normVel = vel2 / combVel;
            //Drag Components
            drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            drag.s2+=((mortMass*9.801f)*-1.0f);
            //Acceleration components
            acc2 = drag/mortMass;

            //Eval 3
            //Euler Velocity
            vel3 = vel2 + (acc2 * (timeStep*0.5f));
            //Euler Position
            pos3 = pos2 + (vel3 * (timeStep*0.5f)) + ((acc2 * pow((timeStep*0.5f), 2))*0.5f);
            
            //Drag and accels
            combVel = sqrt(pow(vel3.s0, 2)+pow(vel3.s1, 2)+pow(vel3.s2, 2));
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            dragForce = mortSigma*1.225f*powr(combVel, 2);
            //Normalise vector
            normVel = vel3 / combVel;
            //Drag Components
            drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            drag.s2+=((mortMass*9.801f)*-1.0f);
            //Acceleration components
            acc3 = drag/mortMass;

            //Eval 4
            //Euler Velocity
            vel4 = vel3 + (acc3 * timeStep);
            //Euler Position
            pos4 = pos3 + (vel4 * timeStep) + ((acc3 * pow(timeStep, 2))*0.5f);
            
            //Drag and accels
            combVel = sqrt(pow(vel4.s0, 2)+pow(vel4.s1, 2)+pow(vel4.s2, 2));
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            dragForce = mortSigma*1.225f*pow(combVel, 2);
            //Normalise vector
            normVel = vel4 / combVel;
            //Drag Components
            drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            drag.s2+=((mortMass*9.801f)*-1.0f);
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
        diffX = fabs(inPos.s0-tgtPos.s0);
        diffY = fabs(inPos.s1-tgtPos.s1);
        dist = sqrt(pow(diffX, 2)+pow(diffY, 2));
        distArr[gid] = dist;
        //printf("%f, %f, %f, %i\n", inPos.s0, inPos.s1, inPos.s2, cnt);
    } else {
        //no run flag found
    }
}




__kernel void shotOptThetaFast(__global float* thetaArr, __global float* distArr, float4 tgtPos, float8 laPos, float d_timeStep)
{
    //Shot Kernel for Optimisation / Brute Forcing the theta shots (parrellelised based on all theta possibilities)
    
    int gid = get_global_id (0);
    
    //Check for run flag - now thinned before
    //if (laPos.s5 == 1.0f){
        
        __private float4 inAcc;
        __private float timeStep = d_timeStep;
        __private float mortSigma;
        __private float mortMass = 3.2f;
        __private float muzVel = 225.1f;
        __private float calibre = 0.081f;
        __private float area;
        __private float dragCoef = 0.15f;
        __private float pi = 3.14159265359f;
        //Calc Mortsigma
        area = pi*(pow((0.5f*calibre), 2));
        //Sigma
        mortSigma = dragCoef*area*0.5f;
        
        //Get Theta from GID
        __private float elRad;
        //Convert to radians
        elRad = ((90.0f-laPos.s3)*0.01745329252f);
        
        //Calc Az
        __private float dX = tgtPos.s0-laPos.s0;
        __private float dY = tgtPos.s1-laPos.s1;
        __private float az = atan2(dY, dX);
        if (az < 0.0f) {
            az+=(2.0f*pi);
        }
        
        //Initialise Velocity
        __private float4 initVel;
        initVel.s0 = muzVel * sin(elRad) * cos(az);
        initVel.s1 = muzVel * sin(elRad) * sin(az);
        initVel.s2 = muzVel * cos(elRad);
        
        //Setup the launch Position
        //inPos is used in the loop as currPos
        //Specifying here to convert from float8 to float4
        __private float4 inPos;
        inPos.s0 = laPos.s0;
        inPos.s1 = laPos.s1;
        inPos.s2 = laPos.s2;
        __private float4 inVel;
        inVel = initVel;
        
        __private int cnt = 0;
        while (inPos.s2 > -1.0f){
            cnt++;
            float4 vel1, vel2, vel3, vel4;
            float4 pos1, pos2, pos3, pos4;
            float4 acc1, acc2, acc3, acc4;
            float4 tayVel;
            float4 tayPos;
            //Forces
            float combVel = 0.0f;
            float dragForce = 0.0f;
            //Normalised Velocity
            float4 normVel;
            //Drag Components
            float4 drag;
            
            //Eval 1
            
            //Euler Velocity
            vel1 = inVel + (inAcc * 0.0f);
            //Euler Position
            pos1 = inPos + (vel1 * 0.0f) + ((inAcc * 0.0f)*0.5f);
            
            //Drag and accels
            //combVel = sqrt(pow(vel1.s0, 2)+pow(vel1.s1, 2)+pow(vel1.s2, 2));
            combVel = half_sqrt(vel1.s0*vel1.s0+vel1.s1*vel1.s1+vel1.s2*vel1.s2);
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            //dragForce = mortSigma*1.225f*pow(combVel, 2);
            dragForce = mortSigma*1.225f*(combVel*combVel);
            //Normalise vector
            normVel = vel1 / combVel;
            //normVel = vel1 * combVel;
            //Drag Components
            //drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            //drag.s2+=((mortMass*9.801f)*-1.0f);
            drag = -normVel * dragForce;
            //Add Gravity force
            drag.s2-=mortMass*9.801f;
            //Acceleration components
            acc1 = drag/mortMass;
            
            //Eval 2
            //Euler Velocity
            vel2 = vel1 + (acc1 * (timeStep*0.5f));
            //Euler Position
            pos2 = pos1 + (vel2 * (timeStep*0.5f)) + ((acc1 * pow((timeStep*0.5f), 2))*0.5f);
            
            //Drag and accels
            //combVel = sqrt(pow(vel2.s0, 2)+pow(vel2.s1, 2)+pow(vel2.s2, 2));
            combVel = half_sqrt(vel2.s0*vel2.s0+vel2.s1*vel2.s1+vel2.s2*vel2.s2);
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            //dragForce = mortSigma*1.225f*pow(combVel, 2.0f);
            dragForce = mortSigma*1.225f*(combVel*combVel);
            //Normalise vector
            normVel = vel2 / combVel;
            //normVel = vel2 * combVel;
            //Drag Components
            //drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            //drag.s2+=((mortMass*9.801f)*-1.0f);
            drag = -normVel * dragForce;
            //Add Gravity force
            drag.s2-=mortMass*9.801f;
            //Acceleration components
            acc2 = drag/mortMass;
            
            //Eval 3
            //Euler Velocity
            vel3 = vel2 + (acc2 * (timeStep*0.5f));
            //Euler Position
            pos3 = pos2 + (vel3 * (timeStep*0.5f)) + ((acc2 * pow((timeStep*0.5f), 2))*0.5f);
            
            //Drag and accels
            //combVel = sqrt(pow(vel3.s0, 2)+pow(vel3.s1, 2)+pow(vel3.s2, 2));
            combVel = half_sqrt(vel3.s0*vel3.s0+vel3.s1*vel3.s1+vel3.s2*vel3.s2);
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            //dragForce = mortSigma*1.225f*powr(combVel, 2);
            dragForce = mortSigma*1.225f*(combVel*combVel);
            //Normalise vector
            normVel = vel3 / combVel;
            //normVel = vel3 * combVel;
            //Drag Components
            //drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            //drag.s2+=((mortMass*9.801f)*-1.0f);
            drag = -normVel * dragForce;
            //Add Gravity force
            drag.s2-=mortMass*9.801f;
            //Acceleration components
            acc3 = drag/mortMass;
            
            //Eval 4
            //Euler Velocity
            vel4 = vel3 + (acc3 * timeStep);
            //Euler Position
            pos4 = pos3 + (vel4 * timeStep) + ((acc3 * pow(timeStep, 2))*0.5f);
            
            //Drag and accels
            //combVel = sqrt(pow(vel4.s0, 2)+pow(vel4.s1, 2)+pow(vel4.s2, 2));
            combVel = half_sqrt(vel4.s0*vel4.s0+vel4.s1*vel4.s1+vel4.s2*vel4.s2);
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            //dragForce = mortSigma*1.225f*pow(combVel, 2);
            dragForce = mortSigma*1.225f*(combVel*combVel);
            //Normalise vector
            normVel = vel4 / combVel;
            //normVel = vel4 * combVel;
            //Drag Components
            //drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            //drag.s2+=((mortMass*9.801f)*-1.0f);
            drag = -normVel * dragForce;
            //Add Gravity force
            drag.s2-=mortMass*9.801f;
            //Acceleration components
            acc4 = drag/mortMass;
            
            //Taylor Expansion
            //tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (1.0f/6.0f);
            //inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (1.0f/6.0f);
            //tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (1.0f/6.0f);
            tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (0.166666f);
            inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (0.166666f);
            tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (0.166666f);
            
            //Swap ready for next iteration
            inPos = inPos + (tayVel * timeStep);
            inVel = inVel + (inAcc * timeStep);
        
            
        }
        //printf("%f, %f, %f\n", inPos.s0, inPos.s1, inPos.s2);
        //Calculate the distance to tgt
        float diffX, diffY, dist;
        diffX = fabs(inPos.s0-tgtPos.s0);
        diffY = fabs(inPos.s1-tgtPos.s1);
        dist = half_sqrt(diffX*diffX+diffY*diffY);
        distArr[gid] = dist;
        //printf("%i, %f\n", cnt, dist);
}



__kernel void shotOptDem(__global float8* demArr, float4 laPos, float d_timeStep)
{
    //Shot Kernel for Optimisation / Brute Forcing the theta shots (parrellelised based on all DEM possibilities)
    
    int gid = get_global_id (0);
    float8 tgtPos = demArr[gid];
    
    //Check for run flag
    if (tgtPos.s5 == 1.0f){
        
        float4 outPos;
        outPos.s0 = 0.0f,
        outPos.s1 = 0.0f,
        outPos.s2 = 0.0f,
        outPos.s3 = 0.0f;
        float4 inPos;
        inPos.s0 = 0.0f,
        inPos.s1 = 0.0f,
        inPos.s2 = 0.0f,
        inPos.s3 = 0.0f;
        float4 initPos;
        initPos.s0 = 0.0f,
        initPos.s1 = 0.0f,
        initPos.s2 = 0.0f,
        initPos.s3 = 0.0f;
        float4 outVel;
        outVel.s0 = 0.0f,
        outVel.s1 = 0.0f,
        outVel.s2 = 0.0f,
        outVel.s3 = 0.0f;
        float4 inVel;
        inVel.s0 = 0.0f,
        inVel.s1 = 0.0f,
        inVel.s2 = 0.0f,
        inVel.s3 = 0.0f;
        float4 initVel;
        initVel.s0 = 0.0f,
        initVel.s1 = 0.0f,
        initVel.s2 = 0.0f,
        initVel.s3 = 0.0f;
        float4 outAcc;
        outAcc.s0 = 0.0f,
        outAcc.s1 = 0.0f,
        outAcc.s2 = 0.0f,
        outAcc.s3 = 0.0f;
        float4 inAcc;
        inAcc.s0 = 0.0f,
        inAcc.s1 = 0.0f,
        inAcc.s2 = 0.0f,
        inAcc.s3 = 0.0f;
        float4 initAcc;
        initAcc.s0 = 0.0f,
        initAcc.s1 = 0.0f,
        initAcc.s2 = 0.0f,
        initAcc.s3 = 0.0f;
        float4 vel1, vel2, vel3, vel4;
        float4 pos1, pos2, pos3, pos4;
        float4 acc1, acc2, acc3, acc4;
        float4 tayVel;
        tayVel.s0 = 0.0f,
        tayVel.s1 = 0.0f,
        tayVel.s2 = 0.0f,
        tayVel.s3 = 0.0f;
        float4 tayPos;
        tayPos.s0 = 0.0f,
        tayPos.s1 = 0.0f,
        tayPos.s2 = 0.0f,
        tayPos.s3 = 0.0f;
        
        float timeStep = d_timeStep;
        float mortSigma;
        float mortMass = 3.2f;
        float muzVel = 225.1f;
        float calibre = 0.081f;
        float area;
        float dragCoef = 0.15f;
        float pi = 3.14159265359f;
        //Calc Mortsigma
        area = pi*(pow((0.5f*calibre), 2));
        //Sigma
        mortSigma = dragCoef*area*0.5f;
        
        //Get Theta from GID
        float elRad;
        //Convert to radians
        elRad = ((90.0f-laPos.s3)*0.01745329252f);
        
        //Calc Az
        float dX = tgtPos.s0-laPos.s0;
        float dY = tgtPos.s1-laPos.s1;
        float az = atan2(dY, dX);
        if (az < 0.0f) {
            az+=(2.0f*pi);
        }
        
        //Forces
        float combVel = 0.0f;
        float dragForce = 0.0f;
        //Normalised Velocity
        float4 normVel;
        normVel.s0 = 0.0f;
        normVel.s1 = 0.0f;
        normVel.s2 = 0.0f;
        normVel.s3 = 0.0f;
        
        //Drag Components
        float4 drag;
        drag.s0 = 0.0f;
        drag.s1 = 0.0f;
        drag.s2 = 0.0f;
        drag.s3 = 0.0f;
        
        //Distance Utils
        float diffX, diffY, dist;
        
        //Initialise Velocity
        initVel.s0 = muzVel * sin(elRad) * cos(az);
        initVel.s1 = muzVel * sin(elRad) * sin(az);
        initVel.s2 = muzVel * cos(elRad);
        
        //Setup the launch Position
        //inPos is used in the loop as currPos
        //Specifying here to convert from float8 to float4
        inPos.s0 = laPos.s0;
        inPos.s1 = laPos.s1;
        inPos.s2 = laPos.s2;
        inVel = initVel;
        
        __private int cnt = 0;
        while (inPos.s2 > -1.0f){
            cnt++;
            //Eval 1
            
            //Euler Velocity
            vel1 = inVel + (inAcc * 0.0f);
            //Euler Position
            pos1 = inPos + (vel1 * 0.0f) + ((inAcc * 0.0f)*0.5f);
            
            //Drag and accels
            combVel = sqrt(pow(vel1.s0, 2)+pow(vel1.s1, 2)+pow(vel1.s2, 2));
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            dragForce = mortSigma*1.225f*pow(combVel, 2);
            //Normalise vector
            normVel = vel1 / combVel;
            //Drag Components
            drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            drag.s2+=((mortMass*9.801f)*-1.0f);
            //Acceleration components
            acc1 = drag/mortMass;
            
            //Eval 2
            //Euler Velocity
            vel2 = vel1 + (acc1 * (timeStep*0.5f));
            //Euler Position
            pos2 = pos1 + (vel2 * (timeStep*0.5f)) + ((acc1 * pow((timeStep*0.5f), 2))*0.5f);
            
            //Drag and accels
            combVel = sqrt(pow(vel1.s0, 2)+pow(vel1.s1, 2)+pow(vel1.s2, 2));
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            dragForce = mortSigma*1.225f*pow(combVel, 2.0f);
            //Normalise vector
            normVel = vel2 / combVel;
            //Drag Components
            drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            drag.s2+=((mortMass*9.801f)*-1.0f);
            //Acceleration components
            acc2 = drag/mortMass;
            
            //Eval 3
            //Euler Velocity
            vel3 = vel2 + (acc2 * (timeStep*0.5f));
            //Euler Position
            pos3 = pos2 + (vel3 * (timeStep*0.5f)) + ((acc2 * pow((timeStep*0.5f), 2))*0.5f);
            
            //Drag and accels
            combVel = sqrt(pow(vel3.s0, 2)+pow(vel3.s1, 2)+pow(vel3.s2, 2));
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            dragForce = mortSigma*1.225f*powr(combVel, 2);
            //Normalise vector
            normVel = vel3 / combVel;
            //Drag Components
            drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            drag.s2+=((mortMass*9.801f)*-1.0f);
            //Acceleration components
            acc3 = drag/mortMass;
            
            //Eval 4
            //Euler Velocity
            vel4 = vel3 + (acc3 * timeStep);
            //Euler Position
            pos4 = pos3 + (vel4 * timeStep) + ((acc3 * pow(timeStep, 2))*0.5f);
            
            //Drag and accels
            combVel = sqrt(pow(vel4.s0, 2)+pow(vel4.s1, 2)+pow(vel4.s2, 2));
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            dragForce = mortSigma*1.225f*pow(combVel, 2);
            //Normalise vector
            normVel = vel4 / combVel;
            //Drag Components
            drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            drag.s2+=((mortMass*9.801f)*-1.0f);
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
        diffX = fabs(inPos.s0-tgtPos.s0);
        diffY = fabs(inPos.s1-tgtPos.s1);
        dist = sqrt(pow(diffX, 2)+pow(diffY, 2));
        //Write to DEMSP
        tgtPos.s4 = dist;
        demArr[gid] = tgtPos;
        //printf("%f, %f, %f, %i\n", inPos.s0, inPos.s1, inPos.s2, cnt);
    }
}



__kernel void shotOptDemFast(__global float8* demArr, float4 laPos, float d_timeStep, __global float4* optimArr)
{
    //Shot Kernel for Optimisation half
    
    
    int gid = get_global_id (0);
    float8 tgtPos = demArr[gid];
    
    //Check for run flag
    if (tgtPos.s5 == 1.0f){
        
        __private float4 inAcc;
        __private float timeStep = d_timeStep;
        __private float mortSigma;
        __private float mortMass = 3.2f;
        __private float muzVel = 225.1f;
        __private float calibre = 0.081f;
        __private float area;
        __private float dragCoef = 0.15f;
        __private float pi = 3.14159265359f;
        //Calc Mortsigma
        area = pi*(pow((0.5f*calibre), 2));
        //Sigma
        mortSigma = dragCoef*area*0.5f;
        
        //Get Theta from GID
        __private float elRad;
        //Convert to radians
        elRad = ((90.0f-laPos.s3)*0.01745329252f);
        
        //Calc Az
        __private float dX = tgtPos.s0-laPos.s0;
        __private float dY = tgtPos.s1-laPos.s1;
        __private float az = atan2(dY, dX);
        if (az < 0.0f) {
            az+=(2.0f*pi);
        }
        
        //Initialise Velocity
        __private float4 initVel;
        initVel.s0 = muzVel * sin(elRad) * cos(az);
        initVel.s1 = muzVel * sin(elRad) * sin(az);
        initVel.s2 = muzVel * cos(elRad);
        
        //Setup the launch Position
        //inPos is used in the loop as currPos
        //Specifying here to convert from float8 to float4
        __private float4 inPos;
        inPos.s0 = laPos.s0;
        inPos.s1 = laPos.s1;
        inPos.s2 = laPos.s2;
        __private float4 inVel;
        inVel = initVel;
        
        __private int cnt = 0;
        while (inPos.s2 > -1.0f){
            cnt++;
            float4 vel1, vel2, vel3, vel4;
            float4 pos1, pos2, pos3, pos4;
            float4 acc1, acc2, acc3, acc4;
            float4 tayVel;
            float4 tayPos;
            //Forces
            float combVel = 0.0f;
            float dragForce = 0.0f;
            //Normalised Velocity
            float4 normVel;
            //Drag Components
            float4 drag;
            
            //Eval 1
            
            //Euler Velocity
            vel1 = inVel + (inAcc * 0.0f);
            //Euler Position
            pos1 = inPos + (vel1 * 0.0f) + ((inAcc * 0.0f)*0.5f);
            
            //Drag and accels
            //combVel = sqrt(pow(vel1.s0, 2)+pow(vel1.s1, 2)+pow(vel1.s2, 2));
            combVel = half_sqrt(vel1.s0*vel1.s0+vel1.s1*vel1.s1+vel1.s2*vel1.s2);
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            //dragForce = mortSigma*1.225f*pow(combVel, 2);
            dragForce = mortSigma*1.225f*(combVel*combVel);
            //Normalise vector
            normVel = vel1 / combVel;
            //normVel = vel1 * combVel;
            //Drag Components
            //drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            //drag.s2+=((mortMass*9.801f)*-1.0f);
            drag = -normVel * dragForce;
            //Add Gravity force
            drag.s2-=mortMass*9.801f;
            //Acceleration components
            acc1 = drag/mortMass;
            
            //Eval 2
            //Euler Velocity
            vel2 = vel1 + (acc1 * (timeStep*0.5f));
            //Euler Position
            pos2 = pos1 + (vel2 * (timeStep*0.5f)) + ((acc1 * pow((timeStep*0.5f), 2))*0.5f);
            
            //Drag and accels
            //combVel = sqrt(pow(vel2.s0, 2)+pow(vel2.s1, 2)+pow(vel2.s2, 2));
            combVel = half_sqrt(vel2.s0*vel2.s0+vel2.s1*vel2.s1+vel2.s2*vel2.s2);
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            //dragForce = mortSigma*1.225f*pow(combVel, 2.0f);
            dragForce = mortSigma*1.225f*(combVel*combVel);
            //Normalise vector
            normVel = vel2 / combVel;
            //normVel = vel2 * combVel;
            //Drag Components
            //drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            //drag.s2+=((mortMass*9.801f)*-1.0f);
            drag = -normVel * dragForce;
            //Add Gravity force
            drag.s2-=mortMass*9.801f;
            //Acceleration components
            acc2 = drag/mortMass;
            
            //Eval 3
            //Euler Velocity
            vel3 = vel2 + (acc2 * (timeStep*0.5f));
            //Euler Position
            pos3 = pos2 + (vel3 * (timeStep*0.5f)) + ((acc2 * pow((timeStep*0.5f), 2))*0.5f);
            
            //Drag and accels
            //combVel = sqrt(pow(vel3.s0, 2)+pow(vel3.s1, 2)+pow(vel3.s2, 2));
            combVel = half_sqrt(vel3.s0*vel3.s0+vel3.s1*vel3.s1+vel3.s2*vel3.s2);
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            //dragForce = mortSigma*1.225f*powr(combVel, 2);
            dragForce = mortSigma*1.225f*(combVel*combVel);
            //Normalise vector
            normVel = vel3 / combVel;
            //normVel = vel3 * combVel;
            //Drag Components
            //drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            //drag.s2+=((mortMass*9.801f)*-1.0f);
            drag = -normVel * dragForce;
            //Add Gravity force
            drag.s2-=mortMass*9.801f;
            //Acceleration components
            acc3 = drag/mortMass;
            
            //Eval 4
            //Euler Velocity
            vel4 = vel3 + (acc3 * timeStep);
            //Euler Position
            pos4 = pos3 + (vel4 * timeStep) + ((acc3 * pow(timeStep, 2))*0.5f);
            
            //Drag and accels
            //combVel = sqrt(pow(vel4.s0, 2)+pow(vel4.s1, 2)+pow(vel4.s2, 2));
            combVel = half_sqrt(vel4.s0*vel4.s0+vel4.s1*vel4.s1+vel4.s2*vel4.s2);
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            //dragForce = mortSigma*1.225f*pow(combVel, 2);
            dragForce = mortSigma*1.225f*(combVel*combVel);
            //Normalise vector
            normVel = vel4 / combVel;
            //normVel = vel4 * combVel;
            //Drag Components
            //drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            //drag.s2+=((mortMass*9.801f)*-1.0f);
            drag = -normVel * dragForce;
            //Add Gravity force
            drag.s2-=mortMass*9.801f;
            //Acceleration components
            acc4 = drag/mortMass;
            
            //Taylor Expansion
            //tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (1.0f/6.0f);
            //inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (1.0f/6.0f);
            //tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (1.0f/6.0f);
            tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (0.166666f);
            inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (0.166666f);
            tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (0.166666f);
            
            //Swap ready for next iteration
            inPos = inPos + (tayVel * timeStep);
            inVel = inVel + (inAcc * timeStep);
            
        }
        //Distance Utils
        float diffX, diffY, dist;
        //Calculate the distance to tgt
        diffX = fabs(inPos.s0-tgtPos.s0);
        diffY = fabs(inPos.s1-tgtPos.s1);
        dist = half_sqrt(diffX*diffX+diffY*diffY);
        //dist = sqrt(pow(diffX, 2)+pow(diffY, 2));
        //Write to DEMSP
        tgtPos.s4 = dist;
        demArr[gid] = tgtPos;
        //Set miss position into Optimisation arr
        __private float4 optData;
        optData = optimArr[gid];
        optData.s0 = inPos.s0;
        optData.s1 = inPos.s1;
        optimArr[gid] = optData;
        //printf("%f, %f, %f, %f\n", inPos.s0, inPos.s1, inPos.s2, dist);
    }
}



__kernel void demIntersectFast(__global float8* demArr, float4 laPos, float d_timeStep,
                               __global float* demZArr, float8 geoTrans, float2 cl_rasterSize)

{
    //Intersect - parallelize based on DEM, check intersect during integrator loop
    //Runs once we've found the optimal theta
    //Breaks at any stage if DEM z > posZ
    
    int gid = get_global_id (0);
    float8 tgtPos = demArr[gid];
    
    //Check for run flag
    if (tgtPos.s5 == 1.0f){
        
        __private float4 inAcc;
        __private float timeStep = d_timeStep;
        __private float mortSigma;
        __private float mortMass = 3.2f;
        __private float muzVel = 225.1f;
        __private float calibre = 0.081f;
        __private float area;
        __private float dragCoef = 0.15f;
        __private float pi = 3.14159265359f;
        //Calc Mortsigma
        area = pi*(pow((0.5f*calibre), 2));
        //Sigma
        mortSigma = dragCoef*area*0.5f;
        
        //Get Theta from GID
        __private float elRad;
        //Convert to radians
        elRad = ((90.0f-tgtPos.s3)*0.01745329252f);
        
        //Calc Az
        __private float dX = tgtPos.s0-laPos.s0;
        __private float dY = tgtPos.s1-laPos.s1;
        __private float az = atan2(dY, dX);
        if (az < 0.0f) {
            az+=(2.0f*pi);
        }
        
        //Initialise Velocity
        __private float4 initVel;
        initVel.s0 = muzVel * sin(elRad) * cos(az);
        initVel.s1 = muzVel * sin(elRad) * sin(az);
        initVel.s2 = muzVel * cos(elRad);
        
        //Setup the launch Position
        //inPos is used in the loop as currPos
        //Specifying here to convert from float8 to float4
        __private float4 inPos;
        inPos.s0 = laPos.s0;
        inPos.s1 = laPos.s1;
        inPos.s2 = laPos.s2;
        __private float4 inVel;
        inVel = initVel;
        
        __private int cnt = 0;
        while (inPos.s2 > -1.0f){
            cnt++;
            float4 vel1, vel2, vel3, vel4;
            float4 pos1, pos2, pos3, pos4;
            float4 acc1, acc2, acc3, acc4;
            float4 tayVel;
            float4 tayPos;
            //Forces
            float combVel = 0.0f;
            float dragForce = 0.0f;
            //Normalised Velocity
            float4 normVel;
            //Drag Components
            float4 drag;
            
            //Eval 1
            
            //Euler Velocity
            vel1 = inVel + (inAcc * 0.0f);
            //Euler Position
            pos1 = inPos + (vel1 * 0.0f) + ((inAcc * 0.0f)*0.5f);
            
            //Drag and accels
            //combVel = sqrt(pow(vel1.s0, 2)+pow(vel1.s1, 2)+pow(vel1.s2, 2));
            combVel = half_sqrt(vel1.s0*vel1.s0+vel1.s1*vel1.s1+vel1.s2*vel1.s2);
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            //dragForce = mortSigma*1.225f*pow(combVel, 2);
            dragForce = mortSigma*1.225f*(combVel*combVel);
            //Normalise vector
            normVel = vel1 / combVel;
            //normVel = vel1 * combVel;
            //Drag Components
            //drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            //drag.s2+=((mortMass*9.801f)*-1.0f);
            drag = -normVel * dragForce;
            //Add Gravity force
            drag.s2-=mortMass*9.801f;
            //Acceleration components
            acc1 = drag/mortMass;
            
            //Eval 2
            //Euler Velocity
            vel2 = vel1 + (acc1 * (timeStep*0.5f));
            //Euler Position
            pos2 = pos1 + (vel2 * (timeStep*0.5f)) + ((acc1 * pow((timeStep*0.5f), 2))*0.5f);
            
            //Drag and accels
            //combVel = sqrt(pow(vel2.s0, 2)+pow(vel2.s1, 2)+pow(vel2.s2, 2));
            combVel = half_sqrt(vel2.s0*vel2.s0+vel2.s1*vel2.s1+vel2.s2*vel2.s2);
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            //dragForce = mortSigma*1.225f*pow(combVel, 2.0f);
            dragForce = mortSigma*1.225f*(combVel*combVel);
            //Normalise vector
            normVel = vel2 / combVel;
            //normVel = vel2 * combVel;
            //Drag Components
            //drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            //drag.s2+=((mortMass*9.801f)*-1.0f);
            drag = -normVel * dragForce;
            //Add Gravity force
            drag.s2-=mortMass*9.801f;
            //Acceleration components
            acc2 = drag/mortMass;
            
            //Eval 3
            //Euler Velocity
            vel3 = vel2 + (acc2 * (timeStep*0.5f));
            //Euler Position
            pos3 = pos2 + (vel3 * (timeStep*0.5f)) + ((acc2 * pow((timeStep*0.5f), 2))*0.5f);
            
            //Drag and accels
            //combVel = sqrt(pow(vel3.s0, 2)+pow(vel3.s1, 2)+pow(vel3.s2, 2));
            combVel = half_sqrt(vel3.s0*vel3.s0+vel3.s1*vel3.s1+vel3.s2*vel3.s2);
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            //dragForce = mortSigma*1.225f*powr(combVel, 2);
            dragForce = mortSigma*1.225f*(combVel*combVel);
            //Normalise vector
            normVel = vel3 / combVel;
            //normVel = vel3 * combVel;
            //Drag Components
            //drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            //drag.s2+=((mortMass*9.801f)*-1.0f);
            drag = -normVel * dragForce;
            //Add Gravity force
            drag.s2-=mortMass*9.801f;
            //Acceleration components
            acc3 = drag/mortMass;
            
            //Eval 4
            //Euler Velocity
            vel4 = vel3 + (acc3 * timeStep);
            //Euler Position
            pos4 = pos3 + (vel4 * timeStep) + ((acc3 * pow(timeStep, 2))*0.5f);
            
            //Drag and accels
            //combVel = sqrt(pow(vel4.s0, 2)+pow(vel4.s1, 2)+pow(vel4.s2, 2));
            combVel = half_sqrt(vel4.s0*vel4.s0+vel4.s1*vel4.s1+vel4.s2*vel4.s2);
            //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
            //dragForce = mortSigma*1.225f*pow(combVel, 2);
            dragForce = mortSigma*1.225f*(combVel*combVel);
            //Normalise vector
            normVel = vel4 / combVel;
            //normVel = vel4 * combVel;
            //Drag Components
            //drag = (normVel * dragForce)*-1.0f;
            //Add Gravity force
            //drag.s2+=((mortMass*9.801f)*-1.0f);
            drag = -normVel * dragForce;
            //Add Gravity force
            drag.s2-=mortMass*9.801f;
            //Acceleration components
            acc4 = drag/mortMass;
            
            //Taylor Expansion
            //tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (1.0f/6.0f);
            //inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (1.0f/6.0f);
            //tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (1.0f/6.0f);
            tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (0.166666f);
            inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (0.166666f);
            tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (0.166666f);
            
            //Swap ready for next iteration
            inPos = inPos + (tayVel * timeStep);
            inVel = inVel + (inAcc * timeStep);
            
            //Check for terrain Intersection
            __private float diffX, diffY, pxX, pxY;
            diffX = inPos.s0-geoTrans.s0;
            diffY = inPos.s1-geoTrans.s3;
            pxX = fabs(diffX/geoTrans.s1);
            pxY = fabs(diffY/geoTrans.s5);
            __private int zIdx;
            zIdx = convert_int((pxY*cl_rasterSize.s0)+pxX);
            //printf("%i\n", zIdx);
            //printf("%f, %f\n", demZArr[zIdx], inPos.s2);
            if (demZArr[zIdx] > inPos.s2) {
                tgtPos.s6 = 1.0f;
            } else {
                tgtPos.s6 = 0.0f;
            }
            
        }
        //Write to DEMSP
        //printf("%f\n", tgtPos.s6);
        demArr[gid] = tgtPos;
        //printf("%f, %f, %f, %i\n", inPos.s0, inPos.s1, inPos.s2, cnt);
        
    } //No Run flag
    
}



__kernel void optimise(__global float8* demArr, __global float4* optimArr, float4 laPos)

    { //Optimise func based on dot product
    }
//      //DEM Arr - contains previous theta and distance from landing pt to tgt
//      //Optim arr contains float4(missPosX, missPosY, seedJump, prevDotP)
//      //Run flag
//        
//        int gid = get_global_id (0);
//        //Get the DEM Pos
//        
//        __private float vec1x, vec1y, vec2x, vec2y, out, dotP, el, seedJump, tmp, tmp2, pxX, pxY;
//        __private float launchX, launchY;
//        __private float2 v1, v2;
//        __private float4 optData;
//        __private float8 demPos;
//        
//        demPos = demArr[gid];
//        if (demPos.s5 == 1.0f){
//            //If the elev is 85.0 its the first shot - alter by seed and continue
//            optData = optimArr[gid];
//            if (demPos.s3 == 85.0f){
//                demPos.s3-=optData.s2;
//                //Do dot P for this position
//                vec1x = laPos.s0 - demPos.s0;
//                vec1y = laPos.s1 - demPos.s1;
//                vec2x = optData.s0 - demPos.s0;
//                vec2y = optData.s1 - demPos.s1;
//                v1 = (vec1x, vec1y);
//                v2 = (vec2x, vec2y);
//                dotP = dot(v1, v2);
//            } else {
//                //Work out DotP between launchPos --> TgtPos and missPos --> tgtPos
//                vec1x = laPos.s0 - demPos.s0;
//                vec1y = laPos.s1 - demPos.s1;
//                vec2x = optData.s0 - demPos.s0;
//                vec2y = optData.s1 - demPos.s1;
//                v1 = (vec1x, vec1y);
//                v2 = (vec2x, vec2y);
//                dotP = dot(v1, v2);
//                //printf("%f, %f, %f, %f, %f, %f, %f\n", vec1x, vec1y, vec2x, vec2y, laPos.s0, optData.s3, dotP);
//                //Check and Move
//                if (dotP < 0.0f) {
//                    if (optData.s3 < 0.0f) {
//                        //Too far and prev was too far (negative dot prod between the vectors)
//                        //Increase elevation, no change in seed
//                        //el = optElArr[gid]+seedJumpArr[gid];
//                        demPos.s3 = demPos.s3+optData.s2;
//                    }
//                    if (optData.s3 > 0.0f) {
//                        //Too far and prev too close:
//                        //alter seedjump and increase el
//                        optData.s2=optData.s2/4.0f;
//                        demPos.s3+=optData.s2;
//                    }
//                }
//                else if (dotP > 0.0f) {
//                    if (optData.s3 > 0.0f) {
//                        //Too close and prev too close
//                        //decrease elevation, no change seed
//                        //el = optElArr[gid]-seedJumpArr[gid];
//                        demPos.s3 = demPos.s3-optData.s2;
//                    }
//                    if (optData.s3 < 0.0f) {
//                        //Too close and prev too far
//                        //Change seed, decrease el
//                        optData.s2=optData.s2/4.0f;
//                        demPos.s3-=optData.s2;
//                    }
//                }
//            }
//            //Set Vals
//            optData.s3 = dotP;
//            optimArr[gid]=optData;
//            demArr[gid] = demPos;
//        }
//    }


__kernel void optEasy(__global float8* demArr, __global float8* optimArr,  float4 laPos,
                      float d_timeStep, float distThreshold, float elevThreshold, float elevMin, float elevMax)

{
    //Shot Kernel for Optimisation without dot product stuff
    //Run three integrations
    //Get closest theta
    //__global float8* optimArr,
    //PRESUMES DATA ALREADY THINNED
    
    
    int gid = get_global_id (0);
    float8 tgtPos = demArr[gid];
    float8 optData = optimArr[gid];
    if (optData.s6 != 0.0f) {
        return;
    };
    //float8 optData = {0.0f, 99999.0f, elevMin, elevMax, elevMin+((elevMax-elevMin)/2.0f), 10.0f, 0.0f, 0.0f};
        
     float4 inAcc;
     float mortSigma, timeStep, mortMass, muzVel, calibre, area, dragCoef, pi, dX, dY, az;
     float diffX, diffY, diffZ, dist;
    timeStep = d_timeStep;
    mortMass = 3.2f;
    muzVel = 192.1f;
    calibre = 0.081f;
    dragCoef = 0.15f;
    pi = 3.14159265359f;
    //Calc Mortsigma
    area = pi*(pow((0.5f*calibre), 2));
    //Sigma
    mortSigma = dragCoef*area*0.5f;
    
    //Calc Az
    dX = tgtPos.s0-laPos.s0;
    dY = tgtPos.s1-laPos.s1;
    az = atan2(dY, dX);
    if (az < 0.0f) {
        az+=(2.0f*pi);
    }
    
    
    //Run three Shots across spread
        for (int i=2; i<5; i++){
            //Setup the launch Position
            //inPos is used in the loop as currPos
            //Specifying here to convert from float8 to float4
             float4 inPos;
            inPos.s0 = laPos.s0;
            inPos.s1 = laPos.s1;
            inPos.s2 = 0.0f;
            tgtPos.s2 = 0.0f;
        
            //Get Theta from elevs
             float elRad;
            //Convert to radians
            elRad = ((90.0f-optData[i])*(pi/180.0f));
            
            //Initialise Velocity
             float4 initVel;
            initVel.s0 = muzVel * sin(elRad) * cos(az);
            initVel.s1 = muzVel * sin(elRad) * sin(az);
            initVel.s2 = muzVel * cos(elRad);
            
             float4 inVel;
            inVel = initVel;
            
             float floorZ;
            floorZ = -1.0f;
              bool floorMove;
            floorMove = 0;
            
                //Run integration for this elev
             int cnt;
            cnt = 0;
            while (inPos.s2 > floorZ){
                cnt++;
                 float4 vel1, vel2, vel3, vel4;
                 float4 pos1, pos2, pos3, pos4;
                 float4 acc1, acc2, acc3, acc4;
                 float4 tayVel;
                 float4 tayPos;
                //Forces
                 float combVel, dragForce;
                dragForce = 0.0f;
                combVel = 0.0f;
                //Normalised Velocity
                 float4 normVel;
                //Drag Components
                 float4 drag;
                
                //Eval 1
                
                //Euler Velocity
                vel1 = inVel + (inAcc * 0.0f);
                //Euler Position
                pos1 = inPos + (vel1 * 0.0f) + ((inAcc * 0.0f)*0.5f);
                
                //Drag and accels
                combVel = sqrt(pow(vel1.s0, 2)+pow(vel1.s1, 2)+pow(vel1.s2, 2));
                //combVel = half_sqrt(vel1.s0*vel1.s0+vel1.s1*vel1.s1+vel1.s2*vel1.s2);
                //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
                dragForce = mortSigma*1.225f*pow(combVel, 2);
                //dragForce = mortSigma*1.225f*(combVel*combVel);
                //Normalise vector
                normVel = vel1 / combVel;
                //normVel = vel1 * combVel;
                //Drag Components
                drag = (normVel * dragForce)*-1.0f;
                //Add Gravity force
                drag.s2+=((mortMass*9.801f)*-1.0f);
                //drag = -normVel * dragForce;
                //Add Gravity force
                //drag.s2-=mortMass*9.801f;
                //Acceleration components
                acc1 = drag/mortMass;
                
                //Eval 2
                //Euler Velocity
                vel2 = vel1 + (acc1 * (timeStep*0.5f));
                //Euler Position
                pos2 = pos1 + (vel2 * (timeStep*0.5f)) + ((acc1 * pow((timeStep*0.5f), 2))*0.5f);
                
                //Drag and accels
                combVel = sqrt(pow(vel2.s0, 2)+pow(vel2.s1, 2)+pow(vel2.s2, 2));
                //combVel = half_sqrt(vel2.s0*vel2.s0+vel2.s1*vel2.s1+vel2.s2*vel2.s2);
                //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
                dragForce = mortSigma*1.225f*pow(combVel, 2.0f);
                //dragForce = mortSigma*1.225f*(combVel*combVel);
                //Normalise vector
                normVel = vel2 / combVel;
                //normVel = vel2 * combVel;
                //Drag Components
                drag = (normVel * dragForce)*-1.0f;
                //Add Gravity force
                drag.s2+=((mortMass*9.801f)*-1.0f);
                //drag = -normVel * dragForce;
                //Add Gravity force
                //drag.s2-=mortMass*9.801f;
                //Acceleration components
                acc2 = drag/mortMass;
                
                //Eval 3
                //Euler Velocity
                vel3 = vel2 + (acc2 * (timeStep*0.5f));
                //Euler Position
                pos3 = pos2 + (vel3 * (timeStep*0.5f)) + ((acc2 * pow((timeStep*0.5f), 2))*0.5f);
                
                //Drag and accels
                combVel = sqrt(pow(vel3.s0, 2)+pow(vel3.s1, 2)+pow(vel3.s2, 2));
                //combVel = half_sqrt(vel3.s0*vel3.s0+vel3.s1*vel3.s1+vel3.s2*vel3.s2);
                //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
                dragForce = mortSigma*1.225f*pow(combVel, 2);
                //dragForce = mortSigma*1.225f*(combVel*combVel);
                //Normalise vector
                normVel = vel3 / combVel;
                //normVel = vel3 * combVel;
                //Drag Components
                drag = (normVel * dragForce)*-1.0f;
                //Add Gravity force
                drag.s2+=((mortMass*9.801f)*-1.0f);
                //drag = -normVel * dragForce;
                //Add Gravity force
                //drag.s2-=mortMass*9.801f;
                //Acceleration components
                acc3 = drag/mortMass;
                
                //Eval 4
                //Euler Velocity
                vel4 = vel3 + (acc3 * timeStep);
                //Euler Position
                pos4 = pos3 + (vel4 * timeStep) + ((acc3 * pow(timeStep, 2))*0.5f);
                
                //Drag and accels
                combVel = sqrt(pow(vel4.s0, 2)+pow(vel4.s1, 2)+pow(vel4.s2, 2));
                //combVel = half_sqrt(vel4.s0*vel4.s0+vel4.s1*vel4.s1+vel4.s2*vel4.s2);
                //motionUtils::drag(netForce, combVel, mortSigma, outPos.s2);
                dragForce = mortSigma*1.225f*pow(combVel, 2);
                //dragForce = mortSigma*1.225f*(combVel*combVel);
                //Normalise vector
                normVel = vel4 / combVel;
                //normVel = vel4 * combVel;
                //Drag Components
                drag = (normVel * dragForce)*-1.0f;
                //Add Gravity force
                drag.s2+=((mortMass*9.801f)*-1.0f);
                //drag = -normVel * dragForce;
                //Add Gravity force
                //drag.s2-=mortMass*9.801f;
                //Acceleration components
                acc4 = drag/mortMass;
                
                //Taylor Expansion
                //tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (1.0f/6.0f);
                //inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (1.0f/6.0f);
                //tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (1.0f/6.0f);
                tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (0.166666f);
                inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (0.166666f);
                tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (0.166666f);
                
                //Swap ready for next iteration
                //Save Z for zen chk
                // float tstZ = inPos.s2;
                inPos = inPos + (tayVel * timeStep);
                inVel = inVel + (inAcc * timeStep);
                
                //Check zenith - currently doing flat shots
//                    if (tstZ > inPos.s2 && floorMove == 0) {
//                        //Test if it passed tgt height
//                        if (tstZ > tgtPos.s2) {
//                            //Reached height - move floor
//                            floorZ = tgtPos.s2;
//                            floorMove = 1;
//                        } else {
//                            //Failed to reach zenith - set outPos as big so distance check for optimisation will fail
//                            inPos.s0 = 0.0f;
//                            inPos.s1 = 0.0f;
//                            inPos.s2 = -9999999999.0f;
//                            inPos.s3 = 0.0f;
//                            break;
//                        }
//                    }
                
            }
            //Get Dist for this elev
            //Calculate the distance to tgt
            diffX = fabs(inPos.s0-tgtPos.s0);
            diffY = fabs(inPos.s1-tgtPos.s1);
            diffZ = fabs(inPos.s2-tgtPos.s2);
            dist = half_sqrt((diffX*diffX)+(diffY*diffY)+(diffZ*diffZ));
            //Store if closer
            if (dist < optData.s1){
                optData.s1 = dist;
                optData.s0 = optData[i];
            }
            //printf("%f, %f, %f, %f, %f, %f, %f\n", dist, diffY, diffX, inPos.s0, tgtPos.s0, inPos.s1, tgtPos.s1);
            //Test distance the shot went
//            float diffX2, diffY2, diffZ2, dist2;
//            diffX2 = fabs(tgtPos.s0-laPos.s0);
//            diffY2 = fabs(tgtPos.s1-laPos.s1);
//            diffZ2 = fabs(tgtPos.s2-laPos.s2);
//            dist2 = half_sqrt((diffX2*diffX2)+(diffY2*diffY2));
//            printf("Distance la to tgt: %f\n", dist2);
            
            
        } //Elevs For
    
        

        //Make optimisation decisions
//        if (optData.s1 <= distThreshold) {
//            //Close enough - set success flag
//            optData.s6 = 1.0f;
//        }
//        if (optData.s0 == elevMin) {
//            if (optData.s5 > elevThreshold) {
//                //Keep trying - special case for minElev
//                optData.s5=optData.s5/2.0f;
//                //min
//                optData.s2=elevMin;
//                //max
//                optData.s3=elevMin+optData.s5;
//                //mid
//                optData.s4=elevMin+(optData.s5/2.0f);
//            } else {
//                //Seems min elev is closest, we havn't got within the correct distance - cant reach the tgt
//                //(0.0 - running, 1.0 - success, 2.0 - failed)
//                optData.s6 = 2.0f;
//                //Set run flag in demarr aswell
//                tgtPos.s5 = 0.0f;
//            }
//        } else if (optData.s0 == elevMax) {
//            if (optData.s5 > elevThreshold) {
//                //Keep trying - special case for maxElev
//                optData.s5=optData.s5/2.0f;
//                //min
//                optData.s2=elevMax-optData.s5;
//                //max
//                optData.s3=elevMax;
//                //mid
//                optData.s4=elevMax-(optData.s5/2.0f);
//            } else {
//                //Seems max elev is closest, we havn't got within the correct distance - cant reach the tgt
//                //(0.0 - running, 1.0 - success, 2.0 - failed)
//                optData.s6 = 2.0f;
//                //Set run flag in demarr aswell
//                tgtPos.s5 = 0.0f;
//            }
//        } else {
//            //Check if we've exceeded elev threshold
//            if (optData.s5< elevThreshold) {
//                //Didnt reach dist and have exceeded elev threshold - break at current
//                //(0.0 - running, 1.0 - success, 2.0 - failed)
//                optData.s6 = 2.0f;
//                //Set run flag in demarr aswell
//                tgtPos.s5 = 0.0f;
//            } else {
//                //Keep trying
//                optData.s5=optData.s5/2.0f;
//                //min
//                optData.s2=optData.s0-(optData.s5/2.0f);
//                //Check Boundary
//                if (optData.s2 < elevMin) {
//                    optData.s2 = elevMin;
//                }
//                //max
//                optData.s3=optData.s0+(optData.s5/2.0f);
//                if (optData.s3 > elevMax) {
//                    optData.s3 = elevMax;
//                }
//                //mid
//                optData.s4=optData.s0;
//            }
//        }
        
//    }// Outside dist threshold while
    
    //TODO: Some other stop conditions
    //Write to DEMSP
    //Theta
    tgtPos.s3 = optData.s0;
    //Dist
    tgtPos.s4 = optData.s1;
    demArr[gid] = tgtPos;
    //Record in OptimArr for posterity
    optimArr[gid] = optData;
    //printf("Done %i\n", gid);
}



//
//
//
//
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
//        diffX = fabs(tgtX-missPosArr[gid].s0);
//        diffY = fabs(tgtY-missPosArr[gid].s1);
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
//            vec2x = missPosArr[gid].s0 - tgtX;
//            vec2y = missPosArr[gid].s1 - tgtY;
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
//    currPos.s0 = prevPos.s0 + ((vel.s0 * timeStep) + (0.5f*acc.s0*(powr(timeStep, 2))));
//    currPos.s1 = prevPos.s1 + ((vel.s1 * timeStep) + (0.5f*acc.s1*(powr(timeStep, 2))));
//    currPos.s2 = prevPos.s2 + ((vel.s2 * timeStep) + (0.5f*acc.s2*(powr(timeStep, 2))));
//    //Check if its passed zenith
//    if (floorMove == false) {
//        if (prevPos.s2 > currPos.s2) {
//            floorZ = ptgtPos.s2;
//            floorMove = true;
//            //Record zenith
//            launchPos.s7 = currPos.s2;
//            //printf("%f, %f, %f, %f, %f, %i, %f, %f\n", launchZ, floorZ, tgtZ, currPosZ, combV, lpCnt, prevPosZ, t);
//        }
//    }
//    lpCnt+=1;
//    //Drag Forces
//    //Combined xyz velocity
//    combV = sqrt(powr(vel.s0, 2)+powr(vel.s1, 2)+powr(vel.s2, 2));
//    //Drag Force xyz
//    dragForce = mortSigma*powr(combV, 2);
//    if (combV > 0.0f){
//        dragForce*=-1.0f;
//    }
//    //This actually gives the combined drag - not a velocity but reusing the class
//    force.s0 = dragForce * sin(elRad) * cos(az);
//    force.s1 = dragForce * sin(elRad) * sin(az);
//    force.s2 = dragForce * cos(elRad);
//    //Gravity Force
//    force.s2-=gForce;
//    //Verlet Integration - TODO REMOVE THIS
//    avgAcc.s0 = 0.5*((force.s0/mortMass)+acc.s0);
//    avgAcc.s1 = 0.5*((force.s1/mortMass)+acc.s1);
//    avgAcc.s2 = 0.5*((force.s2/mortMass)+acc.s2);
//    //Alter Velocity - TODO USE FMA HERE
//    vel.s0+=avgAcc.s0*timeStep;
//    vel.s1+=avgAcc.s1*timeStep;
//    vel.s2+=avgAcc.s2*timeStep;
//    //Reset Vertlet
//    acc.s0 = avgAcc.s0;
//    acc.s1 = avgAcc.s1;
//    acc.s2 = avgAcc.s2;
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
//    //pos.s0, pos.s1, pos.s2
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
//    launchPos.s2+=20.0f;
//    ptgtPos.s2+=20.0f;
//    timeStep = 0.01f;
//    pi = 3.1415926f;
//    //Calc Az
//    diffX = fabs(ptgtPos.s0-launchPos.s0);
//    diffY = fabs(ptgtPos.s1-launchPos.s1);
//    az = atan2(diffY, diffX);
//    if (az < 0.0f) {
//        az+=(2.0f*pi);
//    }
//    //Deg2rad
//    elRad = (90.0f-launchPos.s3) * (pi/180.0f);
//    //Vel components
//    vel.s0 = muzVel * sin(elRad) * cos(az);
//    vel.s1 = muzVel * sin(elRad) * sin(az);
//    vel.s2 = muzVel * cos(elRad);
//    //Acc Components
//    acc.s0 = 0.0f;
//    acc.s1 = 0.0f;
//    acc.s2 = 0.0f;
//    //Setup Floor
//    floorZ = -99999.0f;
//    //Reset Forces
//    force = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
//    dragForce=0.0f;
//    //Setup currpos
//    currPos.s0 = launchPos.s0;
//    currPos.s1 = launchPos.s1;
//    currPos.s2 = launchPos.s2;
//    //Run Shot
//    while (currPos.s2 > floorZ) { //While not below floor
//        //Record Previous Position
//        prevPos = currPos;
//        //Total Time
//        t+=timeStep;
//        //Move Position
//        currPos.s0 = prevPos.s0 + ((vel.s0 * timeStep) + (0.5f*acc.s0*(powr(timeStep, 2))));
//        currPos.s1 = prevPos.s1 + ((vel.s1 * timeStep) + (0.5f*acc.s1*(powr(timeStep, 2))));
//        currPos.s2 = prevPos.s2 + ((vel.s2 * timeStep) + (0.5f*acc.s2*(powr(timeStep, 2))));
//        //Check if its passed zenith
//        if (floorMove == false) {
//            if (prevPos.s2 > currPos.s2) {
//                floorZ = ptgtPos.s2;
//                floorMove = true;
//                //Record zenith
//                launchPos.s7 = currPos.s2;
//                //printf("%f, %f, %f, %f, %f, %i, %f, %f\n", launchZ, floorZ, tgtZ, currPosZ, combV, lpCnt, prevPosZ, t);
//            }
//        }
//        lpCnt+=1;
//        //Drag Forces
//        //Combined xyz velocity
//        combV = sqrt(powr(vel.s0, 2)+powr(vel.s1, 2)+powr(vel.s2, 2));
//        //Drag Force xyz
//        dragForce = mortSigma*powr(combV, 2);
//        if (combV > 0.0f){
//            dragForce*=-1.0f;
//        }
//        //This actually gives the combined drag - not a velocity but reusing the class
//        force.s0 = dragForce * sin(elRad) * cos(az);
//        force.s1 = dragForce * sin(elRad) * sin(az);
//        force.s2 = dragForce * cos(elRad);
//        //Gravity Force
//        force.s2-=gForce;
//        //Verlet Integration - TODO REMOVE THIS
//        avgAcc.s0 = 0.5*((force.s0/mortMass)+acc.s0);
//        avgAcc.s1 = 0.5*((force.s1/mortMass)+acc.s1);
//        avgAcc.s2 = 0.5*((force.s2/mortMass)+acc.s2);
//        //Alter Velocity - TODO USE FMA HERE
//        vel.s0+=avgAcc.s0*timeStep;
//        vel.s1+=avgAcc.s1*timeStep;
//        vel.s2+=avgAcc.s2*timeStep;
//        //Reset Vertlet
//        acc.s0 = avgAcc.s0;
//        acc.s1 = avgAcc.s1;
//        acc.s2 = avgAcc.s2;
//    } //End while not below floor loop
//    //Record Miss point in 2float
//    //Record
//    missPosArr[gid] = (float2)(currPos.s0, currPos.s1);
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

