//
//  simplekern.c
//  gputest1
//
//  Created by Chris Nicholas on 05/06/2015.
//  Copyright (c) 2015 Chris Nicholas. All rights reserved.
//

__kernel void rangeChk(__global float8* worldArr, float tgtX, float tgtY,
                       float minDist, float maxDist)
    {
        //Reduces by Range
        //Input - demArr, min, max ranges
        //Process - distances
        //Output - array of index numbers in demArr that are within range (procIdx)
        int gid = get_global_id (0);

        //Work out the geo coords for this location
        __private float distance, diffX, diffY;
        __private float8 launchPos;
        launchPos = worldArr[gid];
        //Check Ranges
        diffX = fabs(tgtX-launchPos.x);
        diffY = fabs(tgtY-launchPos.y);
        launchPos.s6 = sqrt(powr(diffX, 2)+powr(diffY, 2));
        worldArr[gid] = launchPos;
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


__kernel void shotStep(__global float8* trjArr, float4 totTime, float4 timestep, float4 tgtPos, float mortMass,
                   float muzVel, float mortSigma, float gForce,
                   __global float2* missPosArr)

{ //Single step of the recshot func
    int gid = get_global_id (0);
    //Get curr position from world Arr
    __private float8 currPos;
    __private float8 prevPos;
    //OutPos is the new position
    __private float8 outPos;
    currPos = worldArr[gid];
    
    //Record Previous Position
    prevPos = currPos;
    //Total Time
    totTime+=timeStep;
    //Move Position
    currPos.x = prevPos.x + ((vel.x * timeStep) + (0.5f*acc.x*(powr(timeStep, 2))));
    currPos.y = prevPos.y + ((vel.y * timeStep) + (0.5f*acc.y*(powr(timeStep, 2))));
    currPos.z = prevPos.z + ((vel.z * timeStep) + (0.5f*acc.z*(powr(timeStep, 2))));
    //Check if its passed zenith
    if (floorMove == false) {
        if (prevPos.z > currPos.z) {
            floorZ = ptgtPos.z;
            floorMove = true;
            //Record zenith
            launchPos.s7 = currPos.z;
            //printf("%f, %f, %f, %f, %f, %i, %f, %f\n", launchZ, floorZ, tgtZ, currPosZ, combV, lpCnt, prevPosZ, t);
        }
    }
    lpCnt+=1;
    //Drag Forces
    //Combined xyz velocity
    combV = sqrt(powr(vel.x, 2)+powr(vel.y, 2)+powr(vel.z, 2));
    //Drag Force xyz
    dragForce = mortSigma*powr(combV, 2);
    if (combV > 0.0f){
        dragForce*=-1.0f;
    }
    //This actually gives the combined drag - not a velocity but reusing the class
    force.x = dragForce * sin(elRad) * cos(az);
    force.y = dragForce * sin(elRad) * sin(az);
    force.z = dragForce * cos(elRad);
    //Gravity Force
    force.z-=gForce;
    //Verlet Integration - TODO REMOVE THIS
    avgAcc.x = 0.5*((force.x/mortMass)+acc.x);
    avgAcc.y = 0.5*((force.y/mortMass)+acc.y);
    avgAcc.z = 0.5*((force.z/mortMass)+acc.z);
    //Alter Velocity - TODO USE FMA HERE
    vel.x+=avgAcc.x*timeStep;
    vel.y+=avgAcc.y*timeStep;
    vel.z+=avgAcc.z*timeStep;
    //Reset Vertlet
    acc.x = avgAcc.x;
    acc.y = avgAcc.y;
    acc.z = avgAcc.z;

    
}

__kernel void shot(__global float8* worldArr, float4 tgtPos, float mortMass,
                        float muzVel, float mortSigma, float gForce,
                        __global float2* missPosArr)

{//Runs single shot - returns last position
    int gid = get_global_id (0);
    __private float diffX, diffY;
    __private float az;
    __private float radian;
    __private float4 vel;
    __private float combV;
    __private float4 acc;
    __private float4 avgAcc;
    __private float4 force;
    __private float dragForce;
    __private float4 currPos;
    __private float4 prevPos;
    __private float floorZ;
    __private float elRad;
    //__local float trjArrX[3072];
    //__local float trjArrY[3072];
    //__local float trjArrZ[3072];
    __private float pi;
    __private float timeStep;
    __private float t;
    __private bool floorMove;
    __private bool zenFail;
    __private int lpCnt;
    __private int resVal;
    //Landing Point
    __private float2 missPos;
    missPos = (float2)(0.0f, 0.0f);
    lpCnt=0;
    floorMove = false;
    t = 0.0f;
    zenFail = false;
    //Get float8 from input
    //float8(x, y ,z ,optEl, dotp, prevdot, dist, SPARE)
    //pos.x, pos.y, pos.z
    //pos.s3 - optel
    //pos.s4 - dotp
    //pos.s5 - prevDot
    //pos.s6 - dist land 2 tgt
    __private float8 launchPos;
    launchPos = worldArr[gid];
    __private float4 ptgtPos;
    ptgtPos = tgtPos;

    //Adding obs height to launch Z and tgtZ
    launchPos.z+=20.0f;
    ptgtPos.z+=20.0f;
    timeStep = 0.01f;
    pi = 3.1415926f;
    //Calc Az
    diffX = fabs(ptgtPos.x-launchPos.x);
    diffY = fabs(ptgtPos.y-launchPos.y);
    az = atan2(diffY, diffX);
    if (az < 0.0f) {
        az+=(2.0f*pi);
    }
    //Deg2rad
    elRad = (90.0f-launchPos.s3) * (pi/180.0f);
    //Vel components
    vel.x = muzVel * sin(elRad) * cos(az);
    vel.y = muzVel * sin(elRad) * sin(az);
    vel.z = muzVel * cos(elRad);
    //Acc Components
    acc.x = 0.0f;
    acc.y = 0.0f;
    acc.z = 0.0f;
    //Setup Floor
    floorZ = -99999.0f;
    //Reset Forces
    force = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    dragForce=0.0f;
    //Setup currpos
    currPos.x = launchPos.x;
    currPos.y = launchPos.y;
    currPos.z = launchPos.z;
    //Run Shot
    while (currPos.z > floorZ) { //While not below floor
        //Record Previous Position
        prevPos = currPos;
        //Total Time
        t+=timeStep;
        //Move Position
        currPos.x = prevPos.x + ((vel.x * timeStep) + (0.5f*acc.x*(powr(timeStep, 2))));
        currPos.y = prevPos.y + ((vel.y * timeStep) + (0.5f*acc.y*(powr(timeStep, 2))));
        currPos.z = prevPos.z + ((vel.z * timeStep) + (0.5f*acc.z*(powr(timeStep, 2))));
        //Check if its passed zenith
        if (floorMove == false) {
            if (prevPos.z > currPos.z) {
                floorZ = ptgtPos.z;
                floorMove = true;
                //Record zenith
                launchPos.s7 = currPos.z;
                //printf("%f, %f, %f, %f, %f, %i, %f, %f\n", launchZ, floorZ, tgtZ, currPosZ, combV, lpCnt, prevPosZ, t);
            }
        }
        lpCnt+=1;
        //Drag Forces
        //Combined xyz velocity
        combV = sqrt(powr(vel.x, 2)+powr(vel.y, 2)+powr(vel.z, 2));
        //Drag Force xyz
        dragForce = mortSigma*powr(combV, 2);
        if (combV > 0.0f){
            dragForce*=-1.0f;
        }
        //This actually gives the combined drag - not a velocity but reusing the class
        force.x = dragForce * sin(elRad) * cos(az);
        force.y = dragForce * sin(elRad) * sin(az);
        force.z = dragForce * cos(elRad);
        //Gravity Force
        force.z-=gForce;
        //Verlet Integration - TODO REMOVE THIS
        avgAcc.x = 0.5*((force.x/mortMass)+acc.x);
        avgAcc.y = 0.5*((force.y/mortMass)+acc.y);
        avgAcc.z = 0.5*((force.z/mortMass)+acc.z);
        //Alter Velocity - TODO USE FMA HERE
        vel.x+=avgAcc.x*timeStep;
        vel.y+=avgAcc.y*timeStep;
        vel.z+=avgAcc.z*timeStep;
        //Reset Vertlet
        acc.x = avgAcc.x;
        acc.y = avgAcc.y;
        acc.z = avgAcc.z;
    } //End while not below floor loop
    //Record Miss point in 2float
    //Record
    missPosArr[gid] = (float2)(currPos.x, currPos.y);
}


__kernel void intersect(__global float* worldArr,  __global int* result)
//Run trajectories & profiles
//Input - demArr, procIdx, optEls, weapon Params
//Process - intersection between lines
//Output - int results for blos and tlos
{
}

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

