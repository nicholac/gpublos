//
//  simplekern.c
//  gputest1
//
//  Created by Chris Nicholas on 05/06/2015.
//  Copyright (c) 2015 Chris Nicholas. All rights reserved.
//

__kernel void rangeChk(__global float4* worldArr, float tgtX, float tgtY,
                       float minDist, float maxDist, __global float* distArr)
    {
        //Reduces by Range
        //Input - demArr, min, max ranges
        //Process - distances
        //Output - array of index numbers in demArr that are within range (procIdx)
        int gid = get_global_id (0);

        //Work out the geo coords for this location
        __private float distance, diffX, diffY;
        __private float4 launchPos;
        launchPos = worldArr[gid];
        //Check Ranges
        diffX = fabs(tgtX-launchPos.x);
        diffY = fabs(tgtY-launchPos.y);
        distArr[gid] = sqrt(powr(diffX, 2)+powr(diffY, 2));
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

//__kernel void shot(__global float* dem, float tgtX, float tgtY, float rasterXSize,
//                        float rasterYSize, float gtrans0,
//                        float gtrans1, float gtrans2, float gtrans3, float gtrans4,
//                        float gtrans5, float mortMass,
//                        float muzVel, float mortSigma, float gForce, float tgtZ,
//                        __global float2* coordArr)
//
//{//Runs single shot - returns last position
//    int gid = get_global_id (0);
//    __private float tmp;
//    __private float tmp2;
//    __private float tmp3;
//    __private float pxX;
//    __private float pxY;
//    __private float dist;
//    dist = 0.0;
//    __private float diffX;
//    __private float diffY;
//    __private float az;
//    __private float radian;
//    __private float launchX;
//    __private float launchY;
//    __private float launchZ;
//    __private float velX;
//    __private float velY;
//    __private float velZ;
//    __private float combV;
//    __private float accX;
//    __private float accY;
//    __private float accZ;
//    __private float avgAccX;
//    __private float avgAccY;
//    __private float avgAccZ;
//    __private float forceX;
//    __private float forceY;
//    __private float forceZ;
//    __private float dragForce;
//    __private float currPosX;
//    __private float currPosY;
//    __private float currPosZ;
//    __private float prevPosX;
//    __private float prevPosY;
//    __private float prevPosZ;
//    __private float floorZ;
//    __private float elDeg;
//    __private float elRad;
//    //__local float trjArrX[3072];
//    //__local float trjArrY[3072];
//    //__local float trjArrZ[3072];
//    //OPTIMISE
//    //LaunchAngle
//    elDeg = 55.0f;
//    __private float pi;
//    __private float timeStep;
//    __private float t;
//    __private bool floorMove;
//    __private bool zenFail;
//    __private int lpCnt;
//    __private int resVal;
//    //Landing Point
//    __local float2 missPos;
//    missPos=(0.0f, 0.0f);
//    lpCnt=0;
//    floorMove = false;
//    t = 0.0f;
//    zenFail = false;
//    //i is still the position in the DEM (1d) as thats the input task size for kern
//
//    //Check range check results from procIdx
//    if (procIdx[gid] == 1){
//        //Run Optimise
//        //Adding obs height to launch Z and tgtZ
//        launchZ = dem[gid]+20.0f;
//        //printf("%f\n", launchZ);
//        tgtZ+=20.0f;
//        timeStep = 0.01f;
//        pi = 3.1415926f;
//        //Work out pixel location and launch from 1d index
//        tmp = modf(float(gid)/rasterXSize, &pxY);
//        tmp2 = modf(tmp * rasterXSize, &pxX);
//        launchX = gtrans0+(pxX*gtrans1)+pxY*gtrans2;
//        launchY = gtrans3+(pxX*gtrans4)+pxY*gtrans5;
//        //RECSHOT
//        //Calc Az
//        diffX = fabs(tgtX-launchX);
//        diffY = fabs(tgtY-launchY);
//        az = atan2(diffY, diffX);
//        if (az < 0.0f) {
//            az+=(2.0f*pi);
//        }
//        //Deg2rad
//        elRad = (90.0f-elDeg) * (pi/180.0f);
//        //Vel components
//        velX = muzVel * sin(elRad) * cos(az);
//        velY = muzVel * sin(elRad) * sin(az);
//        velZ = muzVel * cos(elRad);
//        //Acc Components
//        accX = 0.0f;
//        accY = 0.0f;
//        accZ = 0.0f;
//        //Push first elevation
//        //trjArrX[0]=launchX;
//        //trjArrY[0]=launchY;
//        //trjArrZ[0]=dem[gid];
//        //Setup Floor
//        floorZ = -99999.0f;
//        //Reset Force
//        //Reset time
//        forceZ=0.0f;
//        dragForce=0.0f;
//        //Setup currpos
//        currPosX = launchX;
//        currPosY = launchY;
//        currPosZ = launchZ;
//        currPosZ = 0.0f;
//        //printf("%i, %i, %i, %i\n", 1, optLps, resVal, lpCnt);
//        //Next - do rest of recshot (once)
//        while (currPosZ > floorZ) { //While not below floor
//            //Record Previous Position
//            prevPosX = currPosX;
//            prevPosY = currPosY;
//            prevPosZ = currPosZ;
//            //Total Time
//            t = t+timeStep;
//            //Move Position
//            currPosX = prevPosX + ((velX * timeStep) + (0.5f*accX*(powr(timeStep, 2))));
//            currPosY = prevPosY + ((velY * timeStep) + (0.5f*accY*(powr(timeStep, 2))));
//            currPosZ = prevPosZ + ((velZ * timeStep) + (0.5f*accZ*(powr(timeStep, 2))));
//            //Check if its passed zenith
//            if (floorMove == false) {
//                if (prevPosZ > currPosZ) {
//                    //Passed Zenith - Check if we reached tgt height
//                    //if (prevPosZ < tgtZ){
//                        //Didnt reach the tgt height
//                        //zenFail = true;
//                        //printf("%f, %f, %f, %f, %f, %f, %i, %f, %f\n", launchZ, floorZ, tgtZ, prevPosZ, currPosZ, combV, lpCnt, prevPosZ, t);
//                        //break;
//                    //}
//                    floorZ = tgtZ;
//                    floorMove = true;
//                    //printf("%f, %f, %f, %f, %f, %i, %f, %f\n", launchZ, floorZ, tgtZ, currPosZ, combV, lpCnt, prevPosZ, t);
//                }
//                //Check if the zenith failed to reach tgt height
//                //if (zenFail == true){
//                //    //Return null result
//                //    resVal = 1;
//                //}
//            }
//            //TODO: Check if its outside the bounds of DEM
//            lpCnt+=1;
//            //Drag Forces
//            //Combined xyz velocity
//            combV = sqrt(powr(velX, 2)+powr(velY, 2)+powr(velZ, 2));
//            //Drag Force xyz
//            dragForce = mortSigma*powr(combV, 2);
//            if (combV > 0.0f){
//                dragForce*=-1.0f;
//            }
//            //This actually gives the combined drag - not a velocity but reusing the class
//            forceX = dragForce * sin(elRad) * cos(az);
//            forceY = dragForce * sin(elRad) * sin(az);
//            forceZ = dragForce * cos(elRad);
//            //printf("%f, %f, %f, %f\n", velX, velY, velZ, combV);
//            //Gravity Force
//            forceZ-=gForce;
//            //Verlet Integration
//            avgAccX = 0.5f*((forceX/mortMass)+accX);
//            avgAccY = 0.5f*((forceY/mortMass)+accY);
//            avgAccZ = 0.5f*((forceZ/mortMass)+accZ);
//            //Alter Velocity
//            velX+=avgAccX*timeStep;
//            velY+=avgAccY*timeStep;
//            velZ+=avgAccZ*timeStep;
//            //Reset Verlet Accelerations
//            accX = avgAccX;
//            accY = avgAccY;
//            accZ = avgAccZ;
//
//        } //End while not below floor loop
//        //Record Miss point in 2float
//        missPos.xy = (float2)(currPosX, currPosY);
//        //Check Distance from last point to tgt - 2d
////        diffX = fabs(tgtX-currPosX);
////        diffY = fabs(tgtY-currPosY);
////        dist = distance(diffX, diffY);
//        //printf("%f\n", dist);
//    }//If proxIdx end
//    coordArr[gid] = missPos;
//}


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

