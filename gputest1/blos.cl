





__kernel void rangeChk(__global float4* worldArr, float tgtX, float tgtY,
                       float minDist, float maxDist)
{
    //Reduces by Range
    //Input - demArr, min, max ranges
    //Process - distances
    //Output - array of index numbers in demArr that are within range (procIdx)
    int gid = get_global_id (0);
    
    //Work out the geo coords for this location
    __private float diffX, diffY;
    __private float4 launchPos;
    launchPos = worldArr[gid];
    //Check Ranges
    diffX = fabs(tgtX-launchPos.x);
    diffY = fabs(tgtY-launchPos.y);
    launchPos.w = sqrt(powr(diffX, 2)+powr(diffY, 2));
    worldArr[gid] = launchPos;
}


__kernel void init(__global float4* demSPArr, __global float4* velArr, __global float4* accArr, __global float* azRadArr, __global float* elRadArr,
                   __global float4* forceArr, float4 tgtPos, float elDeg, float muzzVel, __global float4* posArr)
{
    int gid = get_global_id (0);
    __private float diffX;
    __private float diffY;
    __private float4 vel;
    __private float pi;
    __private float4 launchPos;
    __private float azRad;
    __private float elRad;
    launchPos = posArr[gid];
    
    pi = 3.1415926f;
    //Calc Az
    diffX = fabs(tgtPos.x-launchPos.x);
    diffY = fabs(tgtPos.y-launchPos.y);
    azRad = atan2(diffY, diffX);
    if (azRad < 0.0f) {
        azRad+=(2.0f*pi);
    }
    //Deg2rad
    elRad = (90.0f-elDeg) * (pi/180.0f);
    //Vel components
    vel.x = muzzVel * sin(elRad) * cos(azRad);
    vel.y = muzzVel * sin(elRad) * sin(azRad);
    vel.z = muzzVel * cos(elRad);
    
    //Dump out
    forceArr[gid] = 0.0f;
    accArr[gid] = 0.0f;
    posArr[gid] = 0.0f;
    velArr[gid]=vel;
    azRadArr[gid]=azRad;
    elRadArr[gid]=elRad;
    
}

__kernel void move(__global float4* posArr, __global float4* velArr, __global float4* accArr,
                   float timeStep)
{
    
    int gid = get_global_id (0);
    __private float4 currPos;
    __private float4 newPos;
    __private float4 vel;
    __private float4 acc;
    //Get the current position
    currPos = posArr[gid];
    vel = velArr[gid];
    acc = accArr[gid];
    
    //Move Position
    newPos.x = currPos.x + ((vel.x * timeStep) + (0.5f*acc.x*(powr(timeStep, 2))));
    newPos.y = currPos.y + ((vel.y * timeStep) + (0.5f*acc.y*(powr(timeStep, 2))));
    newPos.z = currPos.z + ((vel.z * timeStep) + (0.5f*acc.z*(powr(timeStep, 2))));
    
    //Reset position
    posArr[gid] = newPos;
    
}



__kernel void drag(__global float4* velArr, __global float4* forceArr, __global float* azRadArr, __global float* elRadArr,
                   float gForce, float mortSigma)
{
    int gid = get_global_id (0);
    __private float combV;
    __private float dragForce;
    __private float4 vel;
    __private float4 force;
    __private float azRad;
    __private float elRad;
    __private float tst;
    azRad = azRadArr[gid];
    elRad = elRadArr[gid];
    
    vel = velArr[gid];
    //Drag Forces
    //Test - delete
    //tst = sqrt(powr(0.04, 2)+powr(1.2, 2)+powr(1.3, 2));
    //printf("%f \n", tst);
    
    //Combined xyz velocity
    if (vel.x < 0.0f){
        vel.x*=-1.0;
    }
    if (vel.y < 0.0f){
        vel.y*=-1.0;
    }
    if (vel.z < 0.0f){
        vel.z*=-1.0;
    }
    combV = sqrt(powr(vel.x, 2)+powr(vel.y, 2)+powr(vel.z, 2));
    //Drag Force xyz
    dragForce = mortSigma*powr(combV, 2);
    if (combV > 0.0f){
        dragForce*=-1.0f;
    }
    //Combined drag
    force.x = dragForce * sin(elRad) * cos(azRad);
    force.y = dragForce * sin(elRad) * sin(azRad);
    force.z = dragForce * cos(elRad);
    //Gravity Force
    //printf("%f, %f, %f, %f \n", vel.z, force.z, combV, gForce);
    force.z-=gForce;
    //printf("%f, %f, %f \n", force.z, combV, gForce);
    
    //Write to force Arr
    forceArr[gid] = force;

}


__kernel void acc(__global float4* accArr, __global float4* forceArr, __global float4* velArr,
                  float mortMass, float timeStep)
{
    int gid = get_global_id (0);
    __private float4 acc;
    __private float4 outAcc;
    __private float4 force;
    __private float4 vel;
    
    acc = accArr[gid];
    vel = velArr[gid];
    force = forceArr[gid];
    
    //Recalc Acceleration from drag forces
    acc.x = force.x/mortMass;
    acc.y = force.y/mortMass;
    acc.z = force.z/mortMass;
    
    //Alter Velocity - TODO USE FMA HERE
    vel.x+=acc.x*timeStep;
    vel.y+=acc.y*timeStep;
    vel.z+=acc.z*timeStep;
    
//    //Verlet Integration - TODO REMOVE THIS
//    outAcc.x = 0.5*((force.x/mortMass)+acc.x);
//    outAcc.y = 0.5*((force.y/mortMass)+acc.y);
//    outAcc.z = 0.5*((force.z/mortMass)+acc.z);
//    //Alter Velocity - TODO USE FMA HERE
//    vel.x+=outAcc.x*timeStep;
//    vel.y+=outAcc.y*timeStep;
//    vel.z+=outAcc.z*timeStep;
    
    //Dump
    accArr[gid] = acc;
    velArr[gid] = vel;

}






