//Kernels for BlosGPU
//


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
    
    //Calc Az - This needs another flag for the shot direction (in/out) 
    dX = tgtPos.s0-laPos.s0;
    dY = tgtPos.s1-laPos.s1;
    az = atan2(dY, dX);
    //if (az < 0.0f) {
    //    az+=(2.0f*pi);
    //}
    //45deg: 0.7853981633974483
    float absDist;
    diffX = fabs(laPos.s0-tgtPos.s0);
    diffY = fabs(laPos.s1-tgtPos.s1);
    diffZ = fabs(laPos.s2-tgtPos.s2);
    absDist = sqrt((diffX*diffX)+(diffY*diffY));
    
    
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
            tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (1.0f/6.0f);
            inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (1.0f/6.0f);
            tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (1.0f/6.0f);
            //tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (0.166666f);
            //inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (0.166666f);
            //tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (0.166666f);
            
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
        diffX = fabs(inPos.s0-laPos.s0);
        diffY = fabs(inPos.s1-laPos.s1);
        //diffZ = fabs(inPos.s2-tgtPos.s2);
        dist = sqrt((diffX*diffX)+(diffY*diffY));
        //Calc diff between travelled and tgt (absolute)
        dist = fabs(dist - absDist);
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
    
    //TODO: Some other stop conditions
    //Write to DEMSP
    //Theta
    tgtPos.s3 = optData.s0;
    //Dist
    tgtPos.s4 = optData.s1;
    demArr[gid] = tgtPos;
    //Record in OptimArr for posterity
    //TESTING - Record azimuth
    optData.s7 = az;
    optimArr[gid] = optData;
    //printf("Done %i\n", gid);
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
        float elRad;
        //Convert to radians
        elRad = ((90.0f-tgtPos.s3)*0.01745329252f);
        
        //Calc Az
        //float dX = tgtPos.s0-laPos.s0;
        //float dY = tgtPos.s1-laPos.s1;
        //float az = atan2(dY, dX);
        float az;
        az = 0.7853981633974483f;
        
        //Initialise Velocity
        __private float4 initVel;
        initVel.s0 = muzVel * sin(elRad) * cos(az);
        initVel.s1 = muzVel * sin(elRad) * sin(az);
        initVel.s2 = muzVel * cos(elRad);
        
        //Setup the launch Position
        //inPos is used in the loop as currPos
        //Specifying here to convert from float8 to float4
        float4 inPos;
        //inPos.s0 = laPos.s0;
        //inPos.s1 = laPos.s1;
        //inPos.s2 = laPos.s2;
        inPos.s0 = 0.0f;
        inPos.s1 = 0.0f;
        inPos.s2 = 0.0f;
        float4 inVel;
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
            combVel = sqrt(vel1.s0*vel1.s0+vel1.s1*vel1.s1+vel1.s2*vel1.s2);
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
            combVel = sqrt(vel2.s0*vel2.s0+vel2.s1*vel2.s1+vel2.s2*vel2.s2);
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
            combVel = sqrt(vel3.s0*vel3.s0+vel3.s1*vel3.s1+vel3.s2*vel3.s2);
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
            combVel = sqrt(vel4.s0*vel4.s0+vel4.s1*vel4.s1+vel4.s2*vel4.s2);
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
            //convert absolute dist to map coord (as we started from 0.0, 0.0)
            float diffX, diffY, pxX, pxY;
            diffX = (laPos.s0+inPos.s0)-geoTrans.s0;
            diffY = (laPos.s1+inPos.s1)-geoTrans.s3;
            pxX = fabs(diffX/geoTrans.s1);
            pxY = fabs(diffY/geoTrans.s5);
            int zIdx;
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


__kernel void optEasy2d(__global float8* demArr, __global float8* optimArr,  float4 laPos,
                      float d_timeStep, float distThreshold, float elevThreshold, float elevMin, float elevMax)

{
    //Shot Kernel for Optimisation runs 2d
    //All we need to know if the distance to target
    //all just in 2d
    
    int gid = get_global_id (0);
    float8 tgtPos = demArr[gid];
    float8 optData = optimArr[gid];
    if (optData.s6 != 0.0f) {
        return;
    };
    //float8 optData = {0.0f, 99999.0f, elevMin, elevMax, elevMin+((elevMax-elevMin)/2.0f), 10.0f, 0.0f, 0.0f};
    
    float2 inAcc, tgtPos2d;
    float mortSigma, timeStep, mortMass, muzVel, calibre, area, dragCoef, pi, dX, dY, az;
    float diffX, diffY, diffZ, dist;
    timeStep = d_timeStep;
    mortMass = 4.2f;
    muzVel = 192.1f;
    calibre = 0.081f;
    dragCoef = 0.15f;
    pi = 3.14159265359f;
    //Calc Mortsigma
    area = pi*(pow((0.5f*calibre), 2));
    //Sigma
    mortSigma = dragCoef*area*0.5f;
    //Absolute Target Pos
    diffX = fabs(laPos.s0-tgtPos.s0);
    diffY = fabs(laPos.s1-tgtPos.s1);
    tgtPos2d.s0 = sqrt((diffX*diffX)+(diffY*diffY));
    
    
    //Run three Shots across spread
    for (int i=2; i<5; i++){
        //Setup the launch Position
        //inPos is used in the loop as currPos
        //Specifying here to convert from float8 to float4
        float2 inPos;
        inPos.s0 = 0.0f;
        inPos.s1 = 0.0f;
        tgtPos2d.s1 = 0.0f;
        
        //Get Theta from elevs
        float elRad;
        //Convert to radians
        elRad = ((90.0f-optData[i])*(pi/180.0f));
        
        //Initialise Velocity
        float2 initVel;
        //Y = SOH - sin(theta)*hyp
        //X = CAH - cos(theta)*hyp
        initVel.s1 = sin(optData[i]) * muzVel;
        initVel.s0 = cos(optData[i]) * muzVel;
        
        float2 inVel;
        inVel = initVel;
        
        float floorZ;
        floorZ = -1.0f;
        bool floorMove;
        floorMove = 0;
        
        //Run integration for this elev
        int cnt;
        cnt = 0;
        while (inPos.s1 > floorZ){
            cnt++;
            float2 vel1, vel2, vel3, vel4;
            float2 pos1, pos2, pos3, pos4;
            float2 acc1, acc2, acc3, acc4;
            float2 tayVel;
            float2 tayPos;
            //Forces
            float combVel, dragForce;
            dragForce = 0.0f;
            combVel = 0.0f;
            //Normalised Velocity
            float2 normVel;
            //Drag Components
            float2 drag;
            
            //Eval 1
            
            //Euler Velocity
            vel1 = inVel + (inAcc * 0.0f);
            //Euler Position
            pos1 = inPos + (vel1 * 0.0f) + ((inAcc * 0.0f)*0.5f);
            
            //Drag and accels
            combVel = sqrt(pow(vel1.s0, 2)+pow(vel1.s1, 2));
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
            drag.s1+=((mortMass*9.801f)*-1.0f);
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
            combVel = sqrt(pow(vel2.s0, 2)+pow(vel2.s1, 2));
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
            drag.s1+=((mortMass*9.801f)*-1.0f);
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
            combVel = sqrt(pow(vel3.s0, 2)+pow(vel3.s1, 2));
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
            drag.s1+=((mortMass*9.801f)*-1.0f);
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
            combVel = sqrt(pow(vel4.s0, 2)+pow(vel4.s1, 2));
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
            drag.s1+=((mortMass*9.801f)*-1.0f);
            //drag = -normVel * dragForce;
            //Add Gravity force
            //drag.s2-=mortMass*9.801f;
            //Acceleration components
            acc4 = drag/mortMass;
            
            //Taylor Expansion
            tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (1.0f/6.0f);
            inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (1.0f/6.0f);
            tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (1.0f/6.0f);
            //tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (0.166666f);
            //inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (0.166666f);
            //tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (0.166666f);
            
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
        dist = fabs(inPos.s0-tgtPos2d.s0);
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
    
    //TODO: Some other stop conditions
    //Write to DEMSP
    //Theta
    tgtPos.s3 = optData.s0;
    //Dist
    tgtPos.s4 = optData.s1;
    demArr[gid] = tgtPos;
    //Record in OptimArr for posterity
    //TESTING - Record azimuth
    optData.s7 = tgtPos2d.s0;
    optimArr[gid] = optData;
    //printf("Done %i\n", gid);
}



__kernel void optEasy3dSingAz(__global float8* demArr, __global float8* optimArr,  float4 laPos,
                      float d_timeStep, float distThreshold, float elevThreshold, float elevMin, float elevMax)

{
    //Shot Kernel for Optimisation
    //3D, but without azimuth calculation - just uses absolute distance to tgt and starts at 0
    //This works much better than with azimuth calc
    //It seems there is something wrong with the azimuth calc because results are not equal at different bearings
    
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
    //dX = tgtPos.s0-laPos.s0;
    //dY = tgtPos.s1-laPos.s1;
    //az = atan2(dY, dX);
    //if (az < 0.0f) {
    //    az+=(2.0f*pi);
    //}
    az = 0.7853981633974483f;
    //Calc absolute distance la - tgt
    float absDist;
    diffX = fabs(laPos.s0-tgtPos.s0);
    diffY = fabs(laPos.s1-tgtPos.s1);
    diffZ = fabs(laPos.s2-tgtPos.s2);
    absDist = sqrt((diffX*diffX)+(diffY*diffY));
    
    
    //Run three Shots across spread
    for (int i=2; i<5; i++){
        //Setup the launch Position
        //inPos is used in the loop as currPos
        //Specifying here to convert from float8 to float4
        float4 inPos;
        inPos.s0 = 0.0f;
        inPos.s1 = 0.0f;
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
            tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (1.0f/6.0f);
            inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (1.0f/6.0f);
            tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (1.0f/6.0f);
            //tayVel = (vel1+((vel2+vel3)*2.0f)+vel4) * (0.166666f);
            //inAcc = (acc1+((acc2+acc3)*2.0f)+acc4) * (0.166666f);
            //tayPos = (pos1+((pos2+pos3)*2.0f)+pos4) * (0.166666f);
            
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
        //Calculate the distance travelled - we started at 0.0
        dist = sqrt((inPos.s0*inPos.s0)+(inPos.s1*inPos.s1));
        //Calc diff between travelled and tgt (absolute)
        dist = fabs(dist - absDist);
        //Store if closer
        if (dist < optData.s1){
            optData.s1 = dist;
            optData.s0 = optData[i];
        }
        //printf("%i\n", cnt);
        //Test distance the shot went
        //            float diffX2, diffY2, diffZ2, dist2;
        //            diffX2 = fabs(tgtPos.s0-laPos.s0);
        //            diffY2 = fabs(tgtPos.s1-laPos.s1);
        //            diffZ2 = fabs(tgtPos.s2-laPos.s2);
        //            dist2 = half_sqrt((diffX2*diffX2)+(diffY2*diffY2));
        //            printf("Distance la to tgt: %f\n", dist2);
        
        
    } //Elevs For
    
    //TODO: Some other stop conditions
    //Write to DEMSP
    //Theta
    tgtPos.s3 = optData.s0;
    //Dist
    tgtPos.s4 = optData.s1;
    demArr[gid] = tgtPos;
    //Record in OptimArr for posterity
    //TESTING - Record azimuth
    //optData.s7 = az;
    optimArr[gid] = optData;
    //printf("Done %i\n", gid);
}


