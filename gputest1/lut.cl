//Massive LUT implemetation


__kernel void genDists(__global float8* tgtArr, float4 chkPt)
{
    //Generates distance 2 target data from DEM positions and a reference pt
    //tgtArr - (x,y,z,theta,dist2tgt,closesttrjpt,S,S)
    int gid = get_global_id (0);
    
    float diffX, diffY, diffZ, chkDist;
    chkDist = 9999999.0;
    float8 procPt;
    procPt = tgtArr[gid];
    //Check Range 2d
    diffX = fabs(procPt.s0-chkPt.s0);
    diffY = fabs(procPt.s1-chkPt.s1);
    //Dist
    procPt.s4 = half_sqrt(diffX*diffX+diffY*diffY);
    //Write
    tgtArr[gid]=procPt;
}


__kernel void trjRange2(__global float4* trjArr, __global float8* demArr,
                        int trjArrLen, float4 tgtPos)
{
    //INWARD
    //2D !! - dist, z - Calculates closest point in trj array for given point in tgtArr
    //trjArr is the massive array of trajectories - x,y,z,theta
    //demArr is the tasksize: (x,y,z,theta,dist2tgt,closesttrjpt,S,S)
    int gid = get_global_id (0);
    
    float diffX, diffY, diffZ, chkDist, closestDist, closestTheta;
    float8 launchPt;
    int closestIdx;
    launchPt = demArr[gid];
    closestDist = 999999.0;
    closestIdx = 0;
    int i;
    for (i=0; i<trjArrLen; i++){
        //Check Range 2d - target to trajectory points
        diffX = fabs(launchPt.s4-trjArr[i].s0);
        //We need to rebase trj to the launch point height - as it was fired at 0.0 (for now)
        diffZ = fabs(tgtPos.s2-(trjArr[i].s1+launchPt.s2));
        chkDist = half_sqrt(diffX*diffX+diffZ*diffZ);
        if (chkDist < closestDist){
            closestDist = chkDist;
            closestIdx = i;
            closestTheta = trjArr[i].s2;
        }
    }
    //printf("%f, %f :closestDist, closestTheta\n", closestDist, closestTheta);
    //Write theta of closest - TO DEM ARR
    launchPt.s3 = closestTheta;
    launchPt.s5 = closestDist;
    demArr[gid]=launchPt;
}


__kernel void trjIntersect(__global float4* trjArr, __global float8* demArr, __global float* demZArr,
                           int trjArrLen, float8 geoTrans,
                           float2 cl_rasterSize, float4 tgtPos, float maxErr)
{
    //INWARD
    //Intersection after finding theta using trj mode (above - trjRange2)
    //demArr - float8(x,y,z,theta,dist2tgt,closesttrjpt,intersectioncnt,S)
    //trj - (dist, z, theta, S)
    //maxErr - maximum dist between tgt and trjpt for checking if we have a good solution
    //Start point is launch pt
    //GID of tgtArr and demZArr should match...
    int gid = get_global_id (0);
    float8 launchPt;
    launchPt = demArr[gid];
    float floatA, floatB, az, intersect;
    intersect = 0.0f;
    bool bIntersect;
    bIntersect = 0.0;
    //Break out if beyond err
    if (launchPt.s5 > maxErr){
        return;
    }
    //Get the azimuth of this shot
    floatA = tgtPos.s0-launchPt.s0;
    floatB = tgtPos.s1-launchPt.s1;
    az = atan2(floatB, floatA);
    for (int i=0; i<trjArrLen; i++){
        //make sure we skip intersecting if below tgt - trj left over from earlier optim
        if ((trjArr[i].s1+launchPt.s2) < tgtPos.s2){
            continue;
        }
        //Find matching theta labels
        if (trjArr[i].s2 == launchPt.s3) {
            float pxX, pxY, floatD, floatC, fIdx;
            int zIdx;
            //Found one - convert trj pt to UTM at this azimuth
            floatA = launchPt.s0+(trjArr[i].s0 * cos(az)); // X
            floatB = launchPt.s1+(trjArr[i].s0 * sin(az)); // Y
            //Get & check DEM Z at this position - 1D
            pxX = convert_int(fabs((floatA-geoTrans.s0)/geoTrans.s1));
            pxY = convert_int(fabs((floatB-geoTrans.s3)/geoTrans.s5));
            zIdx = (pxY*convert_int(cl_rasterSize.s0))+pxX;
            //Check intersetion - avoiding an if statement here
            //Note rebasing trjZ to where it was launched from
            floatC = convert_float((demZArr[zIdx]-(trjArr[i].s1+launchPt.s2)) > 1.0f);
            demArr[gid].s6+=floatC;
        }
    }
}

//__kernel void trjIntersect(__global float4* trjArr, __global float8* demArr, __global float* demZArr,
//                           int trjArrLen, float8 geoTrans,
//                           float2 cl_rasterSize, float4 tgtPos, float maxErr)
//{
//    //INWARD
//    //Intersection after finding theta using trj mode (above - trjRange2)
//    //demArr - float8(x,y,z,theta,dist2tgt,closesttrjpt,intersectioncnt,S)
//    //trj - (dist, z, theta, S)
//    //maxErr - maximum dist between tgt and trjpt for checking if we have a good solution
//    //Start point is launch pt
//    //GID of tgtArr and demZArr should match...
//    int gid = get_global_id (0);
//    float8 launchPt;
//    launchPt = demArr[gid];
//    float az;
//    //Break out if beyond err
//    if (launchPt.s5 > maxErr){
//        return;
//    }
//    //Get the azimuth of this shot
//    //floatA = tgtPos.s0-launchPt.s0;
//    //floatB = tgtPos.s1-launchPt.s1;
//    az = atan2((tgtPos.s0-launchPt.s0), (tgtPos.s1-launchPt.s1));
//    for (int i=0; i<trjArrLen; i++){
//        //make sure we skip intersecting if below tgt - trj left over from earlier optim
//        if ((trjArr[i].s1+launchPt.s2) < tgtPos.s2){
//            continue;
//        }
//        //Find matching theta labels
//        if (trjArr[i].s2 == launchPt.s3) {
//            float pxX, pxY, floatA, floatB, floatC;
//            int zIdx;
//            //Found one - convert trj pt to UTM at this azimuth
//            floatA = launchPt.s0+(trjArr[i].s0 * cos(az)); // X
//            floatB = launchPt.s1+(trjArr[i].s0 * sin(az)); // Y
//            //Get & check DEM Z at this position - 1D
//            pxX = convert_int(fabs((floatA-geoTrans.s0)/geoTrans.s1));
//            pxY = convert_int(fabs((floatB-geoTrans.s3)/geoTrans.s5));
//            zIdx = (pxY*convert_int(cl_rasterSize.s0))+pxX;
//            //Check intersetion - avoiding an if statement here
//            //Note rebasing trjZ to where it was launched from
//            floatC = convert_float((demZArr[zIdx]-(trjArr[i].s1+launchPt.s2)) > 1.0f);
//            demArr[gid].s6+=floatC;
//        }
//    }
//}





