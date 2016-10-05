//
//  utils.cpp
//  gputest1
//
//  Created by Chris Nicholas on 11/06/2015.
//  Copyright (c) 2015 Chris Nicholas. All rights reserved.
//

#include <iomanip>
#include "utils.h"
#include "position.h"
#include <vector>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>
#include <recShot.h>
#include "gdal_priv.h"
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#include <OpenCL/cl.hpp>
#else
#include <CL/cl.h>
#include <CL/cl.hpp>
#endif



//double Utils::getLinLen(const std::vector<std::shared_ptr<position>>& inTrj) {
//    //Works out the line length (x,y hypot) between trajectory sample points
//    //Work out the hypoteneuse - using position at 2'nd index for safety?
//    double outLen = sqrt((pow(((inTrj.at(2)->x) - (inTrj.at(1)->x)), 2.0) + pow(((inTrj.at(2)->y) - (inTrj.at(1)->y)), 2.0)));
//    return outLen;
//}
//
//
//std::vector<std::shared_ptr<position>> Utils::mkLoS(const std::vector<std::shared_ptr<position>>& worldProf, const double       barrelZ, const double tgtZ) {
//
//        /*
//         * Makes a line of sight to match the points of the input trajectory profile
//         * tgtZ is height of tgt above terrain
//         */
//
//        //TODO:Clean MKLoS up tot make it alot faster - [] instead of .at() etc
//        auto outTLoS = std::vector<std::shared_ptr<position>>();
//        //Get difference between first and last heights & work the amount each step in this profile should be
//        int worldSize = static_cast<int>(worldProf.size());
//
//        //Start and stop heights (observer height above terrain and target height above terrain
//        double obsHeight = worldProf.front()->z + barrelZ;
//        double tgtHeight = worldProf.back()->z + tgtZ;
//
//        //Using std::abs to ensure we get an absolute difference
//        double zStep = std::abs(obsHeight - tgtHeight)/worldProf.size();
//
//        //Create the LoS profile
//        //Add first position without altering
//        outTLoS.push_back( std::make_shared<position>(position { worldProf.front()->x, worldProf.front()->y, obsHeight }));
//
//        for (int i=1; i < (worldSize-1); i++) {
//            //Skipping first and last
//            //Make sure the gradient is correct
//            if (worldProf.front()->z > worldProf.back()->z) {
//                //Downward
//                outTLoS.push_back( std::make_shared<position>(position { worldProf.at(i)->x, worldProf.at(i)->y, obsHeight - (zStep*i)}));
//            }
//            else if ( worldProf.front()->z < worldProf.back()->z ){
//                //Upward
//                outTLoS.push_back( std::make_shared<position>(position { worldProf.at(i)->x, worldProf.at(i)->y, obsHeight + (zStep*i) }));
//            }
//            else {
//                //Flat - pretty unlikely unless in test cases for flat world
//                outTLoS.push_back( std::make_shared<position>(position { worldProf.at(i)->x, worldProf.at(i)->y, obsHeight }));
//            }
//        }
//        //Add last position without altering
//        outTLoS.push_back( std::make_shared<position>(position { worldProf.back()->x, worldProf.back()->y, tgtHeight }));
//        return outTLoS;
//}
//
//
//
//
//std::vector<std::shared_ptr<position>> Utils::buildProfile(const std::shared_ptr <std::vector<std::vector<double>>>& demData,
//                                                                 const std::vector<std::shared_ptr<position>>& inTrj,
//                                                                 const double geoTransform[6], const double& rasterXSize,
//                                                                 const double& rasterYSize) {
//    
//    /*
//     * Makes dem profile based on dem data and inTrj positions
//     */
//    
//    double pxX = 0.0;
//    double pxY = 0.0;
//    double zVal = 0.0;
//    //TODO: Make shared 1 level deeper here
//    auto outProfile = std::vector<std::shared_ptr<position>>();
//    
//    for (auto pos : inTrj) {
//        //Get the x,y and find the corresponding dem value
//        //DEM values are array of raster scanliines (rows) - so Y, X format in vector.
//        
//        //Convert coords to 2d px val
//        coordtoPx2d(pos->x, pos->y, pxX, pxY, geoTransform, rasterXSize, rasterYSize);
//        
//        //Get dem data at px - 2d - Y, X format from earlier load
//        zVal = demData->at(pxY).at(pxX);
//        
//        //record to out data
//        outProfile.push_back(std::make_shared<position>(position { pos->x, pos->y, zVal }));
//        
//    }
//    
//    return outProfile;
//}
//
//std::vector<float> Utils::GDAL2VEC (GDALDataset *poDataset) {
//    
//    /*
//     * Read an entire GDAL DEM into 1D vector shared pointer memory
//     */
//    
//    auto worldZPtr = std::vector<float>();
//    
//    GDALRasterBand  *poBand;
//    poBand = poDataset->GetRasterBand( 1 );
//    int   nXSize = poBand->GetXSize();
//    int   nYSize = poBand->GetYSize();
//    //Options around storing shared ptr's for each row inside the top vector, for now going easy...
//    //TODO: Store each row as a shared ptr - this will be resused and flattened
//    std::vector<double> rowVec;
//    //Using vector instead of CPLMalloc so memory is cleared by compiler...
//    typedef std::vector<float> raster_data_t;
//    raster_data_t scanline(nXSize);
//    
//    //TODO: Improve this to accommodate where we dont have enough memory for the whole ptr...
//    
//    //Start looping and reading rows into the smart ptr
//    //NOTE: This ends in a Y, X array...
//    for (int row = 0; row < nYSize; row++)  {
//        //Read one column from gdal
//        poBand->RasterIO( GF_Read, 0, row, nXSize, 1, &scanline[0], nXSize, 1, GDT_Float32, 0, 0 );
//        //Dump scan line into smart ptr for this row
//        for ( int i = 0; i < nXSize; i++ ) {
//            worldZPtr.push_back(scanline[i]);
//        }
//        //Copy the row vec into the output array
//        //TODO: Implement shared_ptr here
//        
//        //Empty the rowVec for security
//        scanline.clear();
//    }
//    return worldZPtr;
//}

int Utils::optimise(cl_float8* optimArr, int taskSize, float distThreshold, float elevThreshold, float elevMin, float elevMax) {
    //Optimise outside the kernel - eeeek
    //This is the amount the spacing is reduced by each step
    float spacingDiv = 1.5f;
    int workLeft=0;
    for (int i = 0; i < taskSize; i++)  {
        //std::cout << std::setprecision(10) << "Dist: " << optimArr[i].s1 << " / theta: " << optimArr[i].s0 << " / runs: " << optimArr[i].s7 << std::endl;
        if (optimArr[i].s6 < 1.0f){
            workLeft++;
            if (optimArr[i].s1 <= distThreshold) {
                //Close enough - set success flag
                optimArr[i].s6 = 1.0f;
                //Increment the numer of shots - its three per optimisation check in kernel
                optimArr[i].s7+=3.0f;
            } else {
                if (optimArr[i].s0 == elevMin) {
                    if (optimArr[i].s5 > elevThreshold) {
                        //Keep trying - special case for minElev
                        optimArr[i].s5=optimArr[i].s5/spacingDiv;
                        //min
                        optimArr[i].s2=elevMin;
                        //max
                        optimArr[i].s3=elevMin+optimArr[i].s5;
                        //mid
                        optimArr[i].s4=elevMin+(optimArr[i].s5/2.0f);
                    } else {
                        //Seems min elev is closest, we havn't got within the correct distance - cant reach the tgt
                        //(0.0 - running, 1.0 - success, 2.0 - failed)
                        optimArr[i].s6 = 2.0f;
                        //Increment the numer of shots - its three per optimisation check in kernel
                        optimArr[i].s7+=3.0f;
                        //tgtPos.s5 = 0.0f;
                    }
                } else if (optimArr[i].s0 == elevMax) {
                    if (optimArr[i].s5 > elevThreshold) {
                        //Keep trying - special case for maxElev
                        optimArr[i].s5=optimArr[i].s5/spacingDiv;
                        //min
                        optimArr[i].s2=elevMax-optimArr[i].s5;
                        //max
                        optimArr[i].s3=elevMax;
                        //mid
                        optimArr[i].s4=elevMax-(optimArr[i].s5/2.0f);
                    } else {
                        //Seems max elev is closest, we havn't got within the correct distance - cant reach the tgt
                        //(0.0 - running, 1.0 - success, 2.0 - failed)
                        optimArr[i].s6 = 2.0f;
                        //Increment the numer of shots - its three per optimisation check in kernel
                        optimArr[i].s7+=3.0f;
                        //tgtPos.s5 = 0.0f;
                    }
                } else {
                    //Check if we've exceeded elev threshold
                    if (optimArr[i].s5 < elevThreshold) {
                        //Didnt reach dist and have exceeded elev threshold - break at current
                        //(0.0 - running, 1.0 - success, 2.0 - failed)
                        optimArr[i].s6 = 2.0f;
                        //Increment the numer of shots - its three per optimisation check in kernel
                        optimArr[i].s7+=3.0f;
                        //TODO: Set run flag in demarr aswell
                        //tgtPos.s5 = 0.0f;
                    } else {
                        //(0theta, 1mindist, 2elevmin, 3elevmax, 4elevmid, 5spacing, 6optimflag (0.0 - running, 1.0 - success, 2.0 - failed), 7SPARE)
                        //Keep trying
                        optimArr[i].s5=optimArr[i].s5/spacingDiv;
                        //min
                        optimArr[i].s2=optimArr[i].s0-(optimArr[i].s5/2.0f);
                        //Check Boundary
                        if (optimArr[i].s2 < elevMin) {
                            optimArr[i].s2 = elevMin;
                        }
                        //max
                        optimArr[i].s3=optimArr[i].s0+(optimArr[i].s5/2.0f);
                        if (optimArr[i].s3 > elevMax) {
                            optimArr[i].s3 = elevMax;
                        }
                        //mid
                        optimArr[i].s4=optimArr[i].s0;
                        //Incrememnt shots
                        optimArr[i].s7+=3.0f;
                    }
                }
            }
        }
    }
    return workLeft;
}


void Utils::GDAL2FLOAT8 (GDALDataset *poDataset, cl_float8 *worldZPtr) {
    
    /*
     * Read an entire GDAL DEM into float8 vector shared pointer memory!
     * float8(x,y,z,theta,dist2tgt,closesttrjpt,intersectioncnt,S)
     * Also build the x and y coords
     */
    
    
    GDALRasterBand  *poBand;
    poBand = poDataset->GetRasterBand( 1 );
    int   nXSize = poBand->GetXSize();
    int   nYSize = poBand->GetYSize();
    double        adfGeoTransform[6];
    poDataset->GetGeoTransform( adfGeoTransform );
    //Options around storing shared ptr's for each row inside the top vector, for now going easy...
    //TODO: Store each row as a shared ptr - this will be resused and flattened
    std::vector<double> rowVec;
    //Using vector instead of CPLMalloc so memory is cleared by compiler...
    typedef std::vector<float> raster_data_t;
    raster_data_t scanline(nXSize);
    
    //TODO: Improve this to accommodate where we dont have enough memory for the whole ptr...
    
    //Start looping and reading rows into the smart ptr
    //NOTE: This ends in a Y, X array...
    int cnt = 0;
    for (int row = 0; row < nYSize; row++)  {
        //Read one column from gdal
        poBand->RasterIO( GF_Read, 0, row, nXSize, 1, &scanline[0], nXSize, 1, GDT_Float32, 0, 0 );
        //Dump scan line into smart ptr for this row
        for ( int i = 0; i < nXSize; i++ ) {
            double mapX, mapY;
            px2Coord(mapX, mapY, i, row, adfGeoTransform, nXSize, nYSize);
            //initiate dist with high value for reduce and runflag
            //Set no run flag on null data
            //WARNING - USING s6 as TESTING - currently means optim success
            cl_float8 out;
            out = {(float)(mapX), (float)(mapY), (float)(scanline[i]), 0.0, 999999.0f, 0.0f, 0.0f, 0.0f};
            worldZPtr[cnt]=(out);
            cnt++;
        }
        //Copy the row vec into the output array
        //TODO: Implement shared_ptr here
        
        //Empty the rowVec for security
        scanline.clear();
    }
    return;
}

void Utils::GDAL2FLOAT4 (GDALDataset *poDataset, cl_float4 *worldZPtr) {
    
    /*
     * Read an entire GDAL DEM into float8 vector shared pointer memory!
     * float8(x,y,z,theta,dist,runflag,pixelX, pixelY)
     * Also build the x and y coords
     * INCOMPLETE
     */
    
    
    GDALRasterBand  *poBand;
    poBand = poDataset->GetRasterBand( 1 );
    int   nXSize = poBand->GetXSize();
    int   nYSize = poBand->GetYSize();
    double        adfGeoTransform[6];
    poDataset->GetGeoTransform( adfGeoTransform );
    //Options around storing shared ptr's for each row inside the top vector, for now going easy...
    //TODO: Store each row as a shared ptr - this will be resused and flattened
    std::vector<double> rowVec;
    //Using vector instead of CPLMalloc so memory is cleared by compiler...
    typedef std::vector<float> raster_data_t;
    raster_data_t scanline(nXSize);
    
    //TODO: Improve this to accommodate where we dont have enough memory for the whole ptr...
    
    //Start looping and reading rows into the smart ptr
    //NOTE: This ends in a Y, X array...
    int cnt = 0;
    for (int row = 0; row < nYSize; row++)  {
        //Read one column from gdal
        poBand->RasterIO( GF_Read, 0, row, nXSize, 1, &scanline[0], nXSize, 1, GDT_Float32, 0, 0 );
        //Dump scan line into smart ptr for this row
        for ( int i = 0; i < nXSize; i++ ) {
            double mapX, mapY;
            px2Coord(mapX, mapY, i, row, adfGeoTransform, nXSize, nYSize);
            cl_float4 out;
            out = {(float)(mapX), (float)(mapY), (float)(scanline[i]), 0.0};
            worldZPtr[cnt]=(out);
            cnt++;
        }
        //Copy the row vec into the output array
        //TODO: Implement shared_ptr here
        
        //Empty the rowVec for security
        scanline.clear();
    }
    return;
}

std::vector<cl_float4> Utils::genTrjSoup(float thetaStep, float minElev, float maxElev, float timeStep, float muzVel,
                             float mortSigma, float mortMass, float baselineZ) {
    
    //Generate trajectories for all shots with given thetaStep
    //Output is a 1d array
    int idx = 0;
    std::vector<cl_float4> trjVec;
    
    //Loop through launch angles running and recording trajectories
    for (float i=minElev; i<maxElev; i+=thetaStep){
        std::vector<position> outPosVec;
        auto launchPos = std::make_shared<position>( position {0.0, 0.0, baselineZ});
        //Low tgt so we let the trj run for low shots
        //TODO: this can be optimised by seeing the range of Z in DEM
        auto tgtPos = std::make_shared<position>( position {5000.0,5000.0,baselineZ-1000.0});
        recShot::launchAir(outPosVec, launchPos, tgtPos, muzVel, i, timeStep, mortSigma, mortMass);
        //Load to intermediate vector
        for (auto pos : outPosVec) {
            float dist = sqrt(pos.x*pos.x+pos.y*pos.y);
            trjVec.push_back({dist, float(pos.z), i, 0.0});
            idx++;
        }
    }
    return trjVec;
}

//void Utils::getzVal (std::shared_ptr<position>& inPos, const double geoTransform[6],
//                           const std::shared_ptr <std::vector<std::vector<double>>>& demData,
//                           const double& rasterXSize, const double& rasterYSize) {
//    
//    
//    /*Updates in-position with DEM Z Value for the given position - presumes this has already been checked to be on DEM...
//     * Input: World coords in position (z val is barrel height)
//     * Output: world coords in position (dem Z + barrel height)
//     */
//    
//    //Convert coords to 2d px val
//    double pxX = 0.0;
//    double pxY = 0.0;
//    Utils::coordtoPx2d(inPos->x, inPos->y, pxX, pxY, geoTransform, rasterXSize, rasterYSize);
//    //Get dem data at px - 2d - Y, X format from earlier load
//    inPos->z = demData->at(pxY).at(pxX);
//}

//std::vector<cl_float4> Utils::GDAL2VEC2D (GDALDataset *poDataset) {
//    
//    /*
//     * Read an entire GDAL DEM into 2D vector shared pointer memory
//     */
//    
//    std::vector<cl_float4> worldZPtr = std::vector<cl_float4>();
//    double        adfGeoTransform[6];
//    poDataset->GetGeoTransform( adfGeoTransform );
//    int rasterXSize = poDataset->GetRasterXSize();
//    int rasterYSize = poDataset->GetRasterYSize();
//    
//    GDALRasterBand  *poBand;
//    poBand = poDataset->GetRasterBand( 1 );
//    int   nXSize = poBand->GetXSize();
//    int   nYSize = poBand->GetYSize();
//    //Options around storing shared ptr's for each row inside the top vector, for now going easy...
//    //TODO: Store each row as a shared ptr - this will be resused and flattened
//    std::vector<float> rowVec;
//    cl_float4 pos;
//    //Using vector instead of CPLMalloc so memory is cleared by compiler...
//    typedef std::vector<float> raster_data_t;
//    raster_data_t scanline(nXSize);
//    
//    //TODO: Improve this to accommodate where we dont have enough memory for the whole ptr...
//    
//    //Start looping and reading rows into the smart ptr
//    //NOTE: This ends in a Y, X array...
//    for (int row = 0; row < nYSize; row++)  {
//        //Read one column from gdal
//        poBand->RasterIO( GF_Read, 0, row, nXSize, 1, &scanline[0], nXSize, 1, GDT_Float32, 0, 0 );
//        //Dump scan line into smart ptr for this row
//        for ( int i = 0; i < nXSize; i++ ) {
//            px2CoordCL(pos, i, row, adfGeoTransform, rasterXSize, rasterYSize);
//            pos.z = scanline[i];
//            worldZPtr.push_back(pos);
//        }
//        //Copy the row vec into the output array
//        //TODO: Implement shared_ptr here
//        //worldZPtr.push_back(rowVec);
//        //Empty the rowVec for security
//        rowVec.clear();
//        scanline.clear();
//    }
//    return worldZPtr;
//}
//
//
//std::vector<std::shared_ptr<position>> Utils::buildLaunchVec(const double geoTransform[6],
//                                                                   const int& rasterSizeX,
//                                                                   const int& rasterSizeY,
//                                                                   const std::shared_ptr <std::vector<std::vector<double>>>& demData) {
//    
//    /*
//     * Builds an array of Launch positions based on input DEM header and skip value
//     * Input: GDAL Header geotransform, integer skipValue
//     * Output: Vector<std:shared_ptr<position>>
//     */
//    
//    //Create the output vector
//    std::vector<std::shared_ptr<position>> outVec;
//    auto pos = std::make_shared<position>(position {0,0,0});
//    double mapX = 0.0;
//    double mapY = 0.0;
//    
//    //Output launch positions to a text file - Working
//    //std::ofstream output("/Users/dusted-ipro/geodata/launchvec3.csv");
//    
//    //Loop through Cols and Rows
//    for (int col=0; col < rasterSizeX; col++) {
//        for (int row=0; row < rasterSizeY; row++){
//            //convert px to coord
//            Utils::px2Coord(mapX, mapY, col, row, geoTransform, rasterSizeX, rasterSizeY);
//            //Add 0.5 of pixel size for centre of pixel
//            pos->x = mapX+(0.5*geoTransform[1]);
//            pos->y = mapY-(0.5*geoTransform[5]);
//            //Get the Height
//            getzVal (pos, geoTransform, demData, rasterSizeX, rasterSizeY);
//            //Add the position to out array
//            outVec.push_back(std::make_shared<position>( position {pos->x, pos->y, pos->z}));
//            //Output launch positions to a text file - Working
//            //output << std::setprecision(10) << pos->x << "," << pos->y << "," << pos->z << std::endl;
//            pos->z = 0.0;
//        }
//    }
//    //Output launch positions to a text file - Working
//    //output.close();
//    return outVec;
//}

void Utils::coordtoPx2d(const double& x0, const double& y0, double& pxX, double& pxY,
                              const double geoTransform[6], const double& rasterXSize, const double& rasterYSize) {
    
    /*
     * Convert world coords to DEM pixel coord
     */
    
    double diffX = 0.0;
    double diffY = 0.0;
    //Get range from origin (top left) and make sure we map to the centre of a pixel
    diffX = x0-geoTransform[0];
    diffY = y0-geoTransform[3];
    pxX = std::abs(diffX/geoTransform[1]);
    pxY = std::abs(diffY/geoTransform[5]);
}

void Utils::px2Coord(double& mapX, double& mapY, const double& pxX, const double& pxY,
                           const double geoTransform[6], const double& rasterXSize, const double& rasterYSize) {
    
    /*
     * Convert DEM pixel index into world coords (as projection)
     */
    
    //Do affine transform for world coords
    //+0.5*std::abs(adfGeoTransform[1]);
    //-0.5*std::abs(adfGeoTransform[5]);
    //mapX = geoTransform[0]+(pxX*geoTransform[1])+pxY*std::abs(geoTransform[2]);
    //mapY = geoTransform[3]+(pxX*geoTransform[4])+pxY*std::abs(geoTransform[5]);
    mapX = geoTransform[0]+(pxX*geoTransform[1])+pxY*geoTransform[2];
    mapY = geoTransform[3]+(pxX*geoTransform[4])+pxY*geoTransform[5];
}

void Utils::px2CoordCL(cl_float4& map, const int& pxX, const int& pxY,
                     const double geoTransform[6], const double& rasterXSize, const double& rasterYSize) {
    
    /*
     * Convert DEM pixel index into world coords (as projection)
     */
    
    //Do affine transform for world coords
    //+0.5*std::abs(adfGeoTransform[1]);
    //-0.5*std::abs(adfGeoTransform[5]);
    //mapX = geoTransform[0]+(pxX*geoTransform[1])+pxY*std::abs(geoTransform[2]);
    //mapY = geoTransform[3]+(pxX*geoTransform[4])+pxY*std::abs(geoTransform[5]);
    map.x = geoTransform[0]+(pxX*geoTransform[1])+pxY*geoTransform[2];
    map.y = geoTransform[3]+(pxX*geoTransform[4])+pxY*geoTransform[5];
}

GDALDataset* Utils::openDem(const char* demFName) {
    
    /*
     * Open a raster dataset using GDAL
     */
    
    //Value equal to the dataset pointed to by poDataset
    GDALDataset  *poDataset;
    GDALAllRegister();
    poDataset = (GDALDataset *) GDALOpenShared( demFName, GA_ReadOnly);
    
    //TODO: Clearly need to do something better here
    if ( poDataset == NULL ) {
        std::cout << "ERROR - Input data was NULL" << std::endl;
        return poDataset;
    }
    else {
        return poDataset;
    }
}



int Utils::singlefOut(const std::vector<cl_float4>& traj, const char* fName) {
    //Safer file output for testing
    
    std::ofstream output(fName);
    
    if (!output) { // check the file opened OK
        std::cerr << "error: open file for output failed!" << std::endl;
        return 1;
    }
    //cl_float4 first = traj.front()->z;
    
    //Sizes should all be same
    for (unsigned int i=0; i<traj.size(); i++){
        // the ofstream object replaces std::cout here, then its the same use
        output << traj.at(i).z << std::endl;
    }
    
    /*
     * Single Traj Output - x,y,z
     * Working
     for (auto pos : inTrj) {
     // the ofstream object replaces std::cout here, then its the same use
     output << pos->x << "," << pos->y << "," << pos->z << std::endl;
     }
     */
    output.close(); // must close file
    
    return 0;
    
}

int Utils::writeRasterOut (GDALDataset *poDataset, const char* outFName, cl_float8* optimDEMOut) {
    /*
     * Dump results out to a geotiff
     * Input: pointer to DEM data - not thinned version
     * Output: 2 band tiff
     */
    
    const char *pszFormat = "GTiff";
    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    char *pszSRS_WKT = NULL;
    GDALRasterBand *poBand;
    int xSize = poDataset->GetRasterXSize();
    int ySize = poDataset->GetRasterYSize();
    int currRow = 0;
    int currCol = 0;
    double        adfGeoTransform[6];
    
    //Out Gdal dataset
    GDALDataset *poDstDS;
    char **papszOptions = NULL;
    
    //Create the driver for geotiff - 1 band
    poDstDS = poDriver->Create( outFName, xSize, ySize, 1, GDT_Float32, papszOptions );
    //Get the geotrans to the DEM geo trans
    poDataset->GetGeoTransform( adfGeoTransform );
    poDstDS->SetGeoTransform( adfGeoTransform );
    poDstDS->SetProjection( poDataset->GetProjectionRef() );
    
    //Outdata array
    float *abyRaster = new float[ySize*xSize];
    //unsigned char abyRaster[ySize*xSize];
    
    //Output band 1 - BloS
    poBand = poDstDS->GetRasterBand(1);
    //NOTE: Looping across the smart ptr does not work - has to be other way round
    for (int row=0; row < ySize; row++){
        for (int col=0; col < xSize; col++){
            //WARNING - USING s6 as TESTING - currently means theta
            abyRaster[row*xSize+col] = optimDEMOut[row*xSize+col].s3;
        }
    }
    
    //Dump the data to raster
    poBand->RasterIO( GF_Write, 0, 0, xSize, ySize, abyRaster, xSize, ySize, GDT_Float32, 0, 0 );
    
    /* Once we're done, close properly the dataset */
    GDALClose( (GDALDatasetH) poDstDS );
    //Delete data ptr
    delete[] abyRaster;
    return 0;
}

