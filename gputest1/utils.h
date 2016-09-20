//
//  utils.h
//  gputest1
//
//  Created by Chris Nicholas on 11/06/2015.
//  Copyright (c) 2015 Chris Nicholas. All rights reserved.
//

#ifndef __gputest1__utils__
#define __gputest1__utils__

#include <stdio.h>
#include "position.h"
#include <vector>
#include <memory>
#include "gdal_priv.h"
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#include <OpenCL/cl.hpp>
#else
#include <CL/cl.h>
#include <CL/cl.hpp>
#endif


class Utils {
public:
    
//    static std::vector<std::shared_ptr<position>> mkLoS(const std::vector<std::shared_ptr<position>>& worldProf, const double barrelZ, const double tgtZ);
//    
//    static double getLinLen(const std::vector<std::shared_ptr<position>>& inTrj);
//    
//    static std::vector<std::shared_ptr<position>> buildProfile(const std::shared_ptr <std::vector<std::vector<double>>>& demData,
//                                                               const std::vector<std::shared_ptr<position>>& inTrj,
//                                                               const double geoTransform[6], const double& rasterXSize,
//                                                               const double& rasterYSize);
//    
//    static void getzVal (std::shared_ptr<position>& inPos, const double geoTransform[6],
//                         const std::shared_ptr <std::vector<std::vector<double>>>& demData,
//                         const double& rasterXSize, const double& rasterYSize);
//    
//    static std::vector<std::shared_ptr<position>> buildLaunchVec(const double geoTransform[6],
//                                                                 const int& rasterSizeX,
//                                                                 const int& rasterSizeY,
//                                                                 const std::shared_ptr <std::vector<std::vector<double>>>& demData);

    static int optimise (cl_float8 *optimArray, int taskSize, float distThreshold, float elevThreshold, float elevMin, float elevMax);
    
    //static std::vector<float> GDAL2VEC (GDALDataset *poDataset);
    
    static void GDAL2FLOAT8 (GDALDataset *poDataset, cl_float8 *worldZPtr);
    
    static std::vector<cl_float4> GDAL2VEC2D (GDALDataset *poDataset);
    
    static void coordtoPx2d(const double& x0, const double& y0, double& pxX, double& pxY,
                            const double geoTransform[6], const double& rasterXSize, const double& rasterYSize);
    
    static void px2Coord(double& mapX, double& mapY, const double& pxX, const double& pxY,
                         const double geoTransform[6], const double& rasterXSize, const double& rasterYSize);
    
    static void px2CoordCL(cl_float4& map, const int& pxX, const int& pxY,
                                  const double geoTransform[6], const double& rasterXSize, const double& rasterYSize);
    
    static GDALDataset* openDem(const char* demFName);
    
    static int singlefOut(const std::vector<cl_float4>& traj, const char* fName);
    
};

#endif /* defined(__gputest1__utils__) */