//
//  main.cpp
//  gputest1
//
//  Created by Chris Nicholas on 26/05/2015.
//  Copyright (c) 2015 Chris Nicholas. All rights reserved.
//

//#define PROGRAM_FILE "simplekern.cl"
#define __CL_ENABLE_EXCEPTIONS

#include <iostream>
#include <fstream>
#include <iterator>
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#include <OpenCL/cl.hpp>
#else
#include <CL/cl.h>
#include <CL/cl.hpp>
#endif
//#include <OpenCL/opencl.h>

#include <time.h>
#include <cmath>
//#include <position.h>
#include "gdal_priv.h"
#include "utils.h"
//#include <worldParams.h>
//#include <motionUtils.h>
#include <iomanip>
//#include <recShot.h>


int main(int argc, const char * argv[]) {
    
    //Files, vars etc
    char* demFName = "/Users/dusted-dstl/Documents/geodata/mount_chip.tif";
    //char* demFName = "/Users/dusted-dstl/Documents/geodata/mount.dem";
    char* kernelsSrc = "/Users/dusted-dstl/Documents/xcodeworkspace/gpublos/gputest1/simplekern.cl";
    //Linux
    //char* demFName = "/home/ec2-user/gputest1/mount.dem";
    //char* kernelsSrc = "/home/ec2-user/gputest1/simplekern.cl";
    float tgtX = 559783.0;
    float tgtY = 5119823.0;
    double tgtPxX = 0.0;
    double tgtPxY = 0.0;
    int skipVal = 1000;
    float minDist = 600.0;
    float maxDist = 4500.0;
    float gravity = 9.82;
    //TODO: Alter this based on the altitude
    float airDens = 1.225;
    //TODO: Allow setting of this based on weapon vars
    //Using Sreamlined body coefficient
    float dragCoef = 0.04; //Coef for sphere
    float pi = 3.1415926; // pi
    float calibre = 0.081;
    float fArea = pi*(pow((0.5*calibre), 2.0));
    float mortSigma = dragCoef*fArea*0.5*airDens;
    float muzzVel = 225.0;
    float mortMass = 4.2;
    float gForce = mortMass * gravity;
    //float launchZ = 0.0;
    float tgtZ = 0.0;
    char* fName = "/Users/dusted-dstl/Documents/geodata/gpu_out1.csv";
    float minElev = 45.0;
    float maxElev = 85.0;
    //Change and threshold are the same...
    float minElevChange = 0.01;
    float elevThreshold = 0.01;
    float distThreshold = 20.0;
    float thetaStep = 0.01;
    float timeStep = 0.05;
    
    
    
    // For catching cl errors
    cl_int err;
    //get all platforms (drivers)
    std::vector<cl::Platform> platforms;
    std::vector<cl::Device> devices;
    std::vector<cl::Kernel> kernels;
    cl::Platform::get(&platforms);
    if(platforms.size()==0){
        std::cout<<" No platforms found. Check OpenCL installation!\n";
        exit(1);
    }
    std::cout << " All Platforms follow " << std::endl;
    for (int i=0; i<platforms.size(); i++) {
        std::cout << platforms[i].getInfo<CL_PLATFORM_NAME>() << std::endl;
    }
    cl::Platform default_platform=platforms[0];
    std::cout << "Querying platform: "<<default_platform.getInfo<CL_PLATFORM_NAME>()<<"\n";
    
    //get default device of the default platform
    default_platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
    if(devices.size()==0){
        std::cout<<" No devices found. Check OpenCL installation!\n";
        exit(1);
    }
    std::cout << " All Devices follow " << std::endl;
    for (int i=0; i<devices.size(); i++) {
        std::cout << "Device Name: " << devices[i].getInfo<CL_DEVICE_NAME>() << std::endl;
        std::cout << "Device Global Mem max size (mb): " << devices[i].getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>()/1000000 << std::endl;
        std::cout << "Device Local Mem max size: " << devices[i].getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << std::endl;
    }
    
    
    // create platform
    cl::Platform::get(&platforms);
    //CPU Implementation
    //platforms[0].getDevices(CL_DEVICE_TYPE_CPU, &devices);
    //GPU Implementation
    platforms[0].getDevices(CL_DEVICE_TYPE_GPU, &devices);
    //Query device info for memory size
    std::cout << devices[0].getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << std::endl;
    
    //Change this if want to use GPU
    cl::Context context = cl::Context(devices);
    
    //Create Command queue
    //CHECK THIS FOR GPU
    cl::CommandQueue queue(context, devices[0]);
    
    // load opencl source
    std::ifstream cl_file(kernelsSrc);
    std::string cl_string(std::istreambuf_iterator<char>(cl_file), (std::istreambuf_iterator<char>()));
    cl::Program::Sources source(1, std::make_pair(cl_string.c_str(), cl_string.length()+1));
    
    // create program
    cl::Program program(context, source);
    
    try {
        
        // compile opencl source
        program.build(devices);
        std::cout << "built program" << std::endl;
        
    }
    catch (cl::Error &err) {
        //Get the build log for the first device
        std::cerr << "Building failed, " << err.what() << "(" << err.err() << ")"
        << "\nRetrieving build log\n"
        << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0])
        << "\n";
        return -1;
        
    }
    
    std::cout << "Building indata" << std::endl;
    
    //Open the DEM and load data into memory
    //GDAL Vars
    double        adfGeoTransform[6];
    //float rasterXSize = 0.0;
    //float rasterYSize = 0.0;
    GDALDataset  *poDataset;
    poDataset = Utils::openDem(demFName);
    
    //Get geotransform info ready
    poDataset->GetGeoTransform( adfGeoTransform );
    float gt0 = float(adfGeoTransform[0]);
    float gt1 = float(adfGeoTransform[1]);
    float gt2 = float(adfGeoTransform[2]);
    float gt3 = float(adfGeoTransform[3]);
    float gt4 = float(adfGeoTransform[4]);
    float gt5 = float(adfGeoTransform[5]);
    cl_float8 cl_geoTrans = {float(adfGeoTransform[0]), float(adfGeoTransform[1]), float(adfGeoTransform[2]), float(adfGeoTransform[3]),
        float(adfGeoTransform[4]), float(adfGeoTransform[5])};
    std::cout << gt0 << " // " << adfGeoTransform[0] << std::endl;
    int rasterXSize = poDataset->GetRasterXSize();
    int rasterYSize = poDataset->GetRasterYSize();
    cl_float2 cl_rasterSize = {float(rasterXSize), float(rasterYSize)};
    int taskSize = rasterXSize*rasterYSize;
    std::cout << "rasterSizeX: " << rasterXSize << std::endl;
    std::cout << "rasterSizeY: " << rasterYSize << std::endl;
    
    GDALRasterBand  *poBand;
    poBand = poDataset->GetRasterBand( 1 );
    int   nXSize = poBand->GetXSize();
    int   nYSize = poBand->GetYSize();
    
    //std::vector<float> demSP = Utils::GDAL2VEC (poDataset);
    //Load DEM into cl_float8:
    //float8(x,y,z,theta,dist,runflag,blos,tlos)
    cl_float8 *h_demSP = new cl_float8[taskSize];
    cl_float8 *output_demSP = new cl_float8[taskSize];
    Utils::GDAL2FLOAT8(poDataset, h_demSP);
    
    //Load DEMZ's into single Arr
    float *h_demZArr = new float[taskSize];
    for (int i=0; i<taskSize; i++){
        //Only want the z's
        h_demZArr[i] = h_demSP[i].s2;
    }
    
    
    GDALClose( (GDALDatasetH) poDataset );
    std::cout << "Size of DEM in memory Ref (kb): " << (taskSize*sizeof(float))/1000 << std::endl;
    std::cout << "Size of DEM: " << taskSize << " tasksize: " << taskSize << std::endl;
    std::cout << "Size of DEM as Float8 (mb): " << taskSize*sizeof(cl_float8)/1000000 << std::endl;
    std::cout << "Will the whole DEM fit into Global mem @ float8?: " << std::endl;
    bool chkGlobalMem = 0;
    for (int i=0; i<devices.size(); i++) {
        chkGlobalMem = taskSize*sizeof(cl_float8) < devices[i].getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();
        std::cout << "Device Name: " << devices[i].getInfo<CL_DEVICE_NAME>() << std::endl;
        std::cout << "Sizes (dem, global mem) & Chk: " << taskSize*sizeof(cl_float8) << " // " << devices[i].getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << " // " << chkGlobalMem << std::endl;
    }
    if (chkGlobalMem == 0){
        std::cout << "Not enough Global memory for the DEM on this device - aborting " << std::endl;
        return 1;
    }
    
//    //Build the array of possible launch angles
//    std::vector<float> thetaArr;
//    for (float i=minElev; i<maxElev; i+=thetaStep){
//        thetaArr.push_back(i);
//    }
//    float *h_thetaArr = new float[thetaArr.size()];
//    for (int i=0; i<thetaArr.size(); i++){
//        h_thetaArr[i] = thetaArr.at(i);
//    }
    
//    //Build the array of output distances for each theta
//    float *h_distArr = new float[thetaArr.size()];
//    std::cout << "Parrellel Shots per DEM Px: " << thetaArr.size()<< std::endl;
//    std::cout << "Total shots: " << thetaArr.size()*taskSize<< std::endl;
    
    //Build Optimisation Data Arr (theta, mindist, elevmin, elevmax, elevmid, spacing, optimflag (0.0 - running, 1.0 - success, 2.0 - failed), SPARE)
    cl_float8 *h_optimArr = new cl_float8[taskSize];
    for (int i=0; i<taskSize; i++){
        cl_float8 out = {0.0f, 99999.0f, minElev, maxElev, minElev+((maxElev-minElev)/2.0f), 10.0f, 0.0f, 0.0f};
        h_optimArr[i] = out;
    }
    
    //clock_t tt = clock();
    //clock_t tt_real = clock();
    
    try {
        std::cout << "Loading Buffers " << std::endl;
        //THESE TWO FOR CPU - then COMMENT CPU
        //Create a buffer object to this memory
        cl::Buffer d_inBuffDEMFull = cl::Buffer(context, CL_MEM_READ_WRITE, taskSize*sizeof(cl_float8), NULL, &err);
//        cl::Buffer d_outBuff = cl::Buffer(context, CL_MEM_READ_WRITE, taskSize*sizeof(float), NULL, &err);
//        cl::Buffer d_thetaArr = cl::Buffer(context, CL_MEM_READ_ONLY, thetaArr.size()*sizeof(float), NULL, &err);
//        cl::Buffer d_distArr = cl::Buffer(context, CL_MEM_READ_WRITE, thetaArr.size()*sizeof(float), NULL, &err);
        cl::Buffer d_demZArr = cl::Buffer(context, CL_MEM_READ_ONLY, taskSize*sizeof(float), NULL, &err);
        cl::Buffer d_optimArr = cl::Buffer(context, CL_MEM_READ_WRITE, taskSize*sizeof(cl_float8), NULL, &err);
        
        //Out data on the host
//        float *h_output = new float[taskSize];
        
        std::cout << "Writing Buffers " << std::endl;
        //Write Buffers to device
        //Theta Arr
//        queue.enqueueWriteBuffer(d_thetaArr,CL_TRUE,0,thetaArr.size()*sizeof(float), h_thetaArr);
        //DEM Arr
        queue.enqueueWriteBuffer(d_inBuffDEMFull,CL_TRUE,0,taskSize*sizeof(cl_float8), h_demSP);
        //DEM Z Arr
        queue.enqueueWriteBuffer(d_demZArr,CL_TRUE,0,taskSize*sizeof(float), h_demZArr);
        //Optimise Arr
        queue.enqueueWriteBuffer(d_optimArr,CL_TRUE,0,taskSize*sizeof(cl_float4), h_optimArr);
        
        
        //WORKS
        //Build Kernels
        //cl::KernelFunctor rangeChk(cl::Kernel(program,"rangeChk"),queue,cl::NullRange,cl::NDRange(taskSize),cl::NullRange);
        //cl::KernelFunctor shotOpt(cl::Kernel(program,"shotOpt"),queue,cl::NullRange,cl::NDRange(thetaArr.size()),cl::NullRange);
        //cl::KernelFunctor thetaDistReduce(cl::Kernel(program,"thetaDistReduce"),queue,cl::NullRange,cl::NDRange(thetaArr.size()),cl::NullRange);
        
        std::cout << "Building kernels " << std::endl;
        cl::Kernel rangeChk = cl::Kernel(program, "rangeChk", &err);
//        cl::Kernel shotOpt = cl::Kernel(program, "shotOpt", &err);
//        cl::Kernel shotOptThetaFast = cl::Kernel(program, "shotOptThetaFast", &err);
//        cl::Kernel shotOptDem = cl::Kernel(program, "shotOptDem", &err);
//        cl::Kernel shotOptDemFast = cl::Kernel(program, "shotOptDemFast", &err);
//        cl::Kernel thetaDistReduce = cl::Kernel(program, "thetaDistReduce", &err);
        cl::Kernel demIntersectFast = cl::Kernel(program, "demIntersectFast", &err);
//        cl::Kernel optimise = cl::Kernel(program, "optimise", &err);
        cl::Kernel optEasy = cl::Kernel(program, "optEasy", &err);
        
        
        //Dummy target position
        cl_float4 tgtPos = {tgtX, tgtY, 3889.0f, 0.0};
        
        //WORKS - Kernel Func
        //CL - Firstly thin the data on device by range
        //rangeChk(d_inBuffDEMFull, tgtPos, minDist, maxDist);
        
        std::cout << "Enqueue Args - rangeChk " << std::endl;
        err = rangeChk.setArg(0, d_inBuffDEMFull);
        err = rangeChk.setArg(1, tgtPos);
        err = rangeChk.setArg(2, minDist);
        err = rangeChk.setArg(3, maxDist);
        queue.finish();
        
        std::cout << "Execute rangeChk " << std::endl;
        err = queue.enqueueNDRangeKernel(rangeChk, cl::NullRange, cl::NDRange(taskSize), cl::NullRange, NULL, NULL);
        queue.finish();
        
        //Optimisation and reduction
        double tgtPxX = 0.0;
        double tgtPxY = 0.0;
        Utils::coordtoPx2d(tgtPos.x, tgtPos.y, tgtPxX, tgtPxY, adfGeoTransform, rasterXSize, rasterYSize);
        int p = (int)(tgtPxY*rasterXSize+tgtPxX);
        tgtPos.z = h_demSP[p].z;
        //cl_float4 laPos = {0.0, 0.0, 0.0, 0.0};
        float minTheta = 0.0;
        
        //Works - DEM Fast
        //std::cout << "Shot Args" << std::endl;
        //err = shotOptDemFast.setArg(0, d_inBuffDEMFull);
        //err = shotOptDemFast.setArg(1, tgtPos);
        //err = shotOptDemFast.setArg(2, timeStep);
        //queue.finish();
        
//        //Works - DEM Fast
        std::cout << "Intersect Args" << std::endl;
        err = demIntersectFast.setArg(0, d_inBuffDEMFull);
        err = demIntersectFast.setArg(1, tgtPos);
        err = demIntersectFast.setArg(2, timeStep);
        err = demIntersectFast.setArg(3, d_demZArr);
        err = demIntersectFast.setArg(4, cl_geoTrans);
        err = demIntersectFast.setArg(5, cl_rasterSize);
        queue.finish();
        
        
//        //WORKS - Theta Fast
//        std::cout << "Enqueue Args - shotOptThetaFast " << std::endl;
//        //Build the static input args - ShotOpt
//        err = shotOptThetaFast.setArg(0, d_thetaArr);
//        err = shotOptThetaFast.setArg(1, d_distArr);
//        err = shotOptThetaFast.setArg(2, tgtPos);
//        err = shotOptThetaFast.setArg(4, timeStep);
//        queue.finish();
//        
//        //WORKS - ShotOptDemFast
//        std::cout << "Enqueue Args - shotOptDemFast " << std::endl;
//        //Build the static input args - ShotOpt
//        err = shotOptDemFast.setArg(0, d_inBuffDEMFull);
//        err = shotOptDemFast.setArg(1, tgtPos);
//        err = shotOptDemFast.setArg(2, timeStep);
//        err = shotOptDemFast.setArg(3, d_optimArr);
//        queue.finish();
//        
//        std::cout << "Enqueue Args - thetaReduce" << std::endl;
//        //Build the static input args - Dist Reduce
//        err = thetaDistReduce.setArg(0, d_thetaArr);
//        err = thetaDistReduce.setArg(1, d_distArr);
//        err = thetaDistReduce.setArg(2, d_inBuffDEMFull);
//        queue.finish();
//        
//        std::cout << "Enqueue Args - optimise" << std::endl;
//        //Build the static input args - Dist Reduce
//        err = optimise.setArg(0, d_inBuffDEMFull);
//        err = optimise.setArg(1, d_optimArr);
//        err = optimise.setArg(2, tgtPos);
//        queue.finish();
        
        std::cout << "Enqueue Args - optEasy" << std::endl;
        //Build the static input args - Dist Reduce
        err = optEasy.setArg(0, d_inBuffDEMFull);
        err = optEasy.setArg(1, d_optimArr);
        err = optEasy.setArg(2, tgtPos);
        err = optEasy.setArg(3, timeStep);
        err = optEasy.setArg(4, distThreshold);
        err = optEasy.setArg(5, elevThreshold);
        err = optEasy.setArg(6, minElev);
        err = optEasy.setArg(7, maxElev);
        queue.finish();
        
//        //Read the processed DEM data back out to get the number of range thinned shots and to loop for theta optim input
//        err = queue.enqueueReadBuffer(d_inBuffDEMFull, CL_TRUE, 0, taskSize*sizeof(cl_float8), h_demSP, NULL, NULL);
//        // wait for completion
//        queue.finish();
//        
//        int inRngCnt = 0;
//        for (int i=0; i<taskSize; i++){
//            if (h_demSP[i].s5 == 1.0f){
//                inRngCnt++;
//            }
//        }
//
//        std::cout << "Num Valid Shots after range thin: " << inRngCnt << std::endl;
        
        //std::clock_t start;
        double duration;
        
        //start = std::clock();
        
        //CL - Loop the optimisation brute force - Whole DEM
        std::cout << "Starting Optimisation" << std::endl;
        int loops = 5;
        
        //DEM Parrellel Version
        int runCnt = 0;
        int skipCnt = 0;
        
        //Easy Shot Optimisation
        err = queue.enqueueNDRangeKernel(optEasy, cl::NullRange, cl::NDRange(taskSize), cl::NullRange, NULL, NULL);
        queue.finish();
        
        std::cout << "Done Optimisation, copying data" << std::endl;
        
        //Brute force
        //for (int i=0; i<loops; i++){
            
            //DEM Parallel
            //err = queue.enqueueNDRangeKernel(shotOptDem, cl::NullRange, cl::NDRange(taskSize), cl::NullRange, NULL, NULL);
            //Fire Shot
            //err = queue.enqueueNDRangeKernel(shotOptDemFast, cl::NullRange, cl::NDRange(taskSize), cl::NullRange, NULL, NULL);
            //queue.finish();


            
            //WORKS - Theta Parallel
//            if (h_demSP[i].s5 == 1.0f){
//                runCnt++;
//                //Fire Shots
//                err = shotOptThetaFast.setArg(3, h_demSP[i]);
//                err = queue.enqueueNDRangeKernel(shotOptThetaFast, cl::NullRange, cl::NDRange(thetaArr.size()), cl::NullRange, NULL, NULL);
//                queue.finish();
//                
//                //Find Closest - write to DEM Arr on device
//                err = thetaDistReduce.setArg(3, i);
//                err = queue.enqueueNDRangeKernel(thetaDistReduce, cl::NullRange, cl::NDRange(thetaArr.size()), cl::NullRange, NULL, NULL);
//                queue.finish();
//                
//            } else {
//                skipCnt++;
//            }
//            
//            //Some prints
//            if (i%1000 == 0){
//                std::cout << "Runs / Skips / Total: " << runCnt << " / " << skipCnt << " / " << i << std::endl;
//            }
            
            
            //WORKS but Slow - kernel Functors
            //Fire Shots
            //Temporary tgt pos for testing
            //float8(x,y,z,theta,dist,runflag,blos,tlos)
            //cl_float8 pos = {569786.0, 5130048.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
            //shotOpt(d_thetaArr, d_distArr, tgtPos, pos, timeStep);
            //shotOpt(d_thetaArr, d_distArr, tgtPos, h_demSP[i], timeStep);
            //Write the closest into the DEM arr
            //thetaDistReduce(d_thetaArr, d_distArr, d_inBuffDEMFull, i);
        //}
        
        
        //Intersect - Now we pivot and parallelize based on the DEM
        //err = queue.enqueueNDRangeKernel(shotOptDem, cl::NullRange, cl::NDRange(taskSize), cl::NullRange, NULL, NULL);
        std::cout << "Done Theta optimisation, running intersect... " << std::endl;
        err = queue.enqueueNDRangeKernel(demIntersectFast, cl::NullRange, cl::NDRange(taskSize), cl::NullRange, NULL, NULL);
        queue.finish();
        
        
        //Read DEM Arr back out
        err = queue.enqueueReadBuffer(d_inBuffDEMFull, CL_TRUE, 0, taskSize*sizeof(cl_float8), h_demSP, NULL, NULL);
        err = queue.enqueueReadBuffer(d_optimArr, CL_TRUE, 0, taskSize*sizeof(cl_float8), h_optimArr, NULL, NULL);
        // wait for completion
        queue.finish();
//
//        std::cout << "Complete - checking results" << std::endl;
//        
//        //Count the results
//        int blosIntersectCnt = 0;
//        int rangeSkip = 0;
//        int optimSuccess = 0;
//        int optimFail = 0;
//        int optimPass = 0;
//        int totalRuns = 0;
//        for (int i=0; i<taskSize; i++) {
//            if (h_demSP[i].s5 == 1.0f) {
//                if (h_optimArr[i].s5 == 2.0f) {
//                    optimFail++;
//                } else if (h_optimArr[i].s5 == 0.0f) {
//                    optimPass++;
//                } else {
//                    optimSuccess++;
//                    if (h_demSP[i].s6 == 1.0f) {
//                        blosIntersectCnt++;
//                    }
//                    if (h_demSP[i].s5 == 0.0f) {
//                        rangeSkip++;
//                    }
//                    totalRuns+=h_optimArr[i].s7;
//                    //std::cout << std::setprecision(10) << "Tgt Dist: " << h_demSP[i].s4 << " / theta: " << h_demSP[i].s3 << " / runs: " << h_optimArr[i].s7 << std::endl;
//                }
//            }
//        }
//        
//        std::cout << "Total Pixels: " << taskSize << std::endl;
//        std::cout << "Total optim success: " << optimSuccess << std::endl;
//        std::cout << "Total optim runs (total shots): " << totalRuns << std::endl;
//        std::cout << "Total optim fail: " << optimFail << std::endl;
//        std::cout << "Total optim pass: " << optimPass << std::endl;
//        std::cout << "Total Range Skip:" << rangeSkip << std::endl;
//        std::cout << "Total Blos intersect:" << blosIntersectCnt << std::endl;
//        std::cout << "Total BLoS  Clear:" << taskSize-rangeSkip-blosIntersectCnt << std::endl;
        
        //tt = clock() - tt;
        //tt_real = clock() - tt_real;
        //std::cout << "Done CL Proc: " << float(tt)/CLOCKS_PER_SEC << std::endl;
        //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        //std::cout << "Real Secs: " << duration << std::endl;
        
        //Cleanup
//        delete[] h_output;
        
    }
    catch (cl::Error &err) {
        //Get the build log for the first device
        std::cerr << "Building failed, " << err.what() << "(" << err.err() << ")"
        << "\nRetrieving build log\n"
        << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0])
        << "\n";
        return -1;
    }
    
    std::cout << "Done" << std::endl;
    
    return 0;
}









