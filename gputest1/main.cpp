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
#include <Utils.h>
//#include <worldParams.h>
//#include <motionUtils.h>
#include <iomanip>
//#include <recShot.h>


int main(int argc, const char * argv[]) {
    
    //Files, vars etc
    char* demFName = "/Users/dusted-dstl/Documents/geodata/mount_chip.tif";
    //char* demFName = "/Users/dusted-dstl/Documents/geodata/mount.dem";
    char* kernelsSrc = "/Users/dusted-dstl/Documents/xcodeworkspace/gpublos/gputest1/simplekern.cl";
    float tgtX = 559786.0;
    float tgtY = 5120048.0;
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
    for (auto plat : platforms) {
        std::cout << plat.getInfo<CL_PLATFORM_NAME>() << std::endl;
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
    for (auto dev : devices) {
        std::cout << "Device Name: " << dev.getInfo<CL_DEVICE_NAME>() << std::endl;
        std::cout << "Device Global Mem max size (mb): " << dev.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>()/1000000 << std::endl;
        std::cout << "Device Local Mem max size: " << dev.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << std::endl;
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
    std::cout << gt0 << " // " << adfGeoTransform[0] << std::endl;
    int rasterXSize = poDataset->GetRasterXSize();
    int rasterYSize = poDataset->GetRasterYSize();
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
    
    GDALClose( (GDALDatasetH) poDataset );
    std::cout << "Size of DEM in memory Ref (kb): " << (taskSize*sizeof(float))/1000 << std::endl;
    std::cout << "Size of DEM: " << taskSize << " tasksize: " << taskSize << std::endl;
    std::cout << "Size of DEM as Float8 (mb): " << taskSize*sizeof(cl_float8)/1000000 << std::endl;
    std::cout << "Will the whole DEM fit into Global mem @ float8?: " << std::endl;
    bool chkGlobalMem = 0;
    for (auto dev : devices) {
        chkGlobalMem = taskSize*sizeof(cl_float8) < dev.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();
        std::cout << "Device Name: " << dev.getInfo<CL_DEVICE_NAME>() << std::endl;
        std::cout << "Sizes (dem, global mem) & Chk: " << taskSize*sizeof(cl_float8) << " // " << dev.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << " // " << chkGlobalMem << std::endl;
    }
    if (chkGlobalMem == 0){
        std::cout << "Not enough Global memory for the DEM on this device - aborting " << std::endl;
        return 1;
    }
    
    //Build the array of possible launch angles
    std::vector<float> thetaArr;
    for (float i=minElev; i<maxElev; i+=thetaStep){
        thetaArr.push_back(i);
    }
    float *h_thetaArr = new float[thetaArr.size()];
    for (int i=0; i<thetaArr.size(); i++){
        h_thetaArr[i] = thetaArr.at(i);
    }
    
    //Build the array of output distances for each theta
    float *h_distArr = new float[thetaArr.size()];
    std::cout << "Parrellel Shots per DEM Px: " << thetaArr.size()<< std::endl;
    std::cout << "Total shots: " << thetaArr.size()*taskSize<< std::endl;
    
    clock_t tt = clock();
    clock_t tt_real = clock();
    
    try {
        std::cout << "Loading Buffers " << std::endl;
        //THESE TWO FOR CPU - then COMMENT CPU
        //Create a buffer object to this memory
        cl::Buffer d_inBuffDEMFull = cl::Buffer(context, CL_MEM_READ_WRITE, taskSize*sizeof(cl_float8), NULL, &err);
        cl::Buffer d_outBuff = cl::Buffer(context, CL_MEM_READ_WRITE, taskSize*sizeof(float), NULL, &err);
        cl::Buffer d_thetaArr = cl::Buffer(context, CL_MEM_READ_ONLY, thetaArr.size()*sizeof(float), NULL, &err);
        cl::Buffer d_distArr = cl::Buffer(context, CL_MEM_READ_WRITE, thetaArr.size()*sizeof(float), NULL, &err);
        
        //Out data on the host
        float *h_output = new float[taskSize];
        
        std::cout << "Writing Buffers " << std::endl;
        //Write Buffers to device
        //Theta Arr
        queue.enqueueWriteBuffer(d_thetaArr,CL_TRUE,0,thetaArr.size()*sizeof(float), h_thetaArr);
        queue.enqueueWriteBuffer(d_inBuffDEMFull,CL_TRUE,0,taskSize*sizeof(cl_float8), h_demSP);
        
        //WORKS
        //Build Kernels
        //cl::KernelFunctor rangeChk(cl::Kernel(program,"rangeChk"),queue,cl::NullRange,cl::NDRange(taskSize),cl::NullRange);
        //cl::KernelFunctor shotOpt(cl::Kernel(program,"shotOpt"),queue,cl::NullRange,cl::NDRange(thetaArr.size()),cl::NullRange);
        //cl::KernelFunctor thetaDistReduce(cl::Kernel(program,"thetaDistReduce"),queue,cl::NullRange,cl::NDRange(thetaArr.size()),cl::NullRange);
        
        std::cout << "Building kernels " << std::endl;
        cl::Kernel rangeChk = cl::Kernel(program, "rangeChk", &err);
        cl::Kernel shotOpt = cl::Kernel(program, "shotOpt", &err);
        cl::Kernel thetaDistReduce = cl::Kernel(program, "thetaDistReduce", &err);
        
        //Dummy target position
        cl_float4 tgtPos = {559786.0, 5120048.0, 0.0, 0.0};
        
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
        
        std::cout << "Enqueue Args - shotOpt " << std::endl;
        //Build the static input args - ShotOpt
        err = shotOpt.setArg(0, d_thetaArr);
        err = shotOpt.setArg(1, d_distArr);
        err = shotOpt.setArg(2, tgtPos);
        err = shotOpt.setArg(4, timeStep);
        queue.finish();
        
        std::cout << "Enqueue Args - thetaReduce" << std::endl;
        //Build the static input args - Dist Reduce
        err = thetaDistReduce.setArg(0, d_thetaArr);
        err = thetaDistReduce.setArg(1, d_distArr);
        err = thetaDistReduce.setArg(2, d_inBuffDEMFull);
        queue.finish();
        
        
        //CL - Loop the optimisation brute force
        std::cout << "Starting Loop" << std::endl;
        for (int i=0; i<1000; i++){
            err = shotOpt.setArg(3, h_demSP[i]);
            err = queue.enqueueNDRangeKernel(shotOpt, cl::NullRange, cl::NDRange(thetaArr.size()), cl::NullRange, NULL, NULL);
            queue.finish();
            
            err = thetaDistReduce.setArg(3, i);
            err = queue.enqueueNDRangeKernel(thetaDistReduce, cl::NullRange, cl::NDRange(thetaArr.size()), cl::NullRange, NULL, NULL);
            queue.finish();
            
            std::cout << "Done: " << i << std::endl;
            
            //WORKS but Slow - kernel Functors
            //Fire Shots
            //Temporary tgt pos for testing
            //float8(x,y,z,theta,dist,runflag,blos,tlos)
            //cl_float8 pos = {569786.0, 5130048.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
            //shotOpt(d_thetaArr, d_distArr, tgtPos, pos, timeStep);
            //shotOpt(d_thetaArr, d_distArr, tgtPos, h_demSP[i], timeStep);
            //Write the closest into the DEM arr
            //thetaDistReduce(d_thetaArr, d_distArr, d_inBuffDEMFull, i);
        }
        
        //Read the processed DEM data back out
        err = queue.enqueueReadBuffer(d_inBuffDEMFull, CL_TRUE, 0, taskSize*sizeof(cl_float8), h_demSP, NULL, NULL);
        // wait for completion
        queue.finish();
        
        std::cout << "Example output value: " << h_demSP[5].s4 << std::endl;
        //Check how many no run flags we have
        int cnt = 0;
        for (int i=0; i<taskSize; i++){
            if (h_demSP[i].s5 == 0.0) {
                cnt++;
            }
        }
        std::cout << "Total no run flags // Total tasksize: " << cnt << " // " << taskSize << std::endl;
        
        tt = clock() - tt;
        tt_real = clock() - tt_real;
        std::cout << "Done CL Proc: " << float(tt)/CLOCKS_PER_SEC << std::endl;
        std::cout << "Real Clock: " << float(tt_real) << std::endl;
        
        //Cleanup
        delete[] h_output;
        
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









