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
#include <recShot.h>


int main(int argc, const char * argv[]) {

    char* kernelsSrc = "/Users/dusted-dstl/Documents/xcodeworkspace/gpublos/gputest1/lut.cl";
    char* demFName = "/Users/dusted-dstl/Documents/geodata/mount_chip.tif";
    //char* demFName = "/Users/dusted-dstl/Documents/geodata/mount.dem";
    char* outfName = "/Users/dusted-dstl/Documents/geodata/gpu_soup3.tif";
    //Linux
    //char* kernelsSrc = "/home/ec2-user/gputest1/simplekern.cl";
    float tgtX = 559770;
    float tgtY = 5119808.0;
    float tgtZ = 3856.0;
    //Dummy target position
    cl_float4 tgtPos = {tgtX, tgtY, tgtZ, 0.0};
    
    //Dummy find position
    float laX = 100.0;
    float laY = 100.0;
    float laZ = 100.0;
    float dist = sqrt(pow(laX, 2.0)+pow(laY, 2.0));
    cl_float4 fiPos = {dist, laZ, 45.0, 0.0};
    
    int trjArrSize = 3000;
    int findIdx = 2500;
    
    float minElev = 45.0;
    float maxElev = 85.0;
    float thetaStep = 0.05;
    
    float timeStep = 0.05;
    
    float airDens = 1.225;
    float dragCoef = 0.15; //Coef for sphere
    float calibre = 0.081;
    float fArea = pi*(pow((0.5*calibre), 2.0));
    float mortSigma = dragCoef*fArea*0.5*airDens;
    float muzzVel = 225.0;
    float mortMass = 3.2;
    //Maximum error acceptable between trj and tgt (m)
    float maxErr = 20.0;
    
    
    
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
        std::cout << "Compiling Kernels..." << std::endl;
        program.build(devices);
        std::cout << "Done Compiling" << std::endl;
        
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
    cl_float8 cl_geoTrans = {float(adfGeoTransform[0]), float(adfGeoTransform[1]), float(adfGeoTransform[2]), float(adfGeoTransform[3]), float(adfGeoTransform[4]), float(adfGeoTransform[5])};
    
    int rasterXSize = poDataset->GetRasterXSize();
    int rasterYSize = poDataset->GetRasterYSize();
    cl_float2 cl_rasterSize = {float(rasterXSize), float(rasterYSize)};
    
    int taskSize = rasterXSize*rasterYSize;
    std::cout << "Done inData" << std::endl;

    std::cout << "Building trajectories..." << std::endl;
    //Build the trajectory data
    //Compute Air Densities
    worldParams::computeAirDensLUT();
    //tgtZ is a baseline height for computation of air densities
    std::vector<cl_float4> h_trjArrVec = Utils::genTrjSoup(thetaStep, minElev, maxElev,
                                                        timeStep, muzzVel, mortSigma, mortMass, 0.0);
    std::cout << "Done trajectories" << std::endl;
    
    //Get size of the trjArr
    trjArrSize = h_trjArrVec.size();
    //Convert to pointer array (not std vec)
    cl_float4 *h_trjArr = new cl_float4[trjArrSize];
    for (int i=0; i<trjArrSize; i++){
        h_trjArr[i]=h_trjArrVec.at(i);
    }
    
    //Build the tgt data - (x,y,z,theta,dist2tgt,S,S,S)
    cl_float8 *h_tgtArr = new cl_float8[taskSize];
    Utils::GDAL2FLOAT8 (poDataset, h_tgtArr);
    
    //Build 1d Dem Z array
    cl_float *h_demZArr = new cl_float[taskSize];
    for (int i=0; i<taskSize; i++){
        h_demZArr[i]=h_tgtArr[i].s2;
    }
    
    //Do some tests on the data
    double testX = 561010.0;
    double testY = 5120024.0;
    double testZ = 4070.0;
    int pxX, pxY;
    Utils::coordtoPx2d(testX, testY, pxX, pxY, adfGeoTransform, rasterXSize, rasterYSize);
    int zIdx = (pxY*rasterXSize)+pxX;
    double chkX, chkY;
    Utils::px2Coord(chkX, chkY, pxX, pxY, adfGeoTransform, rasterXSize, rasterYSize);
    std::cout << std::setprecision(10) << std::endl;
    std::cout << std::setprecision(10) << "GeoTrans(0-6): " << adfGeoTransform[0] << ", " << adfGeoTransform[1] << ", " << adfGeoTransform[2] << ", " << adfGeoTransform[3] << ", " << adfGeoTransform[4] << ", " << adfGeoTransform[5] << ", " << std::endl;
    std::cout << std::setprecision(10) << "Rastersize: " << rasterXSize << ", " << rasterYSize << std::endl;
    std::cout << std::setprecision(10) << "TestCoord X, Y: " << testX<< ", " << testY << std::endl;
    std::cout << std::setprecision(10) << "Converted Pixel X, Y: " <<  pxX << " " << pxY << std::endl;
    std::cout << std::setprecision(10) << "Converted 1dIdx: " <<  zIdx << std::endl;
    std::cout << std::setprecision(10) << "TgtArr X, Y at Position: " <<  chkX << " " << chkY << std::endl;
    std::cout << std::setprecision(10) << "Test Height: " <<  testZ << std::endl;
    std::cout << std::setprecision(10) << "TgtArrZ at position: " <<  h_tgtArr[zIdx].s2 << std::endl;
    std::cout << std::setprecision(10) << "demZArr at position: " <<  h_demZArr[zIdx] << std::endl;
    std::cout << std::setprecision(10) << std::endl;
    
    
    //Build buffers and write to device
    try {
        std::cout << "Loading Buffers " << std::endl;
        cl::Buffer d_trjArr = cl::Buffer(context, CL_MEM_READ_ONLY, trjArrSize*sizeof(cl_float4), NULL, &err);
        cl::Buffer d_tgtArr = cl::Buffer(context, CL_MEM_READ_WRITE, taskSize*sizeof(cl_float8), NULL, &err);
        cl::Buffer d_demZArr = cl::Buffer(context, CL_MEM_READ_ONLY, taskSize*sizeof(cl_float), NULL, &err);
        
        std::cout << "Writing Buffers " << std::endl;
        //Write Buffers to device
        queue.enqueueWriteBuffer(d_trjArr,CL_TRUE,0,trjArrSize*sizeof(cl_float4), h_trjArr);
        queue.enqueueWriteBuffer(d_tgtArr,CL_TRUE,0,taskSize*sizeof(cl_float8), h_tgtArr);
        queue.enqueueWriteBuffer(d_demZArr,CL_TRUE,0,taskSize*sizeof(cl_float), h_demZArr);
        queue.finish();
        std::cout << "Done Writing Buffers " << std::endl;
        
        
        //Run the distance creation kernel
        std::cout << "Building kernel" << std::endl;
        cl::Kernel genDists = cl::Kernel(program, "genDists", &err);
        cl::Kernel trjRange2 = cl::Kernel(program, "trjRange2", &err);
        cl::Kernel trjIntersect = cl::Kernel(program, "trjIntersect", &err);
        
        
        std::cout << "Enqueue Args - genDists " << std::endl;
        err = genDists.setArg(0, d_tgtArr);
        err = genDists.setArg(1, tgtPos);
        queue.finish();
        

        std::cout << "Enqueue Args - trjRange " << std::endl;
        err = trjRange2.setArg(0, d_trjArr);
        err = trjRange2.setArg(1, d_tgtArr);
        err = trjRange2.setArg(2, trjArrSize);
        err = trjRange2.setArg(3, tgtPos);
        queue.finish();
        
        
        std::cout << "Enqueue Args - trjIntersect " << std::endl;
        err = trjIntersect.setArg(0, d_trjArr);
        err = trjIntersect.setArg(1, d_tgtArr);
        err = trjIntersect.setArg(2, d_demZArr);
        err = trjIntersect.setArg(3, trjArrSize);
        err = trjIntersect.setArg(4, cl_geoTrans);
        err = trjIntersect.setArg(5, cl_rasterSize);
        err = trjIntersect.setArg(6, tgtPos);
        err = trjIntersect.setArg(7, maxErr);
        queue.finish();
        
        //Start CL processing
        clock_t tt = clock();
        
        std::cout << "Running " << taskSize << " Targets..." << std::endl;
        
        std::cout << "Distance..." << std::endl;
        //Generate the distances data from each DEM postion --> tgt
        err = queue.enqueueNDRangeKernel(genDists, cl::NullRange, cl::NDRange(taskSize), cl::NullRange, NULL, NULL);
        queue.finish();
        
        std::cout << "Trj..." << std::endl;
        //Optimise using trj soup
        err = queue.enqueueNDRangeKernel(trjRange2, cl::NullRange, cl::NDRange(taskSize), cl::NullRange, NULL, NULL);
        queue.finish();
        
        //Intersect
        std::cout << "Interescting..." << std::endl;
        err = queue.enqueueNDRangeKernel(trjIntersect, cl::NullRange, cl::NDRange(taskSize), cl::NullRange, NULL, NULL);
        queue.finish();

        //Read results
        err = queue.enqueueReadBuffer(d_tgtArr, CL_TRUE, 0, taskSize*sizeof(cl_float8), h_tgtArr, NULL, NULL);
        queue.finish();
        
        //Check results
        int intCnt = 0;
        int outSideErr = 0;
        for (int i=0; i<taskSize; i++){
            if (h_tgtArr[i].s6 > 0.0){
                intCnt++;
            }
            if (h_tgtArr[i].s5 > maxErr){
                outSideErr++;
            }
        }
        std::cout << "Total tgt Intersects: " << intCnt << std::endl;
        std::cout << "Total outside err: " << outSideErr << std::endl;
        std::cout << "Total good shots: " << taskSize-intCnt-outSideErr << std::endl;
        
        //Write to raster
        std::cout << "Dumping Raster" << std::endl;
        int chk = Utils::writeRasterOut (poDataset, outfName, h_tgtArr);
        //Close base DEM
        GDALClose( (GDALDatasetH) poDataset );
        std::cout << "Write raster result: " << chk << std::endl;
        
        //Cleanup
        delete[] h_tgtArr;
        delete[] h_trjArr;
        
        tt = clock() - tt;
        std::cout << "Done CL Proc: " << float(tt)/CLOCKS_PER_SEC << std::endl;
        
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

