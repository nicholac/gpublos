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
#include <cl.hpp>
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
    char* demFName = "/Users/dusted-ipro/geodata/mount.dem";
    char* intersectKernel = "/Users/dusted-ipro/Documents/code/code/xcode/gputest1/gputest1/simplekern.cl";
    float tgtX = 562683.0;
    float tgtY = 5116767.0;
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
    float tgtSize = 50.0;
    char* fName = "/Users/dusted-ipro/geodata/gpu_out1.csv";
    
    
    
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
        std::cout << dev.getInfo<CL_DEVICE_NAME>() << std::endl;
    }
    
    
    // create platform
    cl::Platform::get(&platforms);
    platforms[0].getDevices(CL_DEVICE_TYPE_GPU, &devices);
    //Query device info for memory size
    std::cout << devices[0].getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << std::endl;
    
    //Change this if want to use GPU
    cl::Context context = cl::Context(devices);
    
    //Create Command queue
    cl::CommandQueue queue(context, devices[0]);
    
    // load opencl source
    std::ifstream cl_file(intersectKernel);
    std::string cl_string(std::istreambuf_iterator<char>(cl_file), (std::istreambuf_iterator<char>()));
    cl::Program::Sources source(1, std::make_pair(cl_string.c_str(), cl_string.length()+1));
    
    // create program
    cl::Program program(context, source);
    
    try {
        
        // compile opencl source
        program.build(devices);
        std::cout << "built kernel" << std::endl;
        
    }
    catch (cl::Error &err) {
        //Get the build log for the first device
        std::cerr << "Building failed, " << err.what() << "(" << err.err() << ")"
        << "\nRetrieving build log\n"
        << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0])
        << "\n";
        return -1;
        
    }
    
    std::cout << "Loading Kern func" << std::endl;
    // load named kernel from opencl source
    cl::Kernel kernel(program, "rangeChk");
    
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
    float rasterXSize = float(poDataset->GetRasterXSize());
    float rasterYSize = float(poDataset->GetRasterYSize());
    int taskSize = rasterXSize*rasterYSize;
    std::cout << std::fixed << "rasterSizeX: " << rasterXSize << std::endl;
    std::cout << std::fixed << "rasterSizeY: " << rasterYSize << std::endl;
    
    GDALRasterBand  *poBand;
    poBand = poDataset->GetRasterBand( 1 );
    int   nXSize = poBand->GetXSize();
    int   nYSize = poBand->GetYSize();
    
    //This will give us a float4 with x, y, z (x and y in world coords)
    std::vector<cl_float4> demSP = Utils::GDAL2VEC2D (poDataset);
    GDALClose( (GDALDatasetH) poDataset );
    std::cout << "Size of DEM in memory Ref (kb): " << (taskSize*sizeof(float))/1000 << std::endl;
    std::cout << "Size of DEM: " << demSP.size() << " tasksize: " << taskSize << std::endl;
    //Get the tgtz ready
    Utils::coordtoPx2d(tgtX, tgtY, tgtPxX, tgtPxY, adfGeoTransform, rasterXSize, rasterYSize);
    //tgtZ = demSP.at(tgtPxY).at(tgtPxX);
    tgtZ = demSP.at(tgtPxY*rasterXSize+tgtPxX).z;
    std::cout << " TgtZ: " << tgtZ << std::endl;
    
    clock_t tt = clock();
    
//RANGE//
    try {
        //Input buffer object - here we are setting useHostPtr to true, so the data is not being copied across to the device
        //Instead it is passed the pointer and data stays on host
        cl::Buffer inBuffDEMFull(context, demSP.begin(), demSP.end(), true, true);
        //Here the data is copied across to the device
        //cl::Buffer inBuffDEMFull(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*demSP.size(), demSP, &err);
        
        
        //Output buffer - tell device to allocate this memory on device THIS FOR CPU AND GPU
        cl::Buffer outDistBuff(context, CL_MEM_WRITE_ONLY, taskSize*sizeof(cl_float), NULL, &err);

        //Set Kernel Args - range
        err = kernel.setArg(0, inBuffDEMFull);
        err = kernel.setArg(1, sizeof(float), &tgtX);
        err = kernel.setArg(2, sizeof(float), &tgtY);
        err = kernel.setArg(3, sizeof(float), &minDist);
        err = kernel.setArg(4, sizeof(float), &maxDist);
        err = kernel.setArg(5, outDistBuff);
        if(err < 0) {
            perror("Couldn't create a kernel argument");
            exit(1);
        }
        // execute kernel
        //This execute just one work item - queue.enqueueTask(kernel);
        cl::NDRange offset(0);
        cl::NDRange local_size(512);
        cl::NDRange global_size(taskSize);
        queue.enqueueNDRangeKernel(kernel, offset, global_size);
        queue.finish();
        
        //Write the outputs
        //GET DATA - THIS FOR GPU - then comment next
        //cl::copy(queue, outBuff, h_outArr.begin(), h_outArr.end());
        
        //GET DATA
        std::vector<cl_float> h_outArr = std::vector<cl_float>(taskSize);
        //Copy the data back out
        cl::copy(queue, outDistBuff, h_outArr.begin(), h_outArr.end());
        
        //Reduce by distance
        
        
        for (int i=0; i<taskSize; i++) {
            std::cout << std::fixed << h_outArr[i] << std::endl;
        }

        
////OPTIMISE//
//        cl::Kernel shotKernel(program, "shot");
//        //Load optimise kernel
//        cl::Kernel optimise(program, "optimise");
//        //Load miss dist kernel
//        cl::Kernel missDist(program, "missDist");
//    
//        //In DEM Buffer (2D)
//        cl::Buffer inBuffDEMFull(context, demSP.begin(), demSP.end(), true, true);
//        
//        //Output Buffer of distances
//        cl::Buffer d_missPosBuff(context, CL_MEM_READ_WRITE, taskSize*sizeof(cl_float2), NULL, &err);
//        
//        //Opt El Buffer
//        cl::Buffer d_optElBuff(context, CL_MEM_READ_WRITE, taskSize*sizeof(cl_float), NULL, &err);
//        
//        //Set Kernel Args - Single Shot
//        err = shotKernel.setArg(0, inBuffDEMFull);
//        err = shotKernel.setArg(1, sizeof(float), &tgtX);
//        err = shotKernel.setArg(2, sizeof(float), &tgtY);
//        err = shotKernel.setArg(3, sizeof(float), &rasterXSize);
//        err = shotKernel.setArg(4, sizeof(float), &rasterYSize);
//        err = shotKernel.setArg(5, sizeof(float), &gt0);
//        err = shotKernel.setArg(6, sizeof(float), &gt1);
//        err = shotKernel.setArg(7, sizeof(float), &gt2);
//        err = shotKernel.setArg(8, sizeof(float), &gt3);
//        err = shotKernel.setArg(9, sizeof(float), &gt4);
//        err = shotKernel.setArg(10, sizeof(float), &gt5);
//        err = shotKernel.setArg(11, sizeof(float), &mortMass);
//        err = shotKernel.setArg(12, sizeof(float), &muzzVel);
//        err = shotKernel.setArg(13, sizeof(float), &mortSigma);
//        err = shotKernel.setArg(14, sizeof(float), &gForce);
//        err = shotKernel.setArg(15, sizeof(float), &tgtZ);
//        err = shotKernel.setArg(16, d_missPosBuff);
//        if(err < 0) {
//            perror("Couldn't create a kernel argument for shot kern");
//            exit(1);
//        }
//        
//        //Set Kernel Args - Miss Dist
//        err = missDist.setArg(0, d_missPosBuff);
//        err = missDist.setArg(1, sizeof(float), &tgtX);
//        err = missDist.setArg(2, sizeof(float), &tgtY);
//        err = missDist.setArg(3, sizeof(float), &tgtSize);
//        if(err < 0) {
//            perror("Couldn't create a kernel argument for miss dist");
//            exit(1);
//        }
//        
//        //Optimise loop - shots from 85.0 - 45.0 at 0.05 intervals
//        for (float i=45.0; i<85.0; i+=20.0) {
//            
//            //Fire Shot (cl)
//            
//            //Get range from land to tgt (cl)
//            
//            //Copy distance data out & check against previous (master)
//        
//        }
    
            
////OPTIMISE//
//
//        std::cout << "Loading Kern funcs" << std::endl;
//        // load setup opt el
//        cl::Kernel setOptElKern(program, "setupOptEl");
//        // load setup seedjump
//        cl::Kernel setupSeedJump(program, "setupSeedJump");
//        // load shot kernel
//        cl::Kernel shotKernel(program, "shot");
//        //Load optimise kernel
//        cl::Kernel optimise(program, "optimise");
//        //Load miss dist kernel
//        cl::Kernel missDist(program, "missDist");
//        //Load Chk otp kernel
//        cl::Kernel checkOpt(program, "checkOpt");
//
//        
//        //Create new input buffer for procIdx (output from previous)
//        //Again we are not copying this into the device
//        cl::Buffer d_inBuffProxIdx(context, h_outArr.begin(), h_outArr.end(), true, true);
//        
//        //Opt El Buffer
//        cl::Buffer d_optElBuff(context, CL_MEM_READ_WRITE, taskSize*sizeof(cl_float), NULL, &err);
// 
//        //Output Buffer of distances
//        cl::Buffer d_missPosBuff(context, CL_MEM_READ_WRITE, taskSize*sizeof(cl_float2), NULL, &err);
//        
//        //Stored buffer of previous dotP's
//        cl::Buffer d_dotPBuff(context, CL_MEM_READ_WRITE, taskSize*sizeof(cl_float), NULL, &err);
//        
//        //Stored buffer of seedJumps
//        cl::Buffer d_seedJumpBuff(context, CL_MEM_READ_WRITE, taskSize*sizeof(cl_float), NULL, &err);
//        
//        //Set Kernel Args - Setup Opt El
//        err = setOptElKern.setArg(0, d_optElBuff);
//        
//        //Set Kernel Args - Setup Seed Jump
//        err = setupSeedJump.setArg(0, d_seedJumpBuff);
//        
//        //Set them up
//        queue.enqueueNDRangeKernel(setOptElKern, offset, global_size);
//        queue.enqueueNDRangeKernel(setupSeedJump, offset, global_size);
//        queue.finish();
//
//        //Set Kernel Args - Single Shot
//        err = shotKernel.setArg(0, inBuffDEMFull);
//        err = shotKernel.setArg(1, sizeof(float), &tgtX);
//        err = shotKernel.setArg(2, sizeof(float), &tgtY);
//        err = shotKernel.setArg(3, sizeof(float), &rasterXSize);
//        err = shotKernel.setArg(4, sizeof(float), &rasterYSize);
//        err = shotKernel.setArg(5, sizeof(float), &gt0);
//        err = shotKernel.setArg(6, sizeof(float), &gt1);
//        err = shotKernel.setArg(7, sizeof(float), &gt2);
//        err = shotKernel.setArg(8, sizeof(float), &gt3);
//        err = shotKernel.setArg(9, sizeof(float), &gt4);
//        err = shotKernel.setArg(10, sizeof(float), &gt5);
//        err = shotKernel.setArg(11, sizeof(float), &mortMass);
//        err = shotKernel.setArg(12, sizeof(float), &muzzVel);
//        err = shotKernel.setArg(13, sizeof(float), &mortSigma);
//        err = shotKernel.setArg(14, sizeof(float), &gForce);
//        err = shotKernel.setArg(15, sizeof(float), &tgtZ);
//        err = shotKernel.setArg(16, d_inBuffProxIdx);
//        err = shotKernel.setArg(17, d_missPosBuff);
//        if(err < 0) {
//            perror("Couldn't create a kernel argument for shot kern");
//            exit(1);
//        }
//        
//        //Optimise Args
//        err = optimise.setArg(0, sizeof(float), &rasterXSize);
//        err = optimise.setArg(1, sizeof(float), &rasterYSize);
//        err = optimise.setArg(2, sizeof(float), &gt0);
//        err = optimise.setArg(3, sizeof(float), &gt1);
//        err = optimise.setArg(4, sizeof(float), &gt2);
//        err = optimise.setArg(5, sizeof(float), &gt3);
//        err = optimise.setArg(6, sizeof(float), &gt4);
//        err = optimise.setArg(7, sizeof(float), &gt5);
//        err = optimise.setArg(8, sizeof(float), &tgtX);
//        err = optimise.setArg(9, sizeof(float), &tgtY);
//        err = optimise.setArg(10, d_missPosBuff);
//        err = optimise.setArg(11, d_dotPBuff);
//        err = optimise.setArg(12, d_seedJumpBuff);
//        err = optimise.setArg(13, d_optElBuff);
//        err = optimise.setArg(14, d_inBuffProxIdx);
//        if(err < 0) {
//            perror("Couldn't create a kernel argument for optimise kern");
//            exit(1);
//        }
//        
//        //Set Kernel Args - Miss Dist
//        err = missDist.setArg(0, d_missPosBuff);
//        err = missDist.setArg(1, d_inBuffProxIdx);
//        err = missDist.setArg(2, sizeof(float), &tgtX);
//        err = missDist.setArg(3, sizeof(float), &tgtY);
//        err = missDist.setArg(4, sizeof(float), &tgtSize);
//        if(err < 0) {
//            perror("Couldn't create a kernel argument for miss dist");
//            exit(1);
//        }
//        
//        //Set Kernel Args - Check Opt
//        err = checkOpt.setArg(0, d_optElBuff);
//        err = checkOpt.setArg(1, d_seedJumpBuff);
//        err = checkOpt.setArg(2, d_inBuffProxIdx);
//        if(err < 0) {
//            perror("Couldn't create a kernel argument for Check Opt");
//            exit(1);
//        }
//
//        
//        std::cout << "Queuing Executing Kernels" << std::endl;
//        // Setup Opt El kernel
//        queue.enqueueNDRangeKernel(setOptElKern, offset, global_size);
//        
//        // execute first shot
//        queue.enqueueNDRangeKernel(shotKernel, offset, global_size);
//        queue.finish();
//        
//        std::cout << "Done First Run" << std::endl;
//        
//        //For loop starts here //
//        for (int i=0; i < 2; i++){
//            std::cout << "Done Opt Loop: " << i << std::endl;
//        
//            //Distance Chk
//            queue.enqueueNDRangeKernel(missDist, offset, global_size);
//            
//            //Optimise
//            queue.enqueueNDRangeKernel(optimise, offset, global_size);
//            
//            //Check Status
//            queue.enqueueNDRangeKernel(checkOpt, offset, global_size);
//            
//            //execute shot
//            queue.enqueueNDRangeKernel(shotKernel, offset, global_size);
//            queue.finish();
//            
//            //Lets see how many proc's we're left with
//            std::vector<cl_int> h_chkProcIdx = std::vector<cl_int>(taskSize);
//            cl::copy(queue, d_inBuffProxIdx, h_chkProcIdx.begin(), h_chkProcIdx.end());
//            int totProc;
//            for (int p : h_chkProcIdx){
//                totProc+=p;
//            }
//            std::cout << "Tot Proc Idx: " << totProc << std::endl;
//        }
//        
//        //End loop here - we now have an optelArr with optimum elevations for shots and procIdx with 0 and 1 for processing
//        
//        //std::cout << "Queuing Copy Output" << std::endl;
//        //Copy the output
//        std::vector<cl_float> h_outOptArr = std::vector<cl_float>(taskSize);
//        cl::copy(queue, d_missPosBuff, h_outOptArr.begin(), h_outOptArr.end());
//        
//        // wait for completion
//        std::cout << "Waiting for completion - Optimise" << std::endl;
        

        //Check some outputs
//        for (auto p : h_outOptArr ){
//                std::cout.precision(15);
//                std::cout << std::fixed << p << std::endl;
//        }
        
////        for (float p : d_optElArr){
////            if (p != 8888.0f && p > 0.0f){
////                std::cout.precision(15);
////                std::cout << " OptEl Success: " << std::fixed << p << std::endl;
////            }
////        }
//        std::cout.precision(15);
//        std::cout << " minDist: " << std::fixed << minDist << std::endl;
        
        
        
        tt = clock() - tt;
        std::cout << "Done CL Proc: " << float(tt)/CLOCKS_PER_SEC << std::endl;
        
        //THIS FOR GPU
        //        int tst = 0;
        //        for (int i; i<taskSize; i++){
        //            if (h_outArr.at(i) > tst){
        //                tst=h_outArr.at(i);
        //                std::cout << "Max Val Out: " << tst << std::endl;
        //            }
        //        }
        
        
        //Find Maximum in outputs - checking
        //int tst = 0;
        //for (int i; i<taskSize; i++){
        //    if (outArr[i] > tst){
        //        tst=outArr[i];
        //        std::cout << "Max Val Out: " << tst << std::endl;
        //    }
        //}
        //std::cout << "Max Val Out: " << tst << std::endl;
        
        /*
         //Drop the output to a textfile
         std::ofstream output(fName);
         if (!output) { // check the file opened OK
         std::cerr << "error: open file for output failed!" << std::endl;
         return 1;
         }
         //Headers
         output << "Position X" << "," << "Position Y" << "," << "pixelX" << "," << "pixelY" << "," << "result" << std::endl;
         //Sizes should all be same
         for (unsigned int i=0; i<taskSize; i++){
         float pxX, pxY;
         //Work out pixel location adn launch from 1d index
         float tmp = std::modf(float(i)/rasterXSize, &pxY);
         float tmp2 = std::modf(tmp * rasterXSize, &pxX);
         //std::cout << gt0 << " // " << gt1 << " // " << gt2 << std::endl;
         float X = gt0+(pxX*gt1)+pxY*gt2;
         float Y = gt3+(pxX*gt4)+pxY*gt5;
         // the ofstream object replaces std::cout here, then its the same use
         output << std::setprecision(10) << X << "," << Y << "," << pxX << "," << pxY << "," << h_outArr[i] << std::endl;
         }
         output.close();
         */
        
        
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
