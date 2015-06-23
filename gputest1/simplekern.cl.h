//
//  simplekern.h
//  gputest1
//
//  Created by Chris Nicholas on 10/06/2015.
//  Copyright (c) 2015 Chris Nicholas. All rights reserved.
//

#ifndef gputest1_simplekern_h
#define gputest1_simplekern_h

// This include pulls in everything you need to develop with OpenCL in OS X.
#include <OpenCL/opencl.h>

extern void (^intersect)(cl_float* input, cl_float* group_result);


#endif
