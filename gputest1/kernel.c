//
//  kernel.c
//  gputest1
//
//  Created by Chris Nicholas on 26/05/2015.
//  Copyright (c) 2015 Chris Nicholas. All rights reserved.
//

#include "kernel.h"

//__kernel int intersect (__global float* x, __global float* y, float a)
//{
    //const int i = get_global_id (0);
    
    //y [i] += a * x [i];


__kernel int intersect (__global float* dem1, __global float* trj1,
                        __global float* dem2, __global float* trj2,
                        int trjSize)
{
    //Dem 1 and Dem 2 are the same but offset by one position to make access work
    // 01234567 Dem 1
    // 12345678 Dem 2
    //So for input list of floats we drop last from first list and drop first from last list = same size
    //This means we can get the appropriate data for each compute
    
    
    //Global ID to input data
    int gID;
    gID = get_global_id(0);
    
    y1 = dem1[gID]
    y2 = dem2[gID]
    y3 = trj1[gID]
    y4 = trj2[gID]
    
    //Check parallel / get delta
    delta = ((x1-x2)*(y3-y4))-((y1-y2)*(x3-x4));
    
    if (delta == 0.0) {
        return 0;
    }
    
    //Calc Sects
    iSectX = ((((x1*y2)-(y1*x2))*(x3-x4))-((x1-x2)*((x3*y4)-(y3*x4))))/delta;
    iSectY = ((((x1*y2)-(y1*x2))*(y3-y4))-((y1-y2)*((x3*y4)-(y3*x4))))/delta;
    
    //Intersect
    if (iSectX >= x1) {
        if (iSectX <= x2) {
            if (iSectY >= std::min(y3, y4)) {
                if (iSectY <= std::max(y3, y4)) {
                    //Intersect
                    return 1;
                }
            }
            
        }
    }
    //Otherwise not
    return 0;
    }
    
}





