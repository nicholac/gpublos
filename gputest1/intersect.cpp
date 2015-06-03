//
//  intersect.cpp
//  gputest1
//
//  Created by Chris Nicholas on 01/06/2015.
//  Copyright (c) 2015 Chris Nicholas. All rights reserved.
//

#include "intersect.h"



int intersectUtils::intersect(const std::vector<std::shared_ptr<position>>& inTrj,
                              const std::vector<std::shared_ptr<position>>& inWorld, const double& linLen)  {
    /*
     Takes trajectory vector and world dem profile and compares for intersect
     TODO: Check to make sure the vector lengths are the same length
     This does the intersecting  for two arrays containing Z values.
     Reduces problem to 2d as the lines are coplanar
     Input:
     WorldZ: array of dem z values
     trjZ: Array of trajectory z values
     linLen: float of line length (hypotenuse of x, y)
     */
    double x1 = 0.0;
    double y1;
    double x2 = linLen;
    double y2;
    double x3 = 0.0;
    double y3;
    double x4 = linLen;
    double y4;
    double delta = 0.0;
    double iSectX = 0.0;
    double iSectY = 0.0;
    
    //TODO: Clean this up as it slows things down
    int trjSize = static_cast<int>(inTrj.size());
    int worldSize = static_cast<int>(inWorld.size());
    
    //Compare vector sizes for safety
    if (trjSize != worldSize) {
        //TODO: Do something sensible here
        printf ("warning: traj and dem different length - intersect will fail\n");
    }
    //TODO: change all project .at to [] for speed?
    for (int i=0; i < (trjSize-1); i++) {
        //std::cout << i << std::endl;
        //Get the z values from vectors
        y1 = inWorld.at(i)->z;
        y2 = inWorld.at(i+1)->z;
        y3 = inTrj.at(i)->z;
        y4 = inTrj.at(i+1)->z;
        
        //check for parallel lines
        intersectUtils::chkParallel(x1, x2, x3, x4, y1, y2, y3, y4, delta);
        
        if (delta == 0.0) {
            continue;
        }
        
        intersectUtils::calcSects(x1, x2, x3, x4, y1, y2, y3, y4, delta, iSectX, iSectY);
        //Now we just check if this point is within the bounds of our line
        bool iSect = intersectUtils::checkIntersect(x1, x2, x3, x4, y1, y2, y3, y4, iSectX, iSectY);
        if ( iSect == true ) {
            return 2;
        }
        
    }// end for
    //No intersect at end of loop if not broken
    return 1;
}