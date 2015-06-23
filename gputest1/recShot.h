//
//  recShot.h
//  gputest1
//
//  Created by Chris Nicholas on 11/06/2015.
//  Copyright (c) 2015 Chris Nicholas. All rights reserved.
//

#ifndef __gputest1__recShot__
#define __gputest1__recShot__

#include <stdio.h>

#include <vector>
#include <memory> // this is for the safe pointer wrapper
#include <position.h>
#include <velocity.h>
#include <worldParams.h>

class recShot {
public:
    
    static std::vector<std::shared_ptr<position>> launchAir(const std::shared_ptr<position>& launchPos,
                                                            const std::shared_ptr<position>& tgtPos,
                                                            const double& v0, const double& elDeg,
                                                            const double& timeStep, const double& mortSigma, const double& mortorMass);
    
    static void recPos(std::vector<std::shared_ptr<position>>& vec, const std::shared_ptr<position>& pos);
    
    static void zenCheck(const std::shared_ptr<position>& currPos, const std::shared_ptr<position>& tgtPos,
                         std::shared_ptr<position>& prevPos, bool& floorMove, double& floorZ, bool& zenFail);
    
    
};//End rec shot class

#endif /* defined(__gputest1__recShot__) */
