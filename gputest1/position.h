//
//  position.h
//  gputest1
//
//  Created by Chris Nicholas on 01/06/2015.
//  Copyright (c) 2015 Chris Nicholas. All rights reserved.
//

#ifndef gputest1_position_h
#define gputest1_position_h

class position {
public:
    double x;
    double y;
    double z;
    position(const double& x = 0, const double& y = 0, const double& z = 0) : x{ x }, y{ y }, z{ z } {}
};

#endif
