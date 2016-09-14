#!/bin/sh

#  build.sh
#  gputest1
#
#  Created by Chris Nicholas on 13/09/2016.
#  Copyright © 2016 Chris Nicholas. All rights reserved.


sudo g++ -c utils.cpp
sudo g++ -c main.cpp
g++ main.o -lOpenCL -lgdal -o gpublosmain utils.o
