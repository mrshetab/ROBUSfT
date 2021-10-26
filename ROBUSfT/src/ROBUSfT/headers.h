#ifndef HEADERS_H_INCLUDED
#define HEADERS_H_INCLUDED

#define _USE_MATH_DEFINES

#include <iostream>
#include <memory>
#include <cmath>
#include <vector>
#include <iterator> 
#include <fstream>
#include <functional>   
#include <algorithm>    
#include <typeinfo>       
#include <chrono>
#include <string>
#include <numeric>
#include <stdlib.h>
#include <math.h>

#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/xfeatures2d.hpp"
#include <opencv2/xfeatures2d/nonfree.hpp>

#include <Eigen/Core>

#include <popsift/common/device_prop.h>
#include <popsift/features.h>
#include <popsift/popsift.h>
#include <popsift/sift_conf.h>
#include <popsift/sift_config.h>
#include <popsift/version.hpp>

using namespace cv;
using namespace cv::xfeatures2d;
using namespace std;

typedef std::vector<double> DoubleVector1D;
typedef std::vector<std::vector<double>> DoubleVector2D;
typedef std::vector<std::vector<std::vector<double>>> DoubleVector3D; 

#endif // HEADERS_H_INCLUDED
