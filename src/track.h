//
// Created by yeyinong on 2022/4/25.
//

#ifndef MUJOCOTEST_TRACK_H
#define MUJOCOTEST_TRACK_H

#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>
//#include "matplotlibcpp.h"

using namespace std;
using namespace cv;

#define NUMPOINTS 1000



//Point2d bezierPoint(vector<Point2d> src, float t){
//    Point2d res = src[0]*(pow(1-t,3)) + 3*src[1]*t*(pow(1-t,2)) + 3*src[2]* pow(t, 2)*(1-t) + src[3]* pow(t, 3);
////    res.x -= t*(src[3].x-src[0].x)/2;
//    return res;
//}

//vector<cv::Point2d> bezierCurve(vector<Point2d> src, int numPoints){
//    vector<Point2d> res;
//    if(src.size()!=4)
//        exit(-1);
//    for (int i = 0; i < numPoints; ++i) {
//        res.push_back(bezierPoint(src, i/double ((numPoints-1))));
//    }
//    return res;
//}
//
//vector<cv::Point2d> backCurve(Point2d begin, int numPoints){
//    vector<Point2d> res;
//    for (int i = 0; i < numPoints; ++i) {
//        res.push_back(Point2d(begin.x+0.3-0.3/numPoints*(i+1), begin.y));
//    }
//    return res;
//}
#endif //MUJOCOTEST_TRACK_H
