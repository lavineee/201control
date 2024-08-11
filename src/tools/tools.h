//
// Created by yeyinong on 2023/3/29.
//

#ifndef MUJOCOTEST3D_TOOLS_H
#define MUJOCOTEST3D_TOOLS_H

#include "vector"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Geometry"
#include "robbie/Robbie.h"
#include "mujoco/mujoco.h"

using namespace std;

vector<double> bezierPoint(double px1, double py1, double px2, double py2, double px3, double py3, double px4, double py4, double tau);

Eigen::Matrix<double, 3, 1> bezierPoint3D(Eigen::Matrix<double, 3, 1> startPoint, Eigen::Matrix<double, 3, 1> controlPoint, Eigen::Matrix<double, 3, 1> endPoint,  double tau);

Eigen::Matrix<double, 3, 1> bezierPoint3D(Eigen::Matrix<double, 3, 1> startPoint, Eigen::Matrix<double, 3, 1> controlPoint1, Eigen::Matrix<double, 3, 1> controlPoint2, Eigen::Matrix<double, 3, 1> endPoint,  double tau);

void getSensorMessage(const mjData *d, vector<double>& floatMessage, vector<double>& motorMessage);
#endif //MUJOCOTEST3D_TOOLS_H
