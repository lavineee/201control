//
// Created by yeyinong on 2023/8/12.
//

#ifndef MUJOCOTEST3D_CALCULATEROBBIE_H
#define MUJOCOTEST3D_CALCULATEROBBIE_H
//
// Created by Administrator on 2023/8/12.
//

using namespace std;
void q2Foot_Left(Eigen::Matrix<double, 3, 1>& leftFoot, const Eigen::Matrix<double, 3, 1>& qLeft);

void q2Base_Left(Eigen::Matrix<double, 3, 1>& leftBase, const Eigen::Matrix<double, 3, 1>& qLeft);

void q2Foot_Right(Eigen::Matrix<double, 3, 1>& rightFoot, const Eigen::Matrix<double, 3, 1>& qRight);

void q2Base_Right(Eigen::Matrix<double, 3, 1>& rightBase, const Eigen::Matrix<double, 3, 1>& qRight);

void foot2Q_Left(Eigen::Matrix<double, 3, 1>& qLeft, const Eigen::Matrix<double, 3, 1>& footLeft);

void foot2Q_Right(Eigen::Matrix<double, 3, 1>& qRight, const Eigen::Matrix<double, 3, 1>& footRight);

void calculateJb_Left(Eigen::Matrix<double, 3, 3>& JbLeft,const Eigen::Matrix<double, 3, 1>& qLeft);

void calculateJb_Right(Eigen::Matrix<double, 3, 3>& JbRight,const Eigen::Matrix<double, 3, 1>& qRight);

void calculateJb_Virtual2Q_Left(Eigen::Matrix<double, 3, 3>& JbLeft, const Eigen::Matrix<double, 3, 1>& qLeft);

void calculateJb_Virtual2Q_Right(Eigen::Matrix<double, 3, 3>& JbRight, const Eigen::Matrix<double, 3, 1>& qRight);
#endif //MUJOCOTEST3D_CALCULATEROBBIE_H
