/*
//
// Created by Administrator on 2023/8/12.
//
#include "eigen3/Eigen/Core"
#include "math.h"
#include "parameters/parameters.h"
#include "robbie/calculateRobBIE.h"

double deg2rad(double deg){
    return deg/180.*M_PI;
}

using namespace std;
void q2Foot_Left(Eigen::Matrix<double, 3, 1>& leftFoot, const Eigen::Matrix<double, 3, 1>& qLeft){
    double q0 = qLeft(0);
    double q1 = qLeft(1);
    double q2 = qLeft(2);

    double x = -Lf* sin(q1) - Lt* sin(q1+q2);
    double zOrg = Lf* cos(q1) + Lt* cos(q1+q2);
    double lY = sqrt(zOrg*zOrg + lShoulder*lShoulder);
    double theta = atan(lShoulder/zOrg);
    double y = lY * sin(q0+theta) + offsetY;
    double z = -lY* cos(q0+theta) - offsetZ;
    leftFoot << x, y, z;
}

void q2Base_Left(Eigen::Matrix<double, 3, 1>& leftBase, const Eigen::Matrix<double, 3, 1>& qLeft){
    q2Foot_Left(leftBase, qLeft);
    leftBase = - leftBase;
}

void q2Foot_Right(Eigen::Matrix<double, 3, 1>& rightFoot, const Eigen::Matrix<double, 3, 1>& qRight){
    double q0 = qRight(0);
    double q1 = qRight(1);
    double q2 = qRight(2);

    double x = -Lf* sin(q1) - Lt* sin(q1+q2);
    double zOrg = Lf* cos(q1) + Lt* cos(q1+q2);
    double lY = sqrt(zOrg*zOrg + lShoulder*lShoulder);
    double theta = atan(lShoulder/zOrg);
    double y = lY * sin(q0-theta) - offsetY;
    double z = -lY* cos(q0-theta) - offsetZ;
    rightFoot << x, y, z;
}

void q2Base_Right(Eigen::Matrix<double, 3, 1>& rightBase, const Eigen::Matrix<double, 3, 1>& qRight){
    q2Foot_Right(rightBase, qRight);
    rightBase = -rightBase;
}

void foot2Q_Left(Eigen::Matrix<double, 3, 1>& qLeft, const Eigen::Matrix<double, 3, 1>& footLeft){
    double x = footLeft(0);
    double y = footLeft(1);
    double z = footLeft(2);

    y =  y-offsetY;
    z = -z - offsetZ;
    double lY = sqrt(y*y + z*z);
    double theta1 = atan(y/z);
    double theta2 = acos(lShoulder/lY);
    double q0 = theta2 - (M_PI/2-theta1);

    double z0 = sqrt(lY*lY -lShoulder*lShoulder);
    double lx = sqrt(z0*z0 + x*x);
    double theta3 = acos((Lf*Lf + Lt*Lt - lx*lx)/(2*Lf*Lt));
    double theta4 = acos((Lf*Lf + lx*lx - Lt*Lt)/(2*Lf*lx));

    double q1 = theta4 - atan(x/z0);
    double q2 = theta3 - M_PI;

    qLeft << q0, q1, q2;
}

void foot2Q_Right(Eigen::Matrix<double, 3, 1>& qRight, const Eigen::Matrix<double, 3, 1>& footRight){
    Eigen::Matrix<double, 3, 1> footTrans = footRight;
    footTrans(1) = -footTrans(1);
    foot2Q_Left(qRight, footTrans);
    qRight(0) = -qRight(0);
}

void calculateJb_Left(Eigen::Matrix<double, 3, 3>& JbLeft,const Eigen::Matrix<double, 3, 1>& qLeft){
    double q0 = qLeft(0);
    double q1 = qLeft(1);
    double q2 = qLeft(2);
    JbLeft << 0,-(Lf*cos(q1)) - Lt*cos(q1 + q2),-(Lt*cos(q1 + q2)),
    sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),
            -((lShoulder*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))*(-(Lf*sin(q1)) - Lt*sin(q1 + q2)))/
              (pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)*(1 + pow(lShoulder,2)/pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)))) +
            ((Lf*cos(q1) + Lt*cos(q1 + q2))*(-(Lf*sin(q1)) - Lt*sin(q1 + q2))*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))/sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)),
            (lShoulder*Lt*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))*sin(q1 + q2))/
            (pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)*(1 + pow(lShoulder,2)/pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))) -
            (Lt*(Lf*cos(q1) + Lt*cos(q1 + q2))*sin(q1 + q2)*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))/sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)),
    sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),
            -(((Lf*cos(q1) + Lt*cos(q1 + q2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))*(-(Lf*sin(q1)) - Lt*sin(q1 + q2)))/sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))) -
            (lShoulder*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*(-(Lf*sin(q1)) - Lt*sin(q1 + q2))*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))/
            (pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)*(1 + pow(lShoulder,2)/pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))),
            (Lt*(Lf*cos(q1) + Lt*cos(q1 + q2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))*sin(q1 + q2))/sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)) +
            (lShoulder*Lt*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q1 + q2)*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))/
            (pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)*(1 + pow(lShoulder,2)/pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)));
}

void calculateJb_Right(Eigen::Matrix<double, 3, 3>& JbRight,const Eigen::Matrix<double, 3, 1>& qRight){
    double q0 = qRight(0);
    double q1 = qRight(1);
    double q2 = qRight(2);

    JbRight <<
        0,-(Lf*cos(q1)) - Lt*cos(q1 + q2),-(Lt*cos(q1 + q2)),
        sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),
            (lShoulder*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))*(-(Lf*sin(q1)) - Lt*sin(q1 + q2)))/
            (pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)*(1 + pow(lShoulder,2)/pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))) +
            ((Lf*cos(q1) + Lt*cos(q1 + q2))*(-(Lf*sin(q1)) - Lt*sin(q1 + q2))*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))/sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)),
            -((lShoulder*Lt*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))*sin(q1 + q2))/
            (pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)*(1 + pow(lShoulder,2)/pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)))) -
            (Lt*(Lf*cos(q1) + Lt*cos(q1 + q2))*sin(q1 + q2)*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))/sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)),
       sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),
            -(((Lf*cos(q1) + Lt*cos(q1 + q2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))*(-(Lf*sin(q1)) - Lt*sin(q1 + q2)))/sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))) +
            (lShoulder*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*(-(Lf*sin(q1)) - Lt*sin(q1 + q2))*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))/
            (pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)*(1 + pow(lShoulder,2)/pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))),
            (Lt*(Lf*cos(q1) + Lt*cos(q1 + q2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))*sin(q1 + q2))/sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)) -
            (lShoulder*Lt*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q1 + q2)*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))/
            (pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)*(1 + pow(lShoulder,2)/pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)));
}


void calculateJb_Virtual2Q_Left(Eigen::Matrix<double, 3, 3>& JbLeft, const Eigen::Matrix<double, 3, 1>& qLeft){
    double q0 = qLeft(0);
    double q1 = qLeft(1);
    double q2 = qLeft(2);
    JbLeft <<
           (sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*(offsetY*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) - offsetZ*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))))/
           sqrt(pow(Lf,2) + pow(lShoulder,2) + pow(Lt,2) + pow(offsetY,2) + pow(offsetZ,2) + 2*Lf*Lt*cos(q2) + 2*offsetZ*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                2*offsetY*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))),
            -((Lf*sin(q1) + Lt*sin(q1 + q2))*(-2*lShoulder*offsetY*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) + Lf*offsetZ*cos(q0 - q1 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                              Lf*offsetZ*cos(q0 + q1 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) + Lt*offsetZ*cos(q0 - q1 - q2 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) + Lt*offsetZ*cos(q0 + q1 + q2 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                              2*lShoulder*offsetZ*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) + Lf*offsetY*sin(q0 - q1 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) + Lf*offsetY*sin(q0 + q1 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                              Lt*offsetY*sin(q0 - q1 - q2 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) + Lt*offsetY*sin(q0 + q1 + q2 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))))/
            (2.*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sqrt(pow(Lf,2) + pow(lShoulder,2) + pow(Lt,2) + pow(offsetY,2) + pow(offsetZ,2) + 2*Lf*Lt*cos(q2) +
                                                                                  2*offsetZ*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                                                                  2*offsetY*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))),
            (Lt*(cos(q1 + q2)*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*(Lf*sin(q1) + Lt*sin(q1 + q2)) -
                 (offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))*sin(q1 + q2)*
                 ((Lf*cos(q1) + Lt*cos(q1 + q2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) + lShoulder*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))) +
                 sin(q1 + q2)*(lShoulder*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) - (Lf*cos(q1) + Lt*cos(q1 + q2))*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))*
                 (offsetY + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))))/
            (sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sqrt(pow(Lf,2) + pow(lShoulder,2) + pow(Lt,2) + pow(offsetY,2) + pow(offsetZ,2) + 2*Lf*Lt*cos(q2) +
                                                                               2*offsetZ*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                                                               2*offsetY*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))),
            (2*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*(sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)) + offsetZ*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                                                             offsetY*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))))/
            (pow(Lf,2) + 2*pow(lShoulder,2) + pow(Lt,2) + 2*pow(offsetY,2) + 2*pow(offsetZ,2) + pow(Lf,2)*cos(2*q1) + 2*Lf*Lt*cos(q2) + pow(Lt,2)*cos(2*(q1 + q2)) + 2*Lf*Lt*cos(2*q1 + q2) +
             4*offsetZ*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
             4*offsetY*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))),
            (2*(Lf*sin(q1) + Lt*sin(q1 + q2))*((lShoulder*offsetZ + Lf*offsetY*cos(q1) + Lt*offsetY*cos(q1 + q2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                               lShoulder*sqrt(pow(lShoulder,2) + pow(Lf,2)*pow(cos(q1),2) + 2*Lf*Lt*cos(q1)*cos(q1 + q2) + pow(Lt,2)*pow(cos(q1 + q2),2))*pow(cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2) +
                                               sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))*(lShoulder*offsetY - Lf*offsetZ*cos(q1) - Lt*offsetZ*cos(q1 + q2) +
                                                                                                         lShoulder*sqrt(pow(lShoulder,2) + pow(Lf,2)*pow(cos(q1),2) + 2*Lf*Lt*cos(q1)*cos(q1 + q2) + pow(Lt,2)*pow(cos(q1 + q2),2))*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))))/
            (sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*(pow(Lf,2) + 2*pow(lShoulder,2) + pow(Lt,2) + 2*pow(offsetY,2) + 2*pow(offsetZ,2) + pow(Lf,2)*cos(2*q1) + 2*Lf*Lt*cos(q2) + pow(Lt,2)*cos(2*(q1 + q2)) +
                                                                           2*Lf*Lt*cos(2*q1 + q2) + 4*offsetZ*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                                                           4*offsetY*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))),
            (2*Lt*sin(q1 + q2)*((lShoulder*offsetZ + Lf*offsetY*cos(q1) + Lt*offsetY*cos(q1 + q2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                lShoulder*sqrt(pow(lShoulder,2) + pow(Lf,2)*pow(cos(q1),2) + 2*Lf*Lt*cos(q1)*cos(q1 + q2) + pow(Lt,2)*pow(cos(q1 + q2),2))*pow(cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2) +
                                sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))*(lShoulder*offsetY - Lf*offsetZ*cos(q1) - Lt*offsetZ*cos(q1 + q2) +
                                                                                          lShoulder*sqrt(pow(lShoulder,2) + pow(Lf,2)*pow(cos(q1),2) + 2*Lf*Lt*cos(q1)*cos(q1 + q2) + pow(Lt,2)*pow(cos(q1 + q2),2))*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))))/
            (sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*(pow(Lf,2) + 2*pow(lShoulder,2) + pow(Lt,2) + 2*pow(offsetY,2) + 2*pow(offsetZ,2) + pow(Lf,2)*cos(2*q1) + 2*Lf*Lt*cos(q2) + pow(Lt,2)*cos(2*(q1 + q2)) +
                                                                           2*Lf*Lt*cos(2*q1 + q2) + 4*offsetZ*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                                                           4*offsetY*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))),
            (sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*(Lf*sin(q1) + Lt*sin(q1 + q2))*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))/
            (pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2)*
             (1 + pow(Lf*sin(q1) + Lt*sin(q1 + q2),2)/pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2))),
            ((Lf*cos(q1) + Lt*cos(q1 + q2))*(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))) +
             (pow(Lf*sin(q1) + Lt*sin(q1 + q2),2)*(Lf*cos(q1)*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) + Lt*cos(q1 + q2)*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                                   lShoulder*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))))/sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)))/
            (pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2)*
             (1 + pow(Lf*sin(q1) + Lt*sin(q1 + q2),2)/pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2))),
            (Lt*cos(q1 + q2)*(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))) +
             (Lt*sin(q1 + q2)*(Lf*sin(q1) + Lt*sin(q1 + q2))*((Lf*cos(q1) + Lt*cos(q1 + q2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) + lShoulder*sin(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))))/
             sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)))/
            (pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2)*
             (1 + pow(Lf*sin(q1) + Lt*sin(q1 + q2),2)/pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 + atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2)));
}

void calculateJb_Virtual2Q_Right(Eigen::Matrix<double, 3, 3>& JbRight, const Eigen::Matrix<double, 3, 1>& qRight){
    double q0 = qRight(0);
    double q1 = qRight(1);
    double q2 = qRight(2);
    JbRight <<
            -((sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*(offsetY*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) + offsetZ*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))))/
              sqrt(pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2) + pow(Lf*sin(q1) + Lt*sin(q1 + q2),2) +
                   pow(offsetY - sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2))),
            ((Lf*sin(q1) + Lt*sin(q1 + q2))*(2*lShoulder*offsetY*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) - Lf*offsetZ*cos(q0 - q1 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) -
                                             Lf*offsetZ*cos(q0 + q1 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) - Lt*offsetZ*cos(q0 - q1 - q2 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) - Lt*offsetZ*cos(q0 + q1 + q2 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                             2*lShoulder*offsetZ*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) + Lf*offsetY*sin(q0 - q1 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) + Lf*offsetY*sin(q0 + q1 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                             Lt*offsetY*sin(q0 - q1 - q2 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) + Lt*offsetY*sin(q0 + q1 + q2 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))))/
            (2.*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sqrt(pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2) +
                                                                                  pow(Lf*sin(q1) + Lt*sin(q1 + q2),2) + pow(offsetY - sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2))),
            (Lt*(cos(q1 + q2)*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*(Lf*sin(q1) + Lt*sin(q1 + q2)) -
                 (offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))*sin(q1 + q2)*
                 ((Lf*cos(q1) + Lt*cos(q1 + q2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) - lShoulder*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))) +
                 sin(q1 + q2)*(lShoulder*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) + (Lf*cos(q1) + Lt*cos(q1 + q2))*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))*
                 (offsetY - sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))))/
            (sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sqrt(pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2) +
                                                                               pow(Lf*sin(q1) + Lt*sin(q1 + q2),2) + pow(offsetY - sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2))),
            (2*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*(sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)) + offsetZ*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) -
                                                                             offsetY*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))))/
            (pow(Lf,2) + 2*pow(lShoulder,2) + pow(Lt,2) + 2*pow(offsetY,2) + 2*pow(offsetZ,2) + pow(Lf,2)*cos(2*q1) + 2*Lf*Lt*cos(q2) + pow(Lt,2)*cos(2*(q1 + q2)) + 2*Lf*Lt*cos(2*q1 + q2) +
             4*offsetZ*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) -
             4*offsetY*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))),
            (-2*(Lf*sin(q1) + Lt*sin(q1 + q2))*((lShoulder*offsetZ + Lf*offsetY*cos(q1) + Lt*offsetY*cos(q1 + q2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                                lShoulder*sqrt(pow(lShoulder,2) + pow(Lf,2)*pow(cos(q1),2) + 2*Lf*Lt*cos(q1)*cos(q1 + q2) + pow(Lt,2)*pow(cos(q1 + q2),2))*pow(cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2) +
                                                sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))*(-(lShoulder*offsetY) + Lf*offsetZ*cos(q1) + Lt*offsetZ*cos(q1 + q2) +
                                                                                                          lShoulder*sqrt(pow(lShoulder,2) + pow(Lf,2)*pow(cos(q1),2) + 2*Lf*Lt*cos(q1)*cos(q1 + q2) + pow(Lt,2)*pow(cos(q1 + q2),2))*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))))/
            (sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*(pow(Lf,2) + 2*pow(lShoulder,2) + pow(Lt,2) + 2*pow(offsetY,2) + 2*pow(offsetZ,2) + pow(Lf,2)*cos(2*q1) + 2*Lf*Lt*cos(q2) + pow(Lt,2)*cos(2*(q1 + q2)) +
                                                                           2*Lf*Lt*cos(2*q1 + q2) + 4*offsetZ*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) -
                                                                           4*offsetY*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))),
            (-2*Lt*sin(q1 + q2)*((lShoulder*offsetZ + Lf*offsetY*cos(q1) + Lt*offsetY*cos(q1 + q2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) +
                                 lShoulder*sqrt(pow(lShoulder,2) + pow(Lf,2)*pow(cos(q1),2) + 2*Lf*Lt*cos(q1)*cos(q1 + q2) + pow(Lt,2)*pow(cos(q1 + q2),2))*pow(cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2) +
                                 sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))*(-(lShoulder*offsetY) + Lf*offsetZ*cos(q1) + Lt*offsetZ*cos(q1 + q2) +
                                                                                           lShoulder*sqrt(pow(lShoulder,2) + pow(Lf,2)*pow(cos(q1),2) + 2*Lf*Lt*cos(q1)*cos(q1 + q2) + pow(Lt,2)*pow(cos(q1 + q2),2))*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))))/
            (sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*(pow(Lf,2) + 2*pow(lShoulder,2) + pow(Lt,2) + 2*pow(offsetY,2) + 2*pow(offsetZ,2) + pow(Lf,2)*cos(2*q1) + 2*Lf*Lt*cos(q2) + pow(Lt,2)*cos(2*(q1 + q2)) +
                                                                           2*Lf*Lt*cos(2*q1 + q2) + 4*offsetZ*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) -
                                                                           4*offsetY*sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))),
            (sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*(Lf*sin(q1) + Lt*sin(q1 + q2))*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))))/
            (pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2)*
             (1 + pow(Lf*sin(q1) + Lt*sin(q1 + q2),2)/pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2))),
            ((Lf*cos(q1) + Lt*cos(q1 + q2))*(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))) +
             (pow(Lf*sin(q1) + Lt*sin(q1 + q2),2)*((Lf*cos(q1) + Lt*cos(q1 + q2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) - lShoulder*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))))/
             sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)))/
            (pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2)*
             (1 + pow(Lf*sin(q1) + Lt*sin(q1 + q2),2)/pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2))),
            (Lt*cos(q1 + q2)*(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))) +
             (Lt*sin(q1 + q2)*(Lf*sin(q1) + Lt*sin(q1 + q2))*((Lf*cos(q1) + Lt*cos(q1 + q2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))) - lShoulder*sin(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2))))))/
             sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2)))/
            (pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2)*
             (1 + pow(Lf*sin(q1) + Lt*sin(q1 + q2),2)/pow(offsetZ + sqrt(pow(lShoulder,2) + pow(Lf*cos(q1) + Lt*cos(q1 + q2),2))*cos(q0 - atan(lShoulder/(Lf*cos(q1) + Lt*cos(q1 + q2)))),2)));


}*/
