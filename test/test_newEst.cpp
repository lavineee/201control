//
// Created by yeyinong on 2023/8/12.
//
#include "Eigen/Core"
#include "iostream"

using namespace std;

double Lf = 0.212;
double Lt = 0.185;
double lShoulder = 0.075;
double offsetY = 0.05;
double offsetZ = 0.08;

double deg2rad(double deg){
    return deg/180.*M_PI;
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

int main(){
    Eigen::Matrix<double, 3, 1> leftFoot, rightFoot;
    Eigen::Matrix<double, 3, 1> qLeft, qRight;
    qLeft << deg2rad(-10), deg2rad(30), deg2rad(-60);
    q2Foot_Right(leftFoot, qLeft);
    cout << leftFoot << endl << endl;
    foot2Q_Right(qLeft, leftFoot);
    cout << qLeft << endl<<endl;

    cout << deg2rad(10) << endl;
    cout << deg2rad(30) << endl;
    cout << deg2rad(-60) << endl;
}