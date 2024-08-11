//
// Created by yeyinong on 2023/2/20.
//
#include "iostream"

#include "math.h"
#include "vector"
#include "stdint.h"
#include "Eigen/Core"

using namespace std;

const double offset_qFem = 3.141592654;
const double offset_qTib = 0;

const double ti = 0;
const double tf = 0.3905;
const double g = 9.81;

const double MT = 3.34496;
const double Mf = 0.90568;
const double Mt = 0.13622;

const double LT = 0.16638;
const double Lf = 0.2115;
const double Lt = 0.225;

const double pzT = 0.02736;
const double pxT = 0.;
const double pzf = 0.0208;
const double pzt = 0.09473;

vector<double> calculateSSP(double y0, double dy0, double z0, double t){
    double lambda = sqrt(g/z0);
    double c1 = 1./2*(y0+dy0/lambda);
    double c2 = 1./2*(y0-dy0/lambda);
    double ySSP = c1* pow(M_E,lambda*t)+c2* pow(M_E,-lambda*t);
    double dySSP = lambda*(c1* pow(M_E,lambda*t)-c2* pow(M_E,-lambda*t));
    return vector<double>{ySSP, dySSP};
}

int main(void){
//    Eigen::Matrix<double, 3, 1> test;
//    test << 1,1,1;
//    cout << -test << endl;
//    double base_state_quat_w = 0.96;
//    double base_state_quat_x = -0.12;
//    double base_state_quat_y = 0.04;
//    double base_state_quat_z = -0.25;
////    double thetaY = asin(2 * base_state_quat_w * base_state_quat_y - base_state_quat_z * base_state_quat_x);
////    double thetaX = atan2(2*(base_state_quat_w*base_state_quat_x+base_state_quat_y*base_state_quat_z), 1-2*(base_state_quat_x*base_state_quat_x+base_state_quat_y*base_state_quat_y));
//    double thetaX = atan2(2*(base_state_quat_w*base_state_quat_x+base_state_quat_y*base_state_quat_z), 1-2*(base_state_quat_x*base_state_quat_x+base_state_quat_y*base_state_quat_y));
//    double thetaY = asin(2 * (base_state_quat_w * base_state_quat_y - base_state_quat_x * base_state_quat_z));
//    cout << "thetaX: " << thetaX/M_PI*180 << ", thetaY: " << thetaY/M_PI*180 << endl;

    // q0:0.96, q1:-0.12, q2:0.04, q3:-0.25 |Roll:-14.65, Pitch:0.87, Yaw:-28.85

}