//
// Created by yeyinong on 2023/8/12.
//

#include "math.h"
#include "parameters.h"

const double gravity = 9.81;

const bool LEFTLEG = 1;
const bool RIGHTLEG = 0;

// parameters of RobBIE mechanical
const double m_RobBIE=15.49;

const double Lf=0.42;
const double Lt=0.42;
const double LT=0.02;
const double lShoulder = 0.;
const double offsetY = 0.;
const double offsetZ = 0.;

const double dx1_L=0.043001, dy1_L=0.15,dz1_L=-0.015843;//0.043001 0.15 -0.015843
const double dx2_L=0.098, dy2_L=0.0,dz2_L=-0.0345;//0.098 0 -0.0345
const double dx3_L=0, dy3_L=0.08,dz3_L=-0.088;
const double dx4_L=0.0171743, dy4_L=0,dz4_L=-0.40964;//0.0171743 0 -0.40964
const double dx5_L=0.0159183, dy5_L=-0.054,dz5_L=-0.4217;//0.0159183 -0.054 -0.4217

/*const double dx1_R=0., dy1_R=-0.095,dz1_R=-0.0325;
const double dx2_R=-0.0775, dy2_R=0.,dz2_R=-0.083;
const double dx3_R=0.0775, dy3_R=-0.0275,dz3_R=0.;
const double dx4_R=0., dy4_R=0.02255,dz4_R=-0.2;
const double dx5_R=0., dy5_R=0.005,dz5_R=-0.25;*/

const double dx1_R=0.043001, dy1_R=-0.15,dz1_R=-0.015843;
const double dx2_R=0.098, dy2_R=0.,dz2_R=-0.0345;
const double dx3_R=0, dy3_R=0.08,dz3_R=-0.088;
const double dx4_R=0.01814, dy4_R=0,dz4_R=-0.4096;
const double dx5_R=0.014952, dy5_R=-0.053,dz5_R=-0.42174;

const double offsetQAbAd = 0;
const double offsetQFem = 0;
const double offsetQTib = 0;
const double offsetY_foot2Com_left = -0.125;
const double offsetY_foot2Com_right = 0.125;
const double offsetX_Com = -0.;

/*
const double kp_Height = 600;
const double kd_Height = 50;
const double kp_Torso_Roll = 200;
const double kd_Torso_Roll = 5;
const double kp_Torso_Pitch = 200;
const double kd_Torso_Pitch = 5;

const double kp_Yaw_Standing = 400;
const double kd_Yaw_Standing = 0.35;
const double kp_Yaw_Swing = 400;
const double kd_Yaw_Swing = 0.35;
const double kp_Abd_Swing = 400;
const double kd_Abd_Swing = 0.35;
const double kp_Hip_Swing = 400;
const double kd_Hip_Swing = 0.35;
const double kp_Knee_Swing = 400;
const double kd_Knee_Swing = 0.35;
*/
const double kp_Height = 600;
const double kd_Height = 50;
const double kp_Torso_Roll = 400;
const double kd_Torso_Roll = 15;
const double kp_Torso_Pitch = 400;
const double kd_Torso_Pitch = 15;

const double kp_Yaw_Standing = 10;
const double kd_Yaw_Standing = 0.35;
const double kp_Yaw_Swing = 10;
const double kd_Yaw_Swing = 0.35;
const double kp_Abd_Swing = 10;
const double kd_Abd_Swing = 0.35;
const double kp_Hip_Swing = 10;
const double kd_Hip_Swing = 0.35;
const double kp_Knee_Swing = 10;
const double kd_Knee_Swing = 0.35;

const double kp_Abd_Swing_cmd = 100;
const double kd_Abd_Swing_cmd = 5;
const double kp_Hip_Swing_cmd = 50;
const double kd_Hip_Swing_cmd = 2;
const double kp_Knee_Swing_cmd = 50;
const double kd_Knee_Swing_cmd = 2;

const double posSWAD_Compensate_Left = 0./180.0*M_PI;
const double posSWAD_Compensate_Right = -0./180.0*M_PI;

const double kp_vx = 0.2;
const double kd_vx = 10;
const double kp_SWY = 0.1;
const double cmdWs = 0.17 * 2;

const double filterVel = 0.2;
const double filterGRF = 0.2;

// joint torque limit
const double torqueLB = -1000;
const double torqueUB = 1000;

// swing foot limit
const double posSWFootXLB = -0.5;
const double posSWFootXUB = 0.5;
const double posSWFootYLB_Right = -0.2;
const double posSWFootYUB_Right = -0.;
const double posSWFootYLB_Left = 0.;
const double posSWFootYUB_Left = 0.2;

const double stepTime = 0.3;

/// Kalman Filter
const double c_Cov_rpx_LRToe_body = 0.01*0.01;
const double c_Cov_rpy_LRToe_body = 0.01*0.01;
const double c_Cov_rpz_LRToe_body = 0.01*0.01;
const double c_Cov_LinearAccelerator_xy = 0.05*0.05;
const double c_Cov_LinearAccelerator_z = 0.05*0.05;

const double c_vx_slide_toe_1 = 2;
const double c_vx_slide_toe_2 = 0.02;
const double c_Tx_slide_toe_2 = 0.05;

const double c_vy_slide_toe_1 = 2;
const double c_vy_slide_toe_2 = 0.02;
const double c_Ty_slide_toe_2 = 0.05;

const double c_vz_slide_toe_1 = 2;
const double c_vz_slide_toe_2 = 0.02;
const double c_Tz_slide_toe_2 = 0.05;
