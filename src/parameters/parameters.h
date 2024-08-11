//
// Created by yeyinong on 2023/8/12.
//

#ifndef MUJOCOTEST3D_PARAMETERS_H
#define MUJOCOTEST3D_PARAMETERS_H

//
// Created by yeyinong on 2023/8/12.
//

#include "parameters.h"

extern const double gravity;// 9.81;

extern const bool LEFTLEG;// 0;
extern const bool RIGHTLEG;// 1;

extern const double m_RobBIE; //=6;

extern const double Lf;
extern const double Lt;
extern const double LT;
extern const double lShoulder;

extern const double dx1_L, dy1_L,dz1_L;
extern const double dx2_L, dy2_L,dz2_L;
extern const double dx3_L, dy3_L,dz3_L;
extern const double dx4_L, dy4_L,dz4_L;
extern const double dx5_L, dy5_L,dz5_L;

extern const double dx1_R, dy1_R,dz1_R;
extern const double dx2_R, dy2_R,dz2_R;
extern const double dx3_R, dy3_R,dz3_R;
extern const double dx4_R, dy4_R,dz4_R;
extern const double dx5_R, dy5_R,dz5_R;

extern const double offsetY;// 0.5;
extern const double offsetZ;// 0.5;

extern const double offsetQAbAd;// 0;
extern const double offsetQFem;// 0;
extern const double offsetQTib;// 0;
extern const double offsetY_foot2Com_left;// -0.108;
extern const double offsetY_foot2Com_right;// 0.108;
extern const double offsetX_Com;// -0.;

extern const double kp_Height;// 1500;
extern const double kd_Height;// 100;
extern const double kp_Torso_Roll;// 500;
extern const double kd_Torso_Roll;// 5;
extern const double kp_Torso_Pitch;// 40;
extern const double kd_Torso_Pitch;// 1;

extern const double kp_Yaw_Standing;
extern const double kd_Yaw_Standing;
extern const double kp_Yaw_Swing;
extern const double kd_Yaw_Swing;
extern const double kp_Abd_Swing;// 15;
extern const double kd_Abd_Swing;// 0.35;
extern const double kp_Hip_Swing;// 10;
extern const double kd_Hip_Swing;// 0.5;
extern const double kp_Knee_Swing;// 10;
extern const double kd_Knee_Swing;// 0.5;

extern const double kp_Abd_Swing_cmd;// 100;
extern const double kd_Abd_Swing_cmd;// 5;
extern const double kp_Hip_Swing_cmd;// 50;
extern const double kd_Hip_Swing_cmd;// 2;
extern const double kp_Knee_Swing_cmd;// 50;
extern const double kd_Knee_Swing_cmd;// 2;

extern const double posSWAD_Compensate_Left;// 0./180.0*M_PI;
extern const double posSWAD_Compensate_Right;// -0./180.0*M_PI;

extern const double kd_vx;
extern const double kp_vx;// 0.4;
extern const double kp_SWY;// 0.05;

extern const double filterVel;// 0.02;
extern const double filterGRF;// 0.2;

// joint torque limit
extern const double torqueLB;// -18;
extern const double torqueUB;// 18;

// swing foot limit
extern const double posSWFootXLB;// -0.2;
extern const double posSWFootXUB;// 0.2;
extern const double posSWFootYLB_Right;// = -0.325;
extern const double posSWFootYUB_Right;// = -0.03;
extern const double posSWFootYLB_Left;// = 0.03;
extern const double posSWFootYUB_Left;// = 0.325;


extern const double stepTime;// 0.35;

extern const double cmdWs;// = 0.125*2;

/// Kalman Filter
extern const double c_Cov_rpx_LRToe_body;
extern const double c_Cov_rpy_LRToe_body;
extern const double c_Cov_rpz_LRToe_body;
extern const double c_Cov_LinearAccelerator_xy;
extern const double c_Cov_LinearAccelerator_z;

extern const double c_vx_slide_toe_1;
extern const double c_vx_slide_toe_2;
extern const double c_Tx_slide_toe_2;

extern const double c_vy_slide_toe_1;
extern const double c_vy_slide_toe_2;
extern const double c_Ty_slide_toe_2;

extern const double c_vz_slide_toe_1;
extern const double c_vz_slide_toe_2;
extern const double c_Tz_slide_toe_2;


#endif //MUJOCOTEST3D_PARAMETERS_H
