//
// Created by yeyinong on 2023/2/23.
//

#ifndef MUJOCOTEST2D_SIMPLEMODELCONTROLLER_H
#define MUJOCOTEST2D_SIMPLEMODELCONTROLLER_H

#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/Geometry"
#include <math.h>
#include "matplotlibcpp.h"
#include "robbie/Robbie.h"
#include "tools/tools.h"

namespace plt = matplotlibcpp;
using namespace std;

const double Lf=0.211;
const double Lt=0.208;
const double mr=5.35;
const double pi=3.1416;
const double g=9.81;
const double mT=3.34496;
const double LT=0.0267419;

//const double offsetQAbAd = M_PI;
//const double offsetQFem = M_PI;
//const double offsetQTib = 0;
//const double offsetY_foot2Com_left = -0.108;
//const double offsetY_foot2Com_right = 0.108;
//const double offsetX_Com = -0.01;

#define JACP_LEFT Lf*cos(qFem_left+thetaY)+Lt*cos(qFem_left+qTib_left+thetaY), Lt*cos(qFem_left+qTib_left+thetaY), \
-Lf*sin(qFem_left+thetaY)-Lt*sin(qFem_left+qTib_left+thetaY), -Lt*sin(qFem_left+qTib_left+thetaY)

#define JACP_RIGHT Lf*cos(qFem_right+thetaY)+Lt*cos(qFem_right+qTib_right+thetaY), Lt*cos(qFem_right+qTib_right+thetaY), \
-Lf*sin(qFem_right+thetaY)-Lt*sin(qFem_right+qTib_right+thetaY), -Lt*sin(qFem_right+qTib_right+thetaY)


Eigen::Matrix<double, 3, 1> posJoint2Base(double thetaX, double thetaY, Eigen::Matrix<double, 3, 1> qMotor){
    double q1 = thetaX + qMotor(0);
    double q2 = thetaY + qMotor(1);
    double q3 = q2 + qMotor(2);
   double x = Lf*sin(q2) + Lt * sin(q3);
    double zOrg = -Lf*cos(q2) - Lt * cos(q3);
    double y = -zOrg* sin(q1);
    double z = zOrg*cos(q1);
    Eigen::Matrix<double, 3, 1> res;
    res << -x, -y, -z;
    return res;
}

Eigen::Matrix<double, 3, 1> velJoint2Base(double thetaX, double thetaY, double dThetaX, double dThetaY, Eigen::Matrix<double, 3, 1> qMotor, Eigen::Matrix<double, 3, 1> dqMotor){
    double q1 = thetaX + qMotor(0);
    double q2 = thetaY + qMotor(1);
    double q3 = q2 + qMotor(2);
    Eigen::Matrix<double, 3, 3> JacobiJ2B;
    JacobiJ2B << 0,-(Lf*cos(q2)),-(Lt*cos(q3)),
                -(cos(q1)*(Lf*cos(q2) + Lt*cos(q3))),Lf*sin(q1)*sin(q2),Lt*sin(q1)*sin(q3),
                (-(Lf*cos(q2)) - Lt*cos(q3))*sin(q1),-(Lf*cos(q1)*sin(q2)),-(Lt*cos(q1)*sin(q3));

    double dq1 = dThetaX + dqMotor(0);
    double dq2 = dThetaY + dqMotor(1);
    double dq3 = dq2 + dqMotor(2);

    Eigen::Matrix<double, 3, 1> dq;
    dq << dq1, dq2, dq3;
    return JacobiJ2B*dq;
}

Eigen::Matrix<double, 3, 1> forceJoint2Foot(double thetaX, double thetaY, Eigen::Matrix<double, 3, 1> qMotor, Eigen::Matrix<double, 3, 1> forceMotor){
    Eigen::Matrix<double, 3, 3> JacobiJ2F;
    JacobiJ2F << 0,Lf*cos(qMotor(1) + thetaY) + Lt*cos(qMotor(1) + qMotor(2) + thetaY),Lt*cos(qMotor(1) + qMotor(2) + thetaY),
                cos(qMotor(0) + thetaX)*(Lf*cos(qMotor(1) + thetaY) + Lt*cos(qMotor(1) + qMotor(2) + thetaY)),sin(qMotor(0) + thetaX)*(-(Lf*sin(qMotor(1) + thetaY)) - Lt*sin(qMotor(1) + qMotor(2) + thetaY)),-(Lt*sin(qMotor(0) + thetaX)*sin(qMotor(1) + qMotor(2) + thetaY)),
                -((-(Lf*cos(qMotor(1) + thetaY)) - Lt*cos(qMotor(1) + qMotor(2) + thetaY))*sin(qMotor(0) + thetaX)),cos(qMotor(0) + thetaX)*(Lf*sin(qMotor(1) + thetaY) + Lt*sin(qMotor(1) + qMotor(2) + thetaY)),Lt*cos(qMotor(0) + thetaX)*sin(qMotor(1) + qMotor(2) + thetaY);
    return JacobiJ2F.transpose().inverse()*forceMotor;
}

Eigen::Matrix<double, 3, 1> posJoint2Foot(double thetaX, double thetaY, Eigen::Matrix<double, 3, 1> qMotor){
    double q1 = qMotor(0) + thetaX;
    double q2 = thetaY + qMotor(1);
    double q3 = q2 + qMotor(2);
    double x = Lf*sin(q2) + Lt * sin(q3);
    double zOrg = -Lf*cos(q2) - Lt * cos(q3);
    double y = -zOrg* sin(q1);
    double z = zOrg*cos(q1);
    Eigen::Matrix<double, 3, 1> res;
    res << x, y, z;
    return res;
}

Eigen::Matrix<double, 3, 1> posFoot2Joint(Eigen::Matrix<double, 3, 1> &posFoot, double thetaX, double thetaY){
    double desiredXSW, desiredYSW, desiredZSW;
    desiredXSW = posFoot(0);
    desiredYSW = posFoot(1);
    desiredZSW = posFoot(2);

    double lVirtualLeg = sqrt(desiredXSW*desiredXSW+desiredYSW*desiredYSW+desiredZSW*desiredZSW);
    double sigma, alpha, beta;
    sigma = pi - asin(desiredXSW/lVirtualLeg );
    alpha = acos((lVirtualLeg * lVirtualLeg + Lf * Lf - Lt * Lt) / (2 * Lf * lVirtualLeg));
    beta = acos((Lf*Lf +Lt*Lt -desiredXSW*desiredXSW - lVirtualLeg * lVirtualLeg) / (2 * Lt * Lf));
    if(isnan(alpha)){
        cout <<"isnan" <<endl;
    }
    double desiredQFemSW = sigma-alpha-thetaY;
    double desiredQTibSW = pi-beta;
    double desiredQAbAd = atan2(-desiredYSW, -desiredZSW) - thetaX + pi;
    Eigen::Matrix<double, 3, 1> desiredPos_motor_SW;
    desiredPos_motor_SW << desiredQAbAd, desiredQFemSW, desiredQTibSW;

    return desiredPos_motor_SW;
}

Eigen::Matrix<double, 3, 1> forceVirtual2Joint(Eigen::Matrix<double, 3, 1> &forceVirtual, double qHip, double qFem, double qTib){
    Eigen::Matrix<double, 3, 3> Jp;
    Eigen::Matrix<double, 3, 1> forceJoint;
    Jp << 1, 0, 0,
          0, 1, -(((Lt*(2*pow(Lf,2) + 2*Lf*Lt*cos(qTib))*sin(qTib))/(2.*pow(pow(Lf,2) + pow(Lt,2) + 2*Lf*Lt*cos(qTib),1.5)) - (Lt*sin(qTib))/sqrt(pow(Lf,2) + pow(Lt,2) + 2*Lf*Lt*cos(qTib)))/sqrt(1 - pow(2*pow(Lf,2) + 2*Lf*Lt*cos(qTib),2)/(4.*pow(Lf,2)*(pow(Lf,2) + pow(Lt,2) + 2*Lf*Lt*cos(qTib))))),
          0, 0, -((Lf*Lt*sin(qTib))/sqrt(pow(Lf,2) + pow(Lt,2) + 2*Lf*Lt*cos(qTib)));
//    cout << "orig forceVirtual: " << forceVirtual.transpose() << endl;
//    cout << "orig Jp: \n" << Jp << endl;
    forceJoint = Jp.transpose()*forceVirtual;
    return forceJoint;
}

vector<double> calculateSSP(double y0, double dy0, double z0, double t){
    if (t<0)
        t=0;
    double lambda = sqrt(g/z0);
    double c1 = 1./2*(y0+dy0/lambda);
    double c2 = 1./2*(y0-dy0/lambda);
    double ySSP = c1* pow(M_E,lambda*t)+c2* pow(M_E,-lambda*t);
    double dySSP = lambda*(c1* pow(M_E,lambda*t)-c2* pow(M_E,-lambda*t));
    return vector<double>{ySSP, dySSP};
}

//-------------------------------------- init_state --------------------------------------------
void init_controller(const mjModel* m, mjData* d)
{
    const int qNum = 12;
    // qTorso = 0.5
//    double initQ[qNum] = {
//            0,
//            0,
//            0,
//            -0.053135,
//            0.5,
//            0,
//            0,
//            1.959312-offsetQFem,
//            1.365517-offsetQTib,
//            0,
//            2.042923-offsetQFem,
//            1.240025-offsetQTib,
//    };
    double initQ[qNum] = {
            0,
            0,
            0,
            0,
            0.,
            0,
            0,
            2.570448-offsetQFem,
            1.241099-offsetQTib,
            0,
            2.600296-offsetQFem,
            1.186003-offsetQTib,
    };
    for (int i = 0; i < qNum; ++i) {
        d->qpos[i] = initQ[i];
    }

    double initDq[qNum] = {
//            -0.005506226,
            -0.056098671,
//            0.008394263,
            -0.56028458,
            1.041142356,
            -0.599870728,
            1.099170618
    };
    for (int i = 0; i < qNum; ++i) {
//        d->qvel[i] = initDq[i];
    }
    d->qvel[0] = 0.;
    d->qvel[1] = 0.385;
    //mj_forward(m,d);
}

//------------------------------------ controllerHandler -------------------------------------------
void controllerHandler(const mjModel *m, mjData *d){
    const double Kp = 2000;
    const double Kd = 50;
    const double Kp_Torso = 200;
    const double Kd_Torso = 10;
    const double Kp_motor = 100;
    const double Kv_motor = 1;
    const double duration = 0.3;
    const double desiredVelocity = 0.;
    const double desiredWs = 0.108*2;
    const double desiredZ0 = 0.35;
    const double Kp_SWY = 0.3;
    const int numPeriod = 1;
    const int samplingRate = 1;
    const int numStep = 20;

    static bool init = false;
    static bool hasPlot = false;

    static bool isLeftStand = true;
    static bool hasFirstPlaned = false;
    static bool hasPlaned[numPeriod] = {false};

    static double stepTime = 0.;
    static double index = 0;
    static double sumTime = 0;

    static double lastTime=0;
    static double time = 0;
    time=d->time;
    double deltaTime = time-lastTime;
    lastTime = time;

    double period[numPeriod] = {0.6};

    if(d->time!=0){
        stepTime+=0.0005;
//        t+=0.0005;
    }else{
        return;
    }
    double tau = stepTime/duration;
    double averageDuration = 0.3;//index==0?duration:(sumTime/index);

    // 计算base方向角
    int sensorIndex_framequat= 0;
    double base_state_quat_w = d->sensordata[sensorIndex_framequat+0];
    double base_state_quat_x = d->sensordata[sensorIndex_framequat+1];
    double base_state_quat_y = d->sensordata[sensorIndex_framequat+2];
    double base_state_quat_z = d->sensordata[sensorIndex_framequat+3];
    Eigen::Quaternion quat = Eigen::Quaternion(base_state_quat_w,base_state_quat_x,base_state_quat_y,base_state_quat_z);
    double thetaY = asin(2 * base_state_quat_w * base_state_quat_y - base_state_quat_z * base_state_quat_x);
    double thetaX = atan2(2*(base_state_quat_w*base_state_quat_x+base_state_quat_y*base_state_quat_z), 1-2*(base_state_quat_x*base_state_quat_x+base_state_quat_y*base_state_quat_y));

    // 获取质心世界坐标系角速度
    int sensorIndex_gyro = 4;
    double dThetaX = d->sensordata[sensorIndex_gyro+0];
    double dThetaY = d->sensordata[sensorIndex_gyro + 1];


    // 获取电机的反馈位置
    int sensorIndex_actuatorPos = 7;
    double qAbAd_left = d->sensordata[sensorIndex_actuatorPos+0] + offsetQAbAd;
    double qFem_left = d->sensordata[sensorIndex_actuatorPos+1] + offsetQFem;
    double qTib_left = d->sensordata[sensorIndex_actuatorPos+2] + offsetQTib;
    double qAbAd_right = d->sensordata[sensorIndex_actuatorPos+3] + offsetQAbAd;
    double qFem_right = d->sensordata[sensorIndex_actuatorPos+4] + offsetQFem;
    double qTib_right = d->sensordata[sensorIndex_actuatorPos+5] + offsetQTib;
    Eigen::Matrix<double, 3, 1> q_motor_left, q_motor_right;
    q_motor_left << qAbAd_left, qFem_left, qTib_left;
    q_motor_right << qAbAd_right, qFem_right, qTib_right;

    // 获取电机的反馈速度
    int sensorIndex_actuatorVel = 13;
    double dqAbAd_left = d->sensordata[sensorIndex_actuatorVel+0];
    double dqFem_left = d->sensordata[sensorIndex_actuatorVel+1];
    double dqTib_left = d->sensordata[sensorIndex_actuatorVel+2];
    double dqAbAd_right = d->sensordata[sensorIndex_actuatorVel+3];
    double dqFem_right = d->sensordata[sensorIndex_actuatorVel+4];
    double dqTib_right = d->sensordata[sensorIndex_actuatorVel+5];
    Eigen::Matrix<double, 3, 1> dq_motor_left, dq_motor_right;
    dq_motor_left << dqAbAd_left, dqFem_left, dqTib_left;
    dq_motor_right << dqAbAd_right, dqFem_right, dqTib_right;

    // 读取电机的反馈扭矩
    int sensorIndex_actuatorFrc = 19;
    double force_AbAd_left = d->sensordata[sensorIndex_actuatorFrc+0];
    double force_Fem_left = d->sensordata[sensorIndex_actuatorFrc+1];
    double force_Tib_left = d->sensordata[sensorIndex_actuatorFrc+2];
    double force_AbAd_right = d->sensordata[sensorIndex_actuatorFrc+3];
    double force_Fem_right = d->sensordata[sensorIndex_actuatorFrc+4];
    double force_Tib_right = d->sensordata[sensorIndex_actuatorFrc+5];
    Eigen::Matrix<double, 3, 1> force_motor_left, force_motor_right;
    force_motor_left << d->sensordata[sensorIndex_actuatorFrc + 0], d->sensordata[sensorIndex_actuatorFrc + 1], d->sensordata[sensorIndex_actuatorFrc + 2];
    force_motor_right << d->sensordata[sensorIndex_actuatorFrc + 3], d->sensordata[sensorIndex_actuatorFrc + 4], d->sensordata[sensorIndex_actuatorFrc + 5];

    // 计算足端反力
    Eigen::Matrix<double, 3, 1> foot_feedback_force_l, foot_feedback_force_r;
    foot_feedback_force_r = forceJoint2Foot(thetaX, thetaY, q_motor_right, force_motor_right);
    foot_feedback_force_l = forceJoint2Foot(thetaX, thetaY, q_motor_left, force_motor_left);
    double forceGround2leftFootZ = -foot_feedback_force_l[2];
    double forceGround2rightFootZ = -foot_feedback_force_r[2];



    // 获取Base相对站立腿足端的位置
    Eigen::Matrix<double, 3, 1> base_state_linear_pos;
    if(isLeftStand){
        base_state_linear_pos = posJoint2Base(thetaX, thetaY, q_motor_left);
    }else{
        base_state_linear_pos = posJoint2Base(thetaX, thetaY, q_motor_right);
    }

    // 获取Base相对站立腿足端速度
    Eigen::Matrix<double, 3, 1> base_state_linear_vel;
    Eigen::Matrix<double, 3, 1> test_base_state_linear_vel;
    if(isLeftStand){
        base_state_linear_vel = velJoint2Base(thetaX, thetaY, dThetaX, dThetaY, q_motor_left, dq_motor_left);
    }else{
        base_state_linear_vel = velJoint2Base(thetaX, thetaY, dThetaX, dThetaY, q_motor_right, dq_motor_right);
    }


    // 获取摆动腿足端相对Base的位置
    Eigen::Matrix<double, 3, 1> pos_foot_left, pos_foot_right;
    pos_foot_right = posJoint2Foot(thetaX, thetaY, q_motor_right);
    pos_foot_left = posJoint2Foot(thetaX, thetaY, q_motor_left);
    double px_foot2Com_left = pos_foot_left[0];
    double px_foot2Com_right = pos_foot_right[0];
    double py_foot2Com_left = pos_foot_left[1];
    double py_foot2Com_right = pos_foot_right[1];
    double pz_foot2Com_left = pos_foot_left[2];
    double pz_foot2Com_right = pos_foot_right[2];

    cout << "stepTime: " << stepTime << endl;
    cout << "orig tau: " << tau << endl;
    if(isLeftStand)
        cout << "standing foot is left!\n" << endl;
    else
        cout << "standing foot is right!\n" << endl;

    printf("thetaX: %f, thetaY: %f\n", thetaX, thetaY);
    printf("dThetaX: %f, dThetaY: %f\n", dThetaX, dThetaY);

    cout << "qMotorLeft: " << q_motor_left.transpose() << endl;
    cout << "qMotorRight: " << q_motor_right.transpose() << endl;
    cout << "dqMotorLeft: " << dq_motor_left.transpose() << endl;
    cout << "dqMotorRight: " << dq_motor_right.transpose() << endl;
    cout << "forceMotorLeft: " << force_motor_left.transpose() << endl;
    cout << "forceMotorRight: " << force_motor_right.transpose() << endl;

    cout << "base_state_linear_pos: " << base_state_linear_pos.transpose() << endl;
    cout << "base_state_linear_vel: " << base_state_linear_vel.transpose() << endl;

    if(isLeftStand){
        cout << "swing foot pos: " << pos_foot_right.transpose() << endl;
        printf("leftFootFeedbackForce: %f \n", forceGround2leftFootZ);
    }else{
        cout << "swing foot pos: " << pos_foot_left.transpose() << endl;
        printf("rightFootFeedbackForce: %f\n", forceGround2rightFootZ);
    }


    // -----------------------------------------判断触地腿-------------------------------------------------------------------
    static int count = 0;
    if(isLeftStand){
        if(forceGround2rightFootZ > 20 && tau > 0.5)
            count++;
        else
            count = 0;
        if(count>10||stepTime>duration){
            for (int i = 0; i < numPeriod; ++i) {
                hasPlaned[i] = false;
            }
            hasFirstPlaned = false;
            count = 0;
            isLeftStand = false;
            tau = 0;
            sumTime+=stepTime;
            stepTime = 0;
            index++;
        }
    }else{
        if(forceGround2leftFootZ > 20 && tau > 0.5)
            count++;
        else
            count = 0;
        if(count>10||stepTime>duration){
            for (int i = 0; i < 4; ++i) {
                hasPlaned[i] = false;
            }
            hasFirstPlaned = false;
            count = 0;
            isLeftStand = true;
            tau = 0;
            sumTime+=stepTime;
            stepTime = 0;
            index++;
        }
    }

    double target_base_state_pos[5] = {0, 0, desiredZ0, 0, 0.};
    double target_base_state_vel[5] = {desiredVelocity, 0, 0, 0, 0};
    // 质心的状态
    double base_state_pos[5] = {base_state_linear_pos[0], base_state_linear_pos[1], base_state_linear_pos[2], thetaX, thetaY};
    double base_state_vel[5] = {base_state_linear_vel[0], base_state_linear_vel[1], base_state_linear_vel[2], dThetaX, dThetaY};
//    printf("base_state_pos: %f %f %f %f %f\n", base_state_pos[0], base_state_pos[1], base_state_pos[2], base_state_pos[3], base_state_pos[4]);
//    printf("base_state_vel: %f %f %f %f %f\n", base_state_vel[0], base_state_vel[1], base_state_vel[2], base_state_vel[3], base_state_vel[4]);

    // 质心当前状态与期望的偏差
    Eigen::Matrix<double, 5, 1> dp_base_state_pos, dp_base_state_vel;
    for (int i = 0; i < 5; ++i) {
        dp_base_state_pos(i, 0) = target_base_state_pos[i] - base_state_pos[i];
        dp_base_state_vel(i, 0) = target_base_state_vel[i] - base_state_vel[i];
    }
//    printf("dp_base_state_pos: %f %f %f %f %f\n", dp_base_state_pos[0], dp_base_state_pos[1], dp_base_state_pos[2], dp_base_state_pos[3], dp_base_state_pos[4]);
//    printf("dp_base_state_vel: %f %f %f %f %f\n", dp_base_state_vel[0], dp_base_state_vel[1], dp_base_state_vel[2], dp_base_state_vel[3], dp_base_state_vel[4]);


    // 站立腿控制
    Eigen::Matrix<double, 3, 1> motorST;
    Eigen::Matrix<double, 3, 1>  forceVirtual;
    if (isLeftStand){
        forceVirtual(0) = (-Kp_Torso*dp_base_state_pos[3]-Kd_Torso*dp_base_state_vel[3] );
        forceVirtual(1) = mT*LT*sin(thetaY)-(Kp_Torso*dp_base_state_pos[4]+Kd_Torso*dp_base_state_vel[4]);
        forceVirtual(2) =  mr*g/base_state_pos[2]* sqrt(px_foot2Com_left*px_foot2Com_left+base_state_pos[2]*base_state_pos[2]+py_foot2Com_left*py_foot2Com_left) + Kp*dp_base_state_pos[2] + Kd*dp_base_state_vel[2];
        motorST = forceVirtual2Joint(forceVirtual, qAbAd_left, qFem_left, qTib_left);
//        d->ctrl[0] = 100*(pi-qAbAd_left)+1*(0-dqAbAd_left);
        d->ctrl[0] = motorST(0);
        d->ctrl[1] = motorST(1);
        d->ctrl[2] = motorST(2);
    }else{
        forceVirtual(0) = (-Kp_Torso*dp_base_state_pos[3]-Kd_Torso*dp_base_state_vel[3] );
        forceVirtual(1) = mT*LT*sin(thetaY)-(Kp_Torso*dp_base_state_pos[4]+Kd_Torso*dp_base_state_vel[4] );
        forceVirtual(2) = mr*g/base_state_pos[2]* sqrt(px_foot2Com_right*px_foot2Com_right+base_state_pos[2]*base_state_pos[2]) + Kp*dp_base_state_pos[2] + Kd*dp_base_state_vel[2];
        motorST = forceVirtual2Joint(forceVirtual, qAbAd_right, qFem_right, qTib_right);
//        d->ctrl[3] = 100*(pi-qAbAd_right)+1*(0-dqAbAd_right);
        d->ctrl[3] = motorST(0);
        d->ctrl[4] = motorST(1);
        d->ctrl[5] = motorST(2);
    }


    //-----------------------------------------摆动腿控制------------------------------------------
    /// 计算估计的SSP临界状态y，dy
    double lambda = sqrt(g/base_state_pos[2]);
    double sigma = lambda* sinh(averageDuration/2*lambda)/ cosh(averageDuration/2*lambda);
    double desiredDySSP = (desiredWs*sigma)/2;
    double y0, dy0;
    vector<double> stateSSP;
    if(isLeftStand){
        y0 = py_foot2Com_left - offsetY_foot2Com_left;
        dy0 = -base_state_vel[1];
        stateSSP = calculateSSP(y0, dy0, base_state_pos[2], averageDuration-tau*duration);
    }else{
        y0 = -py_foot2Com_right + offsetY_foot2Com_right;
        dy0 = base_state_vel[1];
        stateSSP = calculateSSP(y0, dy0, base_state_pos[2], averageDuration-tau*duration);
    }


    static double startXFoot, startYFoot, startZFoot;
    static double pxSWFoot, pySWFoot, pzSWFoot, newPxSWFoot, newPySWFoot, newPzSWFoot;
    static double nextX, nextY, nextZ, k_yx, k_zx;
    double vxBase, vxDes, pzBase;
    Eigen::Matrix<double, 3, 1> beginPoint, endPoint, controlPoint1, controlPoint2;
    static Eigen::Matrix<double, 3, 1> desiredSWFootPos, nextDesiredSWFootPos, dNextDesiredSWFootPos;
    if(tau<period[0]){
        if(!hasFirstPlaned) {
            hasFirstPlaned = true;
            if (isLeftStand) {
                pxSWFoot = px_foot2Com_right;
                pySWFoot = py_foot2Com_right;
                pzSWFoot = pz_foot2Com_right;

                newPySWFoot = -stateSSP[0]+Kp_SWY*(stateSSP[1] - desiredDySSP) + offsetY_foot2Com_right;
            } else {
                pxSWFoot = px_foot2Com_left;
                pySWFoot = py_foot2Com_left;
                pzSWFoot = pz_foot2Com_left;
                newPySWFoot = stateSSP[0] - Kp_SWY*(stateSSP[1] - desiredDySSP) + offsetY_foot2Com_left;
            }
            vxBase = base_state_linear_vel[0];
            vxDes = target_base_state_vel[0];
            pzBase = base_state_pos[2];
            if (vxBase * vxDes < 0)
                vxDes = 0;
            if (vxBase * vxBase < vxDes * vxDes)
                newPxSWFoot = vxDes > 0 ? -0.02 : 0.02;
            else if (vxBase < 0.)
                newPxSWFoot = -sqrt((vxBase * vxBase - vxDes * vxDes) * pzBase / g);
            else
                newPxSWFoot = sqrt((vxBase * vxBase - vxDes * vxDes) * pzBase / g);
            newPzSWFoot = -target_base_state_pos[2];
            newPxSWFoot += offsetX_Com;
        }
        beginPoint << pxSWFoot, pySWFoot, pzSWFoot;
        endPoint << newPxSWFoot, newPySWFoot, newPzSWFoot;
        controlPoint1 << pxSWFoot, pySWFoot, pzSWFoot+0.1;
        controlPoint2 << newPxSWFoot, newPySWFoot, newPzSWFoot+0.1;
        desiredSWFootPos = bezierPoint3D(beginPoint, controlPoint1, controlPoint2, endPoint, tau);
        nextDesiredSWFootPos = bezierPoint3D(beginPoint, controlPoint1, controlPoint2, endPoint, tau+0.0005/duration);
        k_yx = (nextDesiredSWFootPos(1)-desiredSWFootPos(1))/(nextDesiredSWFootPos(0)-desiredSWFootPos(0));
        k_zx = (nextDesiredSWFootPos(2)-desiredSWFootPos(2))/(nextDesiredSWFootPos(0)-desiredSWFootPos(0));
        printf("start point: %f %f %f\n", beginPoint(0), beginPoint(1), beginPoint(2));
        printf("end point: %f %f %f\n", newPxSWFoot, newPySWFoot, newPzSWFoot);
        cout << tau << " time, desired point: " << desiredSWFootPos.transpose() << endl;
    }
    for (int i = 0; i < numPeriod; ++i) {
        double time = period[i];
        if(tau>=time){
            static double px1, py1, pz1, px2, py2, pz2, rate;
            if(!hasPlaned[i]){
                hasPlaned[i] = true;
                if (isLeftStand) {
                    pxSWFoot = px_foot2Com_right;
                    pySWFoot = py_foot2Com_right;
                    pzSWFoot = pz_foot2Com_right;
                    newPySWFoot = -stateSSP[0]+Kp_SWY*(stateSSP[1] - desiredDySSP) + offsetY_foot2Com_right;
                } else {
                    pxSWFoot = px_foot2Com_left;
                    pySWFoot = py_foot2Com_left;
                    pzSWFoot = pz_foot2Com_left;
                    newPySWFoot = stateSSP[0] - Kp_SWY*(stateSSP[1] - desiredDySSP) + offsetY_foot2Com_left;
                }
                vxBase = base_state_linear_vel[0];
                vxDes = target_base_state_vel[0];
                pzBase = base_state_pos[2];
                if(vxBase*vxDes<0)
                    vxDes = 0;
                if(vxBase*vxBase<vxDes*vxDes)
                    newPxSWFoot = vxDes>0?-0.02:0.02;
                else if (vxBase < 0.)
                    newPxSWFoot = -sqrt((vxBase * vxBase - vxDes * vxDes) * pzBase / g);
                else
                    newPxSWFoot = sqrt((vxBase * vxBase - vxDes * vxDes) * pzBase / g);
                newPzSWFoot = -target_base_state_pos[2];
                newPxSWFoot += offsetX_Com;
            }
            Eigen::Matrix<double, 3, 1> beginPoint, endPoint, controlPoint1, controlPoint2;
            beginPoint << nextDesiredSWFootPos(0), nextDesiredSWFootPos(1), nextDesiredSWFootPos(2);
            endPoint << newPxSWFoot, newPySWFoot, newPzSWFoot;
            double dx = endPoint(0) - beginPoint(0);
            double dz = endPoint(2) - beginPoint(2);
            controlPoint1 << beginPoint(0)+dx/4, beginPoint(1)+dx/4*k_yx, beginPoint(2)+dx/4*k_zx;
            controlPoint2 << endPoint(0), endPoint(1), endPoint(2)+dz/4;
//            desiredSWFootPos = bezierPoint3D(beginPoint, controlPoint1, controlPoint2, endPoint, (tau - time) / (1 - time));
            desiredSWFootPos = bezierPoint3D(beginPoint, controlPoint2, endPoint, (tau - time) / (1 - time));
            printf("start point: %f %f %f\n", beginPoint(0), beginPoint(1), beginPoint(2));
            printf("end point: %f %f %f\n", newPxSWFoot, newPySWFoot, newPzSWFoot);
            cout << tau << " time, desired point: " << desiredSWFootPos.transpose() << endl;
        }
    }

    Eigen::Matrix<double, 3, 1> desPos_motor_SW = posFoot2Joint(desiredSWFootPos, thetaX, thetaY);
    double desiredQAbAdSW, desiredQFemSW, desiredQTibSW;
    desiredQAbAdSW = desPos_motor_SW[0];
    desiredQFemSW = desPos_motor_SW[1];
    desiredQTibSW = desPos_motor_SW[2];
    Eigen::Matrix<double, 3, 1> motorSW;
    if(isLeftStand){
        if(isnan(desiredQFemSW)||isnan(desiredQTibSW)|| isnan(desiredQAbAdSW)){
            desiredQAbAdSW = qAbAd_right;
            desiredQFemSW = qFem_right;
            desiredQTibSW = qTib_right;
        }
        motorSW << Kp_motor*(desiredQAbAdSW-qAbAd_right)+(Kv_motor+10)*(0-dqAbAd_right),
                Kp_motor*(desiredQFemSW-qFem_right)+Kv_motor*(0-dqFem_right),
                Kp_motor*(desiredQTibSW-qTib_right)+Kv_motor*(0-dqTib_right);

        d->ctrl[3] = Kp_motor*(desiredQAbAdSW-qAbAd_right)+(Kv_motor+10)*(0-dqAbAd_right);
        d->ctrl[4] = Kp_motor*(desiredQFemSW-qFem_right)+Kv_motor*(0-dqFem_right);
        d->ctrl[5] = Kp_motor*(desiredQTibSW-qTib_right)+Kv_motor*(0-dqTib_right);
    }else{
        if(isnan(desiredQFemSW)||isnan(desiredQTibSW)|| isnan(desiredQAbAdSW)){
            desiredQAbAdSW = qAbAd_left;
            desiredQFemSW = qFem_left;
            desiredQTibSW = qTib_left;
        }
        motorSW << Kp_motor*(desiredQAbAdSW-qAbAd_left)+(Kv_motor+10)*(0-dqAbAd_left),
                Kp_motor*(desiredQFemSW-qFem_left)+Kv_motor*(0-dqFem_left),
                Kp_motor*(desiredQTibSW-qTib_left)+Kv_motor*(0-dqTib_left);

        d->ctrl[0] = Kp_motor*(desiredQAbAdSW-qAbAd_left)+(Kv_motor+10)*(0-dqAbAd_left);
        d->ctrl[1] = Kp_motor*(desiredQFemSW-qFem_left)+Kv_motor*(0-dqFem_left);
        d->ctrl[2] = Kp_motor*(desiredQTibSW-qTib_left)+Kv_motor*(0-dqTib_left);
    }
    cout << "motorST: " << motorST.transpose() << endl;
    cout << "motorSW: " << motorSW.transpose() << endl << endl << endl;

//    cout << "\n\norig lambda, sigma: " << lambda << " " << sigma << endl;
//    cout << "orig averageDuration: " << averageDuration << endl;
//    cout << "orig averageDuration: " << averageDuration << endl;
    cout << "orig desiredDySSP" << desiredDySSP << endl;
    cout << "orig StateSSP: " << stateSSP[0] << " " << stateSSP[1] << endl;
    cout << "orig desiredSWFootPos: " << desiredSWFootPos.transpose() << endl;
    cout << "orig desiredSWMotorPos： " << desPos_motor_SW.transpose() << endl<<endl;



    static vector<double> vector_time, vector_leftFootForceZ;
    static double t = 0.;
    static int countPlot = 0;
    if(countPlot%samplingRate==0){
        countPlot = 0;
        vector_time.push_back(t);
        vector_leftFootForceZ.push_back(foot_feedback_force_l[2]);
    }
    countPlot++;
    if(index==numStep&&!hasPlot){
        hasPlot = true;
        plt::named_plot("ground force 2 left foot", vector_time, vector_leftFootForceZ);
        plt::show();
    }
    t+=0.0005;


    vector<double> floatMessage;
    vector<double> motorMessage;
    vector<double> control;
    floatMessage = {thetaX, thetaY, 0, dThetaX, dThetaY, 0};
    motorMessage = {qAbAd_left, qFem_left, qTib_left, qAbAd_right, qFem_right, qTib_right,
                    dqAbAd_left, dqFem_left, dqTib_left, dqAbAd_right, dqFem_right, dqTib_right,
                    force_AbAd_left, force_Fem_left, force_Tib_left, force_AbAd_right, force_Fem_right, force_Tib_right};

    static Robbie robbie(LEFTLEG, motorMessage, floatMessage);
    if(d->time!=0)
        robbie.update(motorMessage, floatMessage, 0.0005);
    robbie.printRobotMessage();
}

#endif //MUJOCOTEST2D_SIMPLEMODELCONTROLLER_H