//
// Created by yeyinong on 2023/3/29.
//

#ifndef MUJOCOTEST3D_ROBBIE_H
#define MUJOCOTEST3D_ROBBIE_H

#include <vector>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "iostream"
#include "tools/tools.h"
#include "fstream"
#include "parameters/parameters.h"

using namespace std;



class Robbie {
private:
    /// command
    bool hasInit = false;
    bool cmdUpdate_ = false;
    double cmdVelocity_ = 0;
    double cmdVelocityBuf_ = 0;
    double cmdHeight_ = 0.1;
//    double cmdDirection = 1.57;
    double cmdDirection = 0;

    double tau_ = 0;
    double t_ = 0;
    double sumTime_ = 0;
    double averageStepTime = stepTime;

    bool standingLeg_ = LEFTLEG;
    bool isJustChangeLeg = true;

    /// state of robot
    Eigen::Matrix<double, 3, 1> a_robBIE_, a_world_;
    Eigen::Matrix<double, 3, 1> pos_world, vel_world;
    Eigen::Matrix<double, 3, 1> posFloatBase_world;
    Eigen::Matrix<double, 3, 1> qFloatBase_;
    Eigen::Matrix<double, 3, 1> dqFloatBase_;
    Eigen::Matrix<double, 3, 3> rotationMatrix_;

    Eigen::Matrix<double, 4, 1> qMotorLeft_;
    Eigen::Matrix<double, 4, 1> qMotorRight_;
    Eigen::Matrix<double, 4, 1> dqMotorLeft_;
    Eigen::Matrix<double, 4, 1> dqMotorRight_;
    Eigen::Matrix<double, 4, 1> torqueMotorLeft_;
    Eigen::Matrix<double, 4, 1> torqueMotorRight_;

    /// state of robot which need to be calculated
    Eigen::Matrix<double, 4, 4> JacobiJ2F_left_, JacobiJ2F_right_;

    Eigen::Matrix<double, 4, 1> forceGroundLeft_;
    Eigen::Matrix<double, 4, 1> forceGroundLeftFil_;
    Eigen::Matrix<double, 4, 1> forceGroundRight_;
    Eigen::Matrix<double, 4, 1> forceGroundRightFil_;

//    Eigen::Matrix<double, 3, 1> posFloatBase_;
    Eigen::Matrix<double, 3, 1> posFloatBase_RobBIE_;

    Eigen::Matrix<double, 3, 1> posBase2ST_RobBIE_;
    Eigen::Matrix<double, 3, 1> posBase2ST_world_;
    Eigen::Matrix<double, 3, 1> posFloatBase_world_;

    Eigen::Matrix<double, 3, 1> posLeftFoot2Base_RobBIE_;
    Eigen::Matrix<double, 3, 1> posRightFoot2Base_RobBIE_;
    Eigen::Matrix<double, 3, 1> posLeftFoot2Base_world_;
    Eigen::Matrix<double, 3, 1> posRightFoot2Base_world_;

    Eigen::Matrix<double, 3, 1> posSW2Base_RobBIE_;
    Eigen::Matrix<double, 3, 1> posSW2Base_world_;


    Eigen::Matrix<double, 3, 1> velFloatBase_;
    Eigen::Matrix<double, 3, 1> velFloatBaseFil_;

    Eigen::Matrix<double, 3, 1> posSW2Base_world_begein_;
    Eigen::Matrix<double, 3, 1> posSWingFoot_RobBIE_;
    double yawSWFoot_RobBIE_begin_ = 0;
    double yawSWFoot_world_begin_ = 0;
    double yawSTFoot_RobBIE_begin_ = 0;
    double yawSTFoot_world_begin_ = 0;

    Eigen::Matrix<double, 4, 1> desQSwingMotor_;
    Eigen::Matrix<double, 3, 1> desPosSwingFoot_;
    Eigen::Matrix<double, 3, 1> desPosSwingFoot_RobBIE_;

    Eigen::Matrix<double, 4, 1> desTorqueLeftMotor_;
    Eigen::Matrix<double, 4, 1> desTorqueRightMotor_;

    // sspY
    Eigen::Matrix<double, 2, 1> estimatedStateSSP_;
    double estimatedDesDyStateSSP_;
    double y0_, dy0_;
    double lambda_;
    double sigma_;

    /// state HZD
    vector<Eigen::Matrix<double , 4, 6>> alphaHZD_;
    vector<double> velHZD_;
    double s;
    double s_fast;
    double s_slow;
    double SWHipPosDes;
    double SWKneePosDes;
    double STHipVelDes;
    double STKneeVelDes;
    double SWHipVelDes;
    double SWKneeVelDes;
    double dpSWAbd_Begin;
    double dpSWHip_Begin;
    double dpSWKnee_Begin;
    double dpSWFootX_Begin_;
    double dpSWFootY_Begin_;
    double dpSWFootZ_Begin_;

    /// Kalman Filter
    Eigen::Matrix<double,3,1> vOpOpswT_KF_x;
    Eigen::Matrix<double,3,1> vOpOpswT_KF_y;
    Eigen::Matrix<double,3,1> vOpOpswT_KF_z;

    Eigen::Matrix<double,3,3> sigma_vOpOpswT_KF_x;
    Eigen::Matrix<double,3,3> sigma_vOpOpswT_KF_y;
    Eigen::Matrix<double,3,3> sigma_vOpOpswT_KF_z;

    Eigen::Matrix<double, 3, 3> Cov_LRT_r;
    Eigen::Matrix<double, 3, 3> Cov_LinearAccelerator;

    /// state of robot which update per step
    double countStep_ = 0;

    ofstream file_;

    void forward_kinematics_right(const Eigen::Matrix<double, 4, 1>& qMotor, Eigen::Matrix<double, 4, 1>& pFoot);
    void forward_kinematics_left(const Eigen::Matrix<double, 4, 1>& qMotor, Eigen::Matrix<double, 4, 1>& pFoot);

    void jacobianLeft(const Eigen::Matrix<double, 4, 1>& qMotor, Eigen::Matrix<double, 4, 4>& Jacobian);
    void jacobianRight(const Eigen::Matrix<double, 4, 1>& qMotor, Eigen::Matrix<double, 4, 4>& Jacobian);

    void updateMotorMessage(vector<double>& motorMessage);
    void updateFloatBaseMessage(vector<double>& doubleBaseMessage);
    void updateJacLeft();
    void updateJacRight();

    void calculateGroundForce();
    void updateStandingFoot();

    void calculatePosLeftFoot2Base_RobBIE();
    void calculatePosRightFoot2Base_RobBIE();

    void calculatePosFloatBase();
    void calculateVelFloatBase();
    void calculatePosSwingFoot();

    void controlStandingLeg();

    void calculateStateSSP();
    void calculateDesPosSwingFoot();
    void calculateDesQSwingMotor();
    void controlSwingLeg();

    void estimator(double deltaTime);

    /// test
    void controlTest();
public:
    Robbie();
    Robbie(bool standingLeg, vector<double>& motorMessage, vector<double>& doubleBaseMessage);
    void update(vector<double>& motorMessage, vector<double>& doubleBaseMessage, double deltaTime);
    vector<double> controlMotor(vector<double>& control);
    void logData(vector<double> &vectorData);
    void getVelX(double* vel, double* velFil);
};

void q2Foot_Left(Eigen::Matrix<double, 3, 1>& leftFoot, const Eigen::Matrix<double, 3, 1>& qLeft);

void q2Base_Left(Eigen::Matrix<double, 3, 1>& leftBase, const Eigen::Matrix<double, 3, 1>& qLeft);

void q2Foot_Right(Eigen::Matrix<double, 3, 1>& rightFoot, const Eigen::Matrix<double, 3, 1>& qRight);

void q2Base_Right(Eigen::Matrix<double, 3, 1>& rightBase, const Eigen::Matrix<double, 3, 1>& qRight);

void foot2Q_Left(Eigen::Matrix<double, 3, 1>& qLeft, const Eigen::Matrix<double, 3, 1>& footLeft);

void foot2Q_Right(Eigen::Matrix<double, 3, 1>& qRight, const Eigen::Matrix<double, 3, 1>& footRight);

void calculateJb_Left(Eigen::Matrix<double, 3, 3>& JbLeft,const Eigen::Matrix<double, 3, 1>& qLeft);

void calculateJb_Right(Eigen::Matrix<double, 3, 3>& JbRight,const Eigen::Matrix<double, 3, 1>& qRight);
#endif //MUJOCOTEST3D_ROBBIE_H
