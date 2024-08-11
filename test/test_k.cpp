//
// Created by yeyinong on 2024/3/8.
//
#include "parameters/parameters.h"
#include "math.h"
#include "robbie/Robbie.h"
#include "iostream"

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
void forward_kinematics_right(const Eigen::Matrix<double, 4, 1>& qMotor, Eigen::Matrix<double, 4, 1>& pFoot){
    double roll, pitch, yaw, q1, q2, q3,q4;
    double dx1, dy1, dz1, dx2, dy2, dz2 , dx3, dy3, dz3, dx4, dy4, dz4, dx5, dy5, dz5;



    q1 = qMotor(0);
    q2 = qMotor(1);
    q3 = qMotor(2);
    q4 = qMotor(3);

//    cout << roll << " " << pitch << " " << yaw << endl;
//    cout << q1 << " " << q2 << " " << q3 << " " << q4 <<endl;

    dx1 = dx1_R;
    dy1 = dy1_R;
    dz1 = dz1_R;

    dx2 = dx2_R;
    dy2 = dy2_R;
    dz2 = dz2_R;

    dx3 = dx3_R;
    dy3 = dy3_R;
    dz3 = dz3_R;

    dx4 = dx4_R;
    dy4 = dy4_R;
    dz4 = dz4_R;

    dx5 = dx5_R;
    dy5 = dy5_R;
    dz5 = dz5_R;
    Eigen::Matrix<double, 4, 4> temp;
    /*temp << cos(q3 - q4)*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) + (cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3))*sin(q3 - q4),-(cos(q2)*sin(q1)),
            cos(q3 - q4)*(cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3)) - (cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3))*sin(q3 - q4),
            dx1 + dx2*cos(q1) + dx3*cos(q1) - dy2*sin(q1) - dy3*cos(q2)*sin(q1) - dy4*cos(q2)*sin(q1) - dy5*cos(q2)*sin(q1) + dz3*sin(q1)*sin(q2) + dz4*(cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3)) + dx4*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) +
            dx5*(cos(q3 - q4)*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) + (cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3))*sin(q3 - q4)) + dz5*(cos(q3 - q4)*(cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3)) - (cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3))*sin(q3 - q4)),
            cos(q3 - q4)*(cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3)) + (-(cos(q2)*cos(q3)*sin(q2)) + sin(q1)*sin(q3))*sin(q3 - q4),pow(cos(q2),2),
            cos(q3 - q4)*(-(cos(q2)*cos(q3)*sin(q2)) + sin(q1)*sin(q3)) - (cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3))*sin(q3 - q4),
            dy1 + dy2*cos(q2) + dy3*pow(cos(q2),2) + dy4*pow(cos(q2),2) + dy5*pow(cos(q2),2) + dx2*sin(q1) + dx3*sin(q1) - dz3*cos(q2)*sin(q2) + dz4*(-(cos(q2)*cos(q3)*sin(q2)) + sin(q1)*sin(q3)) + dx4*(cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3)) +
            dx5*(cos(q3 - q4)*(cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3)) + (-(cos(q2)*cos(q3)*sin(q2)) + sin(q1)*sin(q3))*sin(q3 - q4)) + dz5*(cos(q3 - q4)*(-(cos(q2)*cos(q3)*sin(q2)) + sin(q1)*sin(q3)) - (cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3))*sin(q3 - q4))
            ,-(cos(q2)*cos(q3 - q4)*sin(q3)) + cos(q2)*cos(q3)*sin(q3 - q4),sin(q2),cos(q2)*cos(q3)*cos(q3 - q4) + cos(q2)*sin(q3)*sin(q3 - q4),
            dz1 + dz2 + dz3*cos(q2) + dz4*cos(q2)*cos(q3) + dy3*sin(q2) + dy4*sin(q2) + dy5*sin(q2) - dx4*cos(q2)*sin(q3) + dx5*(-(cos(q2)*cos(q3 - q4)*sin(q3)) + cos(q2)*cos(q3)*sin(q3 - q4)) + dz5*(cos(q2)*cos(q3)*cos(q3 - q4) + cos(q2)*sin(q3)*sin(q3 - q4)),
            0,0,0,1;*/
    temp << cos(q4)*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) - (cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3))*sin(q4),-(cos(q2)*sin(q1)),
            cos(q4)*(cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3)) + (cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3))*sin(q4),
            dx1 + dx2*cos(q1) + dx3*cos(q1) - dy2*sin(q1) - dy3*cos(q2)*sin(q1) - dy4*cos(q2)*sin(q1) - dy5*cos(q2)*sin(q1) + dz3*sin(q1)*sin(q2) + dz4*(cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3)) +
            dx4*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) + dx5*(cos(q4)*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) - (cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3))*sin(q4)) +
            dz5*(cos(q4)*(cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3)) + (cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3))*sin(q4)),
            cos(q4)*(cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3)) - (-(cos(q2)*cos(q3)*sin(q2)) + sin(q1)*sin(q3))*sin(q4),pow(cos(q2),2),
            cos(q4)*(-(cos(q2)*cos(q3)*sin(q2)) + sin(q1)*sin(q3)) + (cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3))*sin(q4),
            dy1 + dy2*cos(q2) + dy3*pow(cos(q2),2) + dy4*pow(cos(q2),2) + dy5*pow(cos(q2),2) + dx2*sin(q1) + dx3*sin(q1) - dz3*cos(q2)*sin(q2) + dz4*(-(cos(q2)*cos(q3)*sin(q2)) + sin(q1)*sin(q3)) +
            dx4*(cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3)) + dx5*(cos(q4)*(cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3)) - (-(cos(q2)*cos(q3)*sin(q2)) + sin(q1)*sin(q3))*sin(q4)) +
            dz5*(cos(q4)*(-(cos(q2)*cos(q3)*sin(q2)) + sin(q1)*sin(q3)) + (cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3))*sin(q4)),
            -(cos(q2)*cos(q4)*sin(q3)) - cos(q2)*cos(q3)*sin(q4),sin(q2),cos(q2)*cos(q3)*cos(q4) - cos(q2)*sin(q3)*sin(q4),
            dz1 + dz2 + dz3*cos(q2) + dz4*cos(q2)*cos(q3) + dy3*sin(q2) + dy4*sin(q2) + dy5*sin(q2) - dx4*cos(q2)*sin(q3) + dx5*(-(cos(q2)*cos(q4)*sin(q3)) - cos(q2)*cos(q3)*sin(q4)) +
            dz5*(cos(q2)*cos(q3)*cos(q4) - cos(q2)*sin(q3)*sin(q4)),0,0,0,1;
    pFoot << temp(0,3), temp(1,3), temp(2,3), qMotor(0);
}
void jacobianRight(const Eigen::Matrix<double, 4, 1>& qMotor, Eigen::Matrix<double, 4, 4>& Jacobian){
    double roll, pitch, yaw, q1, q2, q3,q4;
    double dx1, dy1, dz1, dx2, dy2, dz2 , dx3, dy3, dz3, dx4, dy4, dz4, dx5, dy5, dz5;


//    roll = 0;
//    pitch = 0;
//    yaw = 0;

    q1 = qMotor(0);
    q2 = qMotor(1);
    q3 = qMotor(2);
    q4 = qMotor(3);

    dx1 = dx1_R;
    dy1 = dy1_R;
    dz1 = dz1_R;

    dx2 = dx2_R;
    dy2 = dy2_R;
    dz2 = dz2_R;

    dx3 = dx3_R;
    dy3 = dy3_R;
    dz3 = dz3_R;

    dx4 = dx4_R;
    dy4 = dy4_R;
    dz4 = dz4_R;

    dx5 = dx5_R;
    dy5 = dy5_R;
    dz5 = dz5_R;

    Jacobian <<   -(dy2*cos(q1)) - dy3*cos(q1)*cos(q2) - dy4*cos(q1)*cos(q2) - dy5*cos(q1)*cos(q2) - dx2*sin(q1) - dx3*sin(q1) + dz3*cos(q1)*sin(q2) + dz4*(cos(q1)*cos(q3)*sin(q2) - sin(q1)*sin(q3)) +
                  dx4*(-(cos(q3)*sin(q1)) - cos(q1)*sin(q2)*sin(q3)) + dx5*(cos(q4)*(-(cos(q3)*sin(q1)) - cos(q1)*sin(q2)*sin(q3)) - (cos(q1)*cos(q3)*sin(q2) - sin(q1)*sin(q3))*sin(q4)) +
                  dz5*(cos(q4)*(cos(q1)*cos(q3)*sin(q2) - sin(q1)*sin(q3)) + (-(cos(q3)*sin(q1)) - cos(q1)*sin(q2)*sin(q3))*sin(q4)),
            dz3*cos(q2)*sin(q1) + dz4*cos(q2)*cos(q3)*sin(q1) + dy3*sin(q1)*sin(q2) + dy4*sin(q1)*sin(q2) + dy5*sin(q1)*sin(q2) - dx4*cos(q2)*sin(q1)*sin(q3) +
            dx5*(-(cos(q2)*cos(q4)*sin(q1)*sin(q3)) - cos(q2)*cos(q3)*sin(q1)*sin(q4)) + dz5*(cos(q2)*cos(q3)*cos(q4)*sin(q1) - cos(q2)*sin(q1)*sin(q3)*sin(q4)),
            dx4*(-(cos(q3)*sin(q1)*sin(q2)) - cos(q1)*sin(q3)) + dz4*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) +
            dz5*(cos(q4)*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) + (-(cos(q3)*sin(q1)*sin(q2)) - cos(q1)*sin(q3))*sin(q4)) +
            dx5*(cos(q4)*(-(cos(q3)*sin(q1)*sin(q2)) - cos(q1)*sin(q3)) - (cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3))*sin(q4)),
            dz5*(cos(q4)*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) - (cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3))*sin(q4)) +
            dx5*(-(cos(q4)*(cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3))) - (cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3))*sin(q4)),
            dx2*cos(q1) + dx3*cos(q1) + dx4*cos(q1)*cos(q3) + dz4*cos(q1)*sin(q3) + dz5*(cos(q1)*cos(q4)*sin(q3) + cos(q1)*cos(q3)*sin(q4)) + dx5*(cos(q1)*cos(q3)*cos(q4) - cos(q1)*sin(q3)*sin(q4)),
            -(dz3*pow(cos(q2),2)) - dy2*sin(q2) - 2*dy3*cos(q2)*sin(q2) - 2*dy4*cos(q2)*sin(q2) - 2*dy5*cos(q2)*sin(q2) + dz3*pow(sin(q2),2) +
            dz4*(-(pow(cos(q2),2)*cos(q3)) + cos(q3)*pow(sin(q2),2)) + dx4*(pow(cos(q2),2)*sin(q3) - pow(sin(q2),2)*sin(q3)) +
            dx5*(cos(q4)*(pow(cos(q2),2)*sin(q3) - pow(sin(q2),2)*sin(q3)) - (-(pow(cos(q2),2)*cos(q3)) + cos(q3)*pow(sin(q2),2))*sin(q4)) +
            dz5*(cos(q4)*(-(pow(cos(q2),2)*cos(q3)) + cos(q3)*pow(sin(q2),2)) + (pow(cos(q2),2)*sin(q3) - pow(sin(q2),2)*sin(q3))*sin(q4)),
            dx4*(cos(q2)*cos(q3)*sin(q2) - sin(q1)*sin(q3)) + dz4*(cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3)) +
            dz5*(cos(q4)*(cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3)) + (cos(q2)*cos(q3)*sin(q2) - sin(q1)*sin(q3))*sin(q4)) +
            dx5*(cos(q4)*(cos(q2)*cos(q3)*sin(q2) - sin(q1)*sin(q3)) - (cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3))*sin(q4)),
            dz5*(cos(q4)*(cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3)) - (-(cos(q2)*cos(q3)*sin(q2)) + sin(q1)*sin(q3))*sin(q4)) +
            dx5*(-(cos(q4)*(-(cos(q2)*cos(q3)*sin(q2)) + sin(q1)*sin(q3))) - (cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3))*sin(q4)),
            0,dy3*cos(q2) + dy4*cos(q2) + dy5*cos(q2) - dz3*sin(q2) - dz4*cos(q3)*sin(q2) + dx4*sin(q2)*sin(q3) + dx5*(cos(q4)*sin(q2)*sin(q3) + cos(q3)*sin(q2)*sin(q4)) +
              dz5*(-(cos(q3)*cos(q4)*sin(q2)) + sin(q2)*sin(q3)*sin(q4)),-(dx4*cos(q2)*cos(q3)) - dz4*cos(q2)*sin(q3) + dz5*(-(cos(q2)*cos(q4)*sin(q3)) - cos(q2)*cos(q3)*sin(q4)) +
                                                                         dx5*(-(cos(q2)*cos(q3)*cos(q4)) + cos(q2)*sin(q3)*sin(q4)),dz5*(-(cos(q2)*cos(q4)*sin(q3)) - cos(q2)*cos(q3)*sin(q4)) + dx5*(-(cos(q2)*cos(q3)*cos(q4)) + cos(q2)*sin(q3)*sin(q4)),
            0,1,0,0;

}
void Robbie::controlStandingLeg() {
//    cout << "height: " << posBase2ST_RobBIE_(2) << endl;
    double dz = cmdHeight_ - posBase2ST_world_(2);

    Eigen::Matrix<double, 4, 1> forceFoot;
    forceFoot(0) = -(kp_Torso_Pitch * qFloatBase_(1) + kd_Torso_Pitch * dqFloatBase_(1)) / posBase2ST_RobBIE_(2);
    forceFoot(1) = (kp_Torso_Roll * qFloatBase_(0) + kd_Torso_Roll * dqFloatBase_(0)) / posBase2ST_RobBIE_(2);
    forceFoot(2) = -m_RobBIE * gravity - kp_Height * (dz) - kd_Height * (0 - velFloatBaseFil_(2));
    if (standingLeg_ == RIGHTLEG)
        forceFoot(3) = kp_Yaw_Standing * (0 - qMotorRight_(1)) + kd_Yaw_Standing * (0 - dqMotorRight_(1));
    else
        forceFoot(3) = kp_Yaw_Standing * (0 - qMotorLeft_(1)) + kd_Yaw_Standing * (0 - dqMotorLeft_(1));

    if (standingLeg_ == LEFTLEG) {
        desTorqueLeftMotor_ = JacobiJ2F_left_.transpose() * forceFoot;
//        cout << "forceFoot: " << forceFoot.transpose() << endl;
//        cout << "desTorque: " << desTorqueLeftMotor_.transpose() << endl;
        cout << "forceFoot: " << forceFoot.transpose() << endl;
        cout << "desTorque: " << desTorqueLeftMotor_.transpose() << endl;
        cout << "qMotor:" << qMotorRight_.transpose() << endl;
        cout << "jacobi: \n" << JacobiJ2F_right_ << endl;
    } else {
        desTorqueRightMotor_ = JacobiJ2F_right_.transpose() * forceFoot;
        cout << "forceFoot: " << forceFoot.transpose() << endl;
        cout << "desTorque: " << desTorqueRightMotor_.transpose() << endl;
        cout << "qMotor:" << qMotorRight_.transpose() << endl;
        cout << "jacobi: \n" << JacobiJ2F_right_ << endl;
    }
}

   /* void calculateGroundForce()
    {forceGroundLeft_ = JacobiJ2F_left_.transpose().inverse() * (-torqueMotorLeft_);
        forceGroundRight_ = JacobiJ2F_right_.transpose().inverse() * (-torqueMotorRight_);

        forceGroundLeftFil_ = forceGroundLeftFil_*(1-filterGRF) + forceGroundLeft_*filterGRF;
        forceGroundRightFil_ = forceGroundRightFil_*(1-filterGRF) + forceGroundRight_*filterGRF;
        cout << "forceGroundLeftFil_  :"<< forceGroundLeftFil_.transpose()  << endl ;
        cout  << " forceGroundRightFil_ : "<< forceGroundRightFil_.transpose()  << endl;}*/



    int main(){
    Eigen::Matrix<double, 4, 1> qMotor;
    Eigen::Matrix<double,4,1> pfoot;
    qMotor << 0,0,-0.5,1.1;
    forward_kinematics_right(qMotor,pfoot);
    std::cout << pfoot << std::endl;

       Eigen::Matrix<double, 4, 4> &Jacobian = a;

   // calculateGroundForce()




}