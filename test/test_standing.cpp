//
// Created by yeyinong on 2024/3/8.
//
#include "parameters/parameters.h"
#include "math.h"
#include "eigen3/Eigen/Core"
#include "iostream"
#include "eigen3/Eigen/LU"
void forward_kinematics_left(const Eigen::Matrix<double, 4, 1> qMotor, Eigen::Matrix<double, 4, 1> pFoot){
    //该方法通过读取基座姿态和左腿电机的角度，结合各个关节的偏移量和三角函数关系，计算左脚在机器人坐标系中的位置。这个过程涉及正向运动学的计算，通过变换矩阵实现。
    double roll, pitch, yaw, q1, q2, q3,q4;
    double dx1, dy1, dz1, dx2, dy2, dz2 , dx3, dy3, dz3, dx4, dy4, dz4, dx5, dy5, dz5;





    q1 = qMotor(0);
    q2 = qMotor(1);
    q3 = qMotor(2);
    q4 = qMotor(3);

//    cout << roll << " " << pitch << " " << yaw << endl;
//    cout << q1 << " " << q2 << " " << q3 << " " << q4 <<endl;

    dx1 = dx1_L;
    dy1 = dy1_L;
    dz1 = dz1_L;

    dx2 = dx2_L;
    dy2 = dy2_L;
    dz2 = dz2_L;

    dx3 = dx3_L;
    dy3 = dy3_L;
    dz3 = dz3_L;

    dx4 = dx4_L;
    dy4 = dy4_L;
    dz4 = dz4_L;

    dx5 = dx5_L;
    dy5 = dy5_L;
    dz5 = dz5_L;
    Eigen::Matrix<double, 4, 4> temp;
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
    pFoot << temp(0,3), temp(1,3), temp(2,3), qMotor(1);

}
int main(){
   /* double roll, pitch, yaw, q1, q2, q3,q4;
    double dx1, dy1, dz1, dx2, dy2, dz2 , dx3, dy3, dz3, dx4, dy4, dz4, dx5, dy5, dz5;
    Eigen::Matrix<double, 4, 1> desTorque;
    Eigen::Matrix<double, 4, 4> Jacobi;

    q1 = 1;
    q2 = 1;
    q3 = -25./180*M_PI;
    q4 = 30./180.*M_PI - q3;

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

    Jacobi <<
           -(dy2*cos(q1)) - dy3*cos(q1)*cos(q2) - dy4*cos(q1)*cos(q2) - dy5*cos(q1)*cos(q2) - dx2*sin(q1) - dx3*sin(q1) + dz3*cos(q1)*sin(q2) + dz4*(cos(q1)*cos(q3)*sin(q2) - sin(q1)*sin(q3)) + dx4*(-(cos(q3)*sin(q1)) - cos(q1)*sin(q2)*sin(q3)) +
           dx5*(cos(q4)*(-(cos(q3)*sin(q1)) - cos(q1)*sin(q2)*sin(q3)) - (cos(q1)*cos(q3)*sin(q2) - sin(q1)*sin(q3))*sin(q4)) + dz5*(cos(q4)*(cos(q1)*cos(q3)*sin(q2) - sin(q1)*sin(q3)) + (-(cos(q3)*sin(q1)) - cos(q1)*sin(q2)*sin(q3))*sin(q4)),
            dz3*cos(q2)*sin(q1) + dz4*cos(q2)*cos(q3)*sin(q1) + dy3*sin(q1)*sin(q2) + dy4*sin(q1)*sin(q2) + dy5*sin(q1)*sin(q2) - dx4*cos(q2)*sin(q1)*sin(q3) + dx5*(-(cos(q2)*cos(q4)*sin(q1)*sin(q3)) - cos(q2)*cos(q3)*sin(q1)*sin(q4)) +
            dz5*(cos(q2)*cos(q3)*cos(q4)*sin(q1) - cos(q2)*sin(q1)*sin(q3)*sin(q4)),dx4*(-(cos(q3)*sin(q1)*sin(q2)) - cos(q1)*sin(q3)) + dz4*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) +
                                                                                    dz5*(cos(q4)*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) + (-(cos(q3)*sin(q1)*sin(q2)) - cos(q1)*sin(q3))*sin(q4)) + dx5*(cos(q4)*(-(cos(q3)*sin(q1)*sin(q2)) - cos(q1)*sin(q3)) - (cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3))*sin(q4)),
            dz5*(cos(q4)*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) - (cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3))*sin(q4)) + dx5*(-(cos(q4)*(cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3))) - (cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3))*sin(q4)),

            dx2*cos(q1) + dx3*cos(q1) + dx4*cos(q1)*cos(q3) + dz4*cos(q1)*sin(q3) + dz5*(cos(q1)*cos(q4)*sin(q3) + cos(q1)*cos(q3)*sin(q4)) + dx5*(cos(q1)*cos(q3)*cos(q4) - cos(q1)*sin(q3)*sin(q4)),
            -(dz3*pow(cos(q2),2)) - dy2*sin(q2) - 2*dy3*cos(q2)*sin(q2) - 2*dy4*cos(q2)*sin(q2) - 2*dy5*cos(q2)*sin(q2) + dz3*pow(sin(q2),2) + dz4*(-(pow(cos(q2),2)*cos(q3)) + cos(q3)*pow(sin(q2),2)) +
            dx4*(pow(cos(q2),2)*sin(q3) - pow(sin(q2),2)*sin(q3)) + dx5*(cos(q4)*(pow(cos(q2),2)*sin(q3) - pow(sin(q2),2)*sin(q3)) - (-(pow(cos(q2),2)*cos(q3)) + cos(q3)*pow(sin(q2),2))*sin(q4)) +
            dz5*(cos(q4)*(-(pow(cos(q2),2)*cos(q3)) + cos(q3)*pow(sin(q2),2)) + (pow(cos(q2),2)*sin(q3) - pow(sin(q2),2)*sin(q3))*sin(q4)),
            dx4*(cos(q2)*cos(q3)*sin(q2) - sin(q1)*sin(q3)) + dz4*(cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3)) + dz5*(cos(q4)*(cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3)) + (cos(q2)*cos(q3)*sin(q2) - sin(q1)*sin(q3))*sin(q4)) +
            dx5*(cos(q4)*(cos(q2)*cos(q3)*sin(q2) - sin(q1)*sin(q3)) - (cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3))*sin(q4)),
            dz5*(cos(q4)*(cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3)) - (-(cos(q2)*cos(q3)*sin(q2)) + sin(q1)*sin(q3))*sin(q4)) + dx5*(-(cos(q4)*(-(cos(q2)*cos(q3)*sin(q2)) + sin(q1)*sin(q3))) - (cos(q3)*sin(q1) + cos(q2)*sin(q2)*sin(q3))*sin(q4)),

            0,dy3*cos(q2) + dy4*cos(q2) + dy5*cos(q2) - dz3*sin(q2) - dz4*cos(q3)*sin(q2) + dx4*sin(q2)*sin(q3) + dx5*(cos(q4)*sin(q2)*sin(q3) + cos(q3)*sin(q2)*sin(q4)) + dz5*(-(cos(q3)*cos(q4)*sin(q2)) + sin(q2)*sin(q3)*sin(q4)),
            -(dx4*cos(q2)*cos(q3)) - dz4*cos(q2)*sin(q3) + dz5*(-(cos(q2)*cos(q4)*sin(q3)) - cos(q2)*cos(q3)*sin(q4)) + dx5*(-(cos(q2)*cos(q3)*cos(q4)) + cos(q2)*sin(q3)*sin(q4)),
            dz5*(-(cos(q2)*cos(q4)*sin(q3)) - cos(q2)*cos(q3)*sin(q4)) + dx5*(-(cos(q2)*cos(q3)*cos(q4)) + cos(q2)*sin(q3)*sin(q4)),

            1,0,0,0;

    desTorque << 0,0,-1.3,-7;
    //std::cout << Jacobi.transpose().inverse()*desTorque;*/
    Eigen::Matrix<double, 4, 1> qMotor;
    Eigen::Matrix<double, 4, 1> p;
    qMotor << 0,0,0,0;
    forward_kinematics_left(qMotor,p);
    std::cout << p;



}