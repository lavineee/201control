//
// Created by yeyinong on 2023/3/29.
//

#include "robbie/Robbie.h"

vector<double> sf_coefficient{0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
vector<double> sl_coefficient{0,0,1,1,1,1};

void inverseCurve_2pts(double x1, double y1, double x2, double y2, double *a, double *b){
    *a = -y1*y2*(x1-x2)/(y1-y2);
    *b = -(y1*x1 - y2*x2)/(y1-y2);
}

double SlideToeCovariance(double t, char direction){
    switch (direction) {
        double cov_1;
        double cov_2;
        case 'x':
            cov_1 = pow(c_vx_slide_toe_1, 2);
            cov_2 = pow(c_vx_slide_toe_2, 2);
            if(t < c_Tx_slide_toe_2){
                double a, b;
                inverseCurve_2pts(0, cov_1, c_Tx_slide_toe_2, cov_2, &a, &b);
                return a/(t+b);
            }else{
                return cov_2;
            }
        case 'y':
            cov_1 = pow(c_vy_slide_toe_1, 2);
            cov_2 = pow(c_vy_slide_toe_2, 2);
            if(t < c_Tx_slide_toe_2){
                double a, b;
                inverseCurve_2pts(0, cov_1, c_Ty_slide_toe_2, cov_2, &a, &b);
                return a/(t+b);
            }else{
                return cov_2;
            }
        case 'z':
            cov_1 = pow(c_vz_slide_toe_1, 2);
            cov_2 = pow(c_vz_slide_toe_2, 2);
            if(t < c_Tz_slide_toe_2){
                double a, b;
                inverseCurve_2pts(0, cov_1, c_Tz_slide_toe_2, cov_2, &a, &b);
                return a/(t+b);
            }else{
                return cov_2;
            }
    }
}

double firstOrderFilter(double prev, double now, double para){
    return prev*(1-para) + now*para;
}

double limitRange(double value, double lb, double ub){
    return value<lb?lb:(value>ub?ub:value);
}

double calculateC(double n, double m){
    double a = 1;
    double b = 1;
    double temp = n-m;
    for (; n > temp; n--) {
        a *= n;
    }
    for (;  m>0 ; m--) {
        b *= m;
    }
    return a/b;
}

double bezierPoint(vector<double> controlPoints,  double tau){

    int size = controlPoints.size();
    int m = size -1;
    double res = 0;
    for (int i = 0; i <= m; ++i) {
        if (i == 0)
            res += calculateC(m, i) * pow(1-tau, m-i) * controlPoints[i];
        else if (i == m)
            res += calculateC(m, i) * pow(tau, i)*controlPoints[i];
        else
            res += calculateC(m, i)*pow(1-tau, m-i)* pow(tau, i)*controlPoints[i];
    }
    return res;
}

Eigen::Matrix<double, 4, 1> bezierVal(Eigen::Matrix<double, 4, 6> alpha, double s){
    //1 5 10 10 5 1
    vector<int> mul = vector<int>{1, 5, 10, 10, 5, 1};
    Eigen::Matrix<double, 4, 1> result;
    result << 0,0,0,0;
    for (int i = 0; i < 6; ++i) {
        result(0) += alpha(0,i)*mul[i]*pow(s, i)* pow(1-s, 5-i);
        result(1) += alpha(1,i)*mul[i]*pow(s, i)* pow(1-s, 5-i);
        result(2) += alpha(2,i)*mul[i]*pow(s, i)* pow(1-s, 5-i);
        result(3) += alpha(3,i)*mul[i]*pow(s, i)* pow(1-s, 5-i);
    }
    return result;
}

Eigen::Matrix<double, 4, 1> bezierVal(Eigen::Matrix<double, 4, 5> alpha, double s){
    //1 5 10 10 5 1
    vector<int> mul = vector<int>{1, 4, 6, 4, 1};
    Eigen::Matrix<double, 4, 1> result;
    result << 0,0,0,0;
    for (int i = 0; i < 5; ++i) {
        result(0) += alpha(0,i)*mul[i]*pow(s, i)* pow(1-s, 5-i);
        result(1) += alpha(1,i)*mul[i]*pow(s, i)* pow(1-s, 5-i);
        result(2) += alpha(2,i)*mul[i]*pow(s, i)* pow(1-s, 5-i);
        result(3) += alpha(3,i)*mul[i]*pow(s, i)* pow(1-s, 5-i);
    }
    return result;
}

Eigen::Matrix<double, 4, 5> diff(Eigen::Matrix<double, 4, 6> alpha){
    Eigen::Matrix<double, 4, 5> newAlpha;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 5; ++j) {
            newAlpha(i, j) = alpha(i, j+1) - alpha(i, j);
        }
    }
    return newAlpha;
}

void Robbie::forward_kinematics_left(const Eigen::Matrix<double, 4, 1>& qMotor, Eigen::Matrix<double, 4, 1>& pFoot){
    //该方法通过读取基座姿态和左腿电机的角度，结合各个关节的偏移量和三角函数关系，计算左脚在机器人坐标系中的位置。这个过程涉及正向运动学的计算，通过变换矩阵实现。
    double roll, pitch, yaw, q1, q2, q3,q4;
    double dx1, dy1, dz1, dx2, dy2, dz2 , dx3, dy3, dz3, dx4, dy4, dz4, dx5, dy5, dz5;

    dx1 = dx1_L;

    roll = qFloatBase_(0);
    pitch = qFloatBase_(1);
    yaw = qFloatBase_(2);

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

void Robbie::forward_kinematics_right(const Eigen::Matrix<double, 4, 1>& qMotor, Eigen::Matrix<double, 4, 1>& pFoot){
    double roll, pitch, yaw, q1, q2, q3,q4;
    double dx1, dy1, dz1, dx2, dy2, dz2 , dx3, dy3, dz3, dx4, dy4, dz4, dx5, dy5, dz5;

    roll = qFloatBase_(0);
    pitch = qFloatBase_(1);
    yaw = qFloatBase_(2);

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
    pFoot << temp(0,3), temp(1,3), temp(2,3), qMotor(1);
}

void Robbie::jacobianLeft(const Eigen::Matrix<double, 4, 1>& qMotor, Eigen::Matrix<double, 4, 4>& Jacobian){
    double roll, pitch, yaw, q1, q2, q3,q4;
    double dx1, dy1, dz1, dx2, dy2, dz2 , dx3, dy3, dz3, dx4, dy4, dz4, dx5, dy5, dz5;

    dx1 = dx1_L;

    roll = qFloatBase_(0);
    pitch = qFloatBase_(1);
    yaw = qFloatBase_(2);

//    roll = 0;
//    pitch = 0;
//    yaw = 0;

    q1 = qMotor(0);
    q2 = qMotor(1);
    q3 = qMotor(2);
    q4 = qMotor(3);

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
    Jacobian <<    -(dy2*cos(q1)) - dy3*cos(q1)*cos(q2) - dy4*cos(q1)*cos(q2) - dy5*cos(q1)*cos(q2) - dx2*sin(q1) - dx3*sin(q1) + dz3*cos(q1)*sin(q2) + dz4*(cos(q1)*cos(q3)*sin(q2) - sin(q1)*sin(q3)) +
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

void Robbie::jacobianRight(const Eigen::Matrix<double, 4, 1>& qMotor, Eigen::Matrix<double, 4, 4>& Jacobian){
    double roll, pitch, yaw, q1, q2, q3,q4;
    double dx1, dy1, dz1, dx2, dy2, dz2 , dx3, dy3, dz3, dx4, dy4, dz4, dx5, dy5, dz5;

    roll = qFloatBase_(0);
    pitch = qFloatBase_(1);
    yaw = qFloatBase_(2);

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

Robbie::Robbie(){
    velFloatBaseFil_ << 0,0,0;

    standingLeg_ = true;
    vector<double> doubleBaseMessage(10, 0);
    vector<double> motorMessage(24, 0);
    updateMotorMessage(motorMessage);
    updateFloatBaseMessage(doubleBaseMessage);
    calculateGroundForce();
    updateStandingFoot();
    calculatePosFloatBase();
    calculateVelFloatBase();
    calculatePosSwingFoot();
}

Robbie::Robbie(bool standingLeg, vector<double> &motorMessage, vector<double> &doubleBaseMessage) {
    velFloatBaseFil_ << 0,0,0;
    standingLeg_ = standingLeg;

    updateMotorMessage(motorMessage);
    updateFloatBaseMessage(doubleBaseMessage);
    updateJacLeft();
    updateJacRight();
    calculateGroundForce();
    updateStandingFoot();

    calculatePosLeftFoot2Base_RobBIE();
    calculatePosRightFoot2Base_RobBIE();
    calculatePosFloatBase();
    calculateVelFloatBase();

}

void Robbie::update(vector<double> &motorMessage, vector<double> &doubleBaseMessage, double deltaTime) {
    tau_ += deltaTime / stepTime;
    t_ += deltaTime;

    updateMotorMessage(motorMessage);
    updateFloatBaseMessage(doubleBaseMessage);
    updateJacLeft();
    updateJacRight();
    calculateGroundForce();
    updateStandingFoot();

    calculatePosLeftFoot2Base_RobBIE();
    calculatePosRightFoot2Base_RobBIE();
    calculatePosFloatBase();
    calculateVelFloatBase();


//    estimator(deltaTime);
    calculatePosSwingFoot();
//    updateHZD();
    controlStandingLeg();
    controlSwingLeg();
    if(isJustChangeLeg)
        isJustChangeLeg = false;

}

void Robbie::updateMotorMessage(vector<double> &motorMessage) {
    int motorNumPerLeg = 4;
    for (int i = 0; i < motorNumPerLeg; ++i) {
        qMotorLeft_(i) = motorMessage[i];
        qMotorRight_(i) = motorMessage[i + motorNumPerLeg*1];
        dqMotorLeft_(i) = motorMessage[i + motorNumPerLeg*2];
        dqMotorRight_(i) = motorMessage[i + motorNumPerLeg*3];
        torqueMotorLeft_(i) = motorMessage[i + motorNumPerLeg*4];
        torqueMotorRight_(i) = motorMessage[i + motorNumPerLeg*5];
    }

    cout << "qMotorLeft: " << qMotorLeft_.transpose() << endl;
    cout << "qMotorRight: " << qMotorRight_.transpose() << endl;
}

void Robbie::updateFloatBaseMessage(vector<double> &floatBaseMessage) {

    int sensorIndex_framequat= 0;
    double base_state_quat_w = floatBaseMessage[0];
    double base_state_quat_x = floatBaseMessage[1];
    double base_state_quat_y = floatBaseMessage[2];
    double base_state_quat_z = floatBaseMessage[3];
    Eigen::Quaternion quat = Eigen::Quaternion(base_state_quat_w,base_state_quat_x,base_state_quat_y,base_state_quat_z);

    qFloatBase_(0) = atan2(2*(base_state_quat_w*base_state_quat_x+base_state_quat_y*base_state_quat_z), 1-2*(base_state_quat_x*base_state_quat_x+base_state_quat_y*base_state_quat_y));
    qFloatBase_(1) = asin(2 * (base_state_quat_w * base_state_quat_y - base_state_quat_z * base_state_quat_x));
    qFloatBase_(2) = atan2(2*(base_state_quat_w*base_state_quat_z+base_state_quat_x*base_state_quat_y), 1.0-2.0*(base_state_quat_y*base_state_quat_y+base_state_quat_z*base_state_quat_z));

    for (int i = 0; i < 3; ++i) {
        dqFloatBase_(i) = floatBaseMessage[i+4];
        a_robBIE_(i) = floatBaseMessage[i+7];
    }
    Eigen::Matrix<double, 3, 3> temp;
    temp = Eigen::AngleAxisd(qFloatBase_[2], Eigen::Vector3d::UnitZ())*Eigen::AngleAxisd(qFloatBase_[1], Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(qFloatBase_[0], Eigen::Vector3d::UnitX());
    temp = Eigen::AngleAxisd(0, Eigen::Vector3d::UnitZ())*Eigen::AngleAxisd(qFloatBase_[1], Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(qFloatBase_[0], Eigen::Vector3d::UnitX());
//    rotationMatrix_ = temp;
    rotationMatrix_ = Eigen::Matrix3d::Identity();
    a_world_ = rotationMatrix_*a_robBIE_;
}

void Robbie::updateJacLeft() {
    double roll, pitch, yaw, q1, q2, q3,q4;
    double dx1, dy1, dz1, dx2, dy2, dz2 , dx3, dy3, dz3, dx4, dy4, dz4, dx5, dy5, dz5;

    dx1 = dx1_L;

    roll = qFloatBase_(0);
    pitch = qFloatBase_(1);
    yaw = qFloatBase_(2);

    q1 = qMotorLeft_(0);
    q2 = qMotorLeft_(1);
    q3 = qMotorLeft_(2);
    q4 = qMotorLeft_(3);

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

    JacobiJ2F_left_ <<-(dy2*cos(q1)) - dy3*cos(q1)*cos(q2) - dy4*cos(q1)*cos(q2) - dy5*cos(q1)*cos(q2) - dx2*sin(q1) - dx3*sin(q1) + dz3*cos(q1)*sin(q2) + dz4*(cos(q1)*cos(q3)*sin(q2) - sin(q1)*sin(q3)) +
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

void Robbie::updateJacRight()  {
    double roll, pitch, yaw, q1, q2, q3,q4;
    double dx1, dy1, dz1, dx2, dy2, dz2 , dx3, dy3, dz3, dx4, dy4, dz4, dx5, dy5, dz5;

    roll = qFloatBase_(0);
    pitch = qFloatBase_(1);
    yaw = qFloatBase_(2);

    q1 = qMotorRight_(0);
    q2 = qMotorRight_(1);
    q3 = qMotorRight_(2);
    q4 = qMotorRight_(3);

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

    JacobiJ2F_right_ << -(dy2*cos(q1)) - dy3*cos(q1)*cos(q2) - dy4*cos(q1)*cos(q2) - dy5*cos(q1)*cos(q2) - dx2*sin(q1) - dx3*sin(q1) + dz3*cos(q1)*sin(q2) + dz4*(cos(q1)*cos(q3)*sin(q2) - sin(q1)*sin(q3)) +
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

void Robbie::calculateGroundForce() {

    forceGroundLeft_ = JacobiJ2F_left_.transpose().inverse() * (-torqueMotorLeft_);
    forceGroundRight_ = JacobiJ2F_right_.transpose().inverse() * (-torqueMotorRight_);

    forceGroundLeftFil_ = forceGroundLeftFil_*(1-filterGRF) + forceGroundLeft_*filterGRF;
    forceGroundRightFil_ = forceGroundRightFil_*(1-filterGRF) + forceGroundRight_*filterGRF;
    cout << "forceGroundLeftFil_  :"<< forceGroundLeftFil_.transpose()  << endl ;
    cout  << " forceGroundRightFil_ : "<< forceGroundRightFil_.transpose()  << endl;
}

void Robbie::updateStandingFoot() {
    static int count = 0;
    if(standingLeg_ == LEFTLEG){
        if(forceGroundRight_[2]>20 && tau_ > 0.5)
            count++;
        else
            count=0;
        if(count>=3 || tau_ > 1.1){
            count = 0;
            isJustChangeLeg = true;
            sumTime_ = t_;
            tau_ = 0;
            countStep_++;
            standingLeg_ = RIGHTLEG;
            if(cmdUpdate_){
                cmdVelocity_ = cmdVelocityBuf_;
                cmdUpdate_ = false;
            }
        }
    }else{
        if(forceGroundLeft_[2]>20 && tau_ > 0.5)
            count++;
        else
            count=0;
        if(count>=3 || tau_ > 1.1){
            count = 0;
            isJustChangeLeg = true;
            sumTime_ = t_;
            tau_ = 0;
            countStep_++;
            standingLeg_ = LEFTLEG;
            if(cmdUpdate_){
                cmdVelocity_ = cmdVelocityBuf_;
                cmdUpdate_ = false;
            }
        }
    }
}

void Robbie::calculatePosLeftFoot2Base_RobBIE(){
    double roll, pitch, yaw, q1, q2, q3,q4;
    double dx1, dy1, dz1, dx2, dy2, dz2 , dx3, dy3, dz3, dx4, dy4, dz4, dx5, dy5, dz5;

    dx1 = dx1_L;

    roll = qFloatBase_(0);
    pitch = qFloatBase_(1);
    yaw = qFloatBase_(2);

    q1 = qMotorLeft_(0);
    q2 = qMotorLeft_(1);
    q3 = qMotorLeft_(2);
    q4 = qMotorLeft_(3);
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
    temp <<cos(q4)*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) - (cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3))*sin(q4),-(cos(q2)*sin(q1)),
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
    /*cos(q3 - q4)*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) + (cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3))*sin(q3 - q4),-(cos(q2)*sin(q1)),
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
    posLeftFoot2Base_RobBIE_ << temp(0,3), temp(1,3), temp(2,3);
    posLeftFoot2Base_world_ = rotationMatrix_*posLeftFoot2Base_RobBIE_;
//    cout << "left foot: " << posLeftFoot2Base_world_.transpose() << endl;
}

void Robbie::calculatePosRightFoot2Base_RobBIE(){
    double roll, pitch, yaw, q1, q2, q3,q4;
    double dx1, dy1, dz1, dx2, dy2, dz2 , dx3, dy3, dz3, dx4, dy4, dz4, dx5, dy5, dz5;

    roll = qFloatBase_(0);
    pitch = qFloatBase_(1);
    yaw = qFloatBase_(2);

    q1 = qMotorRight_(0);
    q2 = qMotorRight_(1);
    q3 = qMotorRight_(2);
    q4 = qMotorRight_(3);

//    cout << "qMotorRight: " << qMotorRight_.transpose() << endl;

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
    temp <<cos(q4)*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) - (cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3))*sin(q4),-(cos(q2)*sin(q1)),
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
    /*cos(q3 - q4)*(cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3)) + (cos(q3)*sin(q1)*sin(q2) + cos(q1)*sin(q3))*sin(q3 - q4),-(cos(q2)*sin(q1)),
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
    posRightFoot2Base_RobBIE_ << temp(0,3), temp(1,3), temp(2,3);
    posRightFoot2Base_world_ = rotationMatrix_*posRightFoot2Base_RobBIE_;
//    cout << "right foot: " << posRightFoot2Base_world_.transpose() << endl;
}

void Robbie::calculatePosFloatBase() {
    if(standingLeg_ == LEFTLEG){
        posBase2ST_RobBIE_ = -posLeftFoot2Base_RobBIE_;
        posBase2ST_world_ = -posLeftFoot2Base_world_;
    }else{
        posBase2ST_RobBIE_ = -posRightFoot2Base_RobBIE_;
        posBase2ST_world_ = -posRightFoot2Base_world_;
    }
    return;
}

void Robbie::calculateVelFloatBase() {
        Eigen::Matrix<double, 4, 1> dqMotor;
        dqMotor = standingLeg_ == LEFTLEG ? dqMotorLeft_ : dqMotorRight_;

        Eigen::Matrix<double, 4, 4> JacobiJ2B;
        if (standingLeg_ == LEFTLEG)
            JacobiJ2B = -JacobiJ2F_left_;
        else
            JacobiJ2B = -JacobiJ2F_right_;

        Eigen::Matrix<double, 3, 3> omegaMatrix;
        omegaMatrix <<
                    0, -dqFloatBase_(2), dqFloatBase_(1),
                dqFloatBase_(2), 0, -dqFloatBase_(0),
                -dqFloatBase_(1), dqFloatBase_(0), 0;

        Eigen::Matrix<double, 4, 1> temp1;
        temp1 = JacobiJ2B * dqMotor;
        Eigen::Matrix<double, 3, 1> temp2;
        temp2 << temp1(0), temp1(1), temp1(2);

        velFloatBase_ << omegaMatrix * rotationMatrix_ * posBase2ST_RobBIE_ + rotationMatrix_ * temp2;
        velFloatBaseFil_ = velFloatBaseFil_ * (1 - filterVel) + velFloatBase_ * filterVel;
    return;
}

void Robbie::calculatePosSwingFoot() {
    if(standingLeg_ == LEFTLEG){
        posSW2Base_RobBIE_ = posRightFoot2Base_RobBIE_;
        posSW2Base_world_ = posRightFoot2Base_world_;
    }else{
        posSW2Base_RobBIE_ = posLeftFoot2Base_RobBIE_;
        posSW2Base_world_ = posLeftFoot2Base_world_;
    }
    if(isJustChangeLeg){
        posSW2Base_world_begein_ << posSW2Base_world_;

        if(standingLeg_ == LEFTLEG){
            yawSWFoot_RobBIE_begin_ = qMotorRight_(1);
            yawSTFoot_RobBIE_begin_ = qMotorLeft_(1);
            yawSWFoot_world_begin_ = yawSWFoot_RobBIE_begin_ + qFloatBase_(2);
            yawSTFoot_world_begin_ = yawSTFoot_RobBIE_begin_ + qFloatBase_(2);
        }else{
            yawSWFoot_RobBIE_begin_ = qMotorLeft_(1);
            yawSTFoot_RobBIE_begin_ = qMotorRight_(1);
            yawSWFoot_world_begin_ = yawSWFoot_RobBIE_begin_ + qFloatBase_(2);
            yawSTFoot_world_begin_ = yawSTFoot_RobBIE_begin_ + qFloatBase_(2);
        }
    }
    return;

}

/*void Robbie::estimator(double deltaTime) {
    if(!hasInit){
        Eigen::Matrix<double, 3, 3> Cov_LRT;
        Cov_LRT.setZero();
        Cov_LRT(0,0) = c_Cov_rpx_LRToe_body;
        Cov_LRT(1, 1) = c_Cov_rpy_LRToe_body;
        Cov_LRT(2, 2) = c_Cov_rpz_LRToe_body;

        Cov_LRT_r = rotationMatrix_*Cov_LRT*rotationMatrix_.transpose();

        vOpOpswT_KF_x.setZero();
        vOpOpswT_KF_y.setZero();
        vOpOpswT_KF_z.setZero();
        vOpOpswT_KF_x(1) = posFloatBase_world(0);
        vOpOpswT_KF_y(1) = posFloatBase_world(1);
        vOpOpswT_KF_z(1) = posFloatBase_world(2);

        Cov_LinearAccelerator.setZero();
        Cov_LinearAccelerator(0,0) = c_Cov_LinearAccelerator_xy;
        Cov_LinearAccelerator(1,1) = c_Cov_LinearAccelerator_xy;
        Cov_LinearAccelerator(2,2) = c_Cov_LinearAccelerator_z;
        sigma_vOpOpswT_KF_x << 1000, 0, 0, 0, 1000, 0, 0, 0, 1000;
        sigma_vOpOpswT_KF_y << 1000, 0, 0, 0, 1000, 0, 0, 0, 1000;
        sigma_vOpOpswT_KF_z << 1000, 0, 0, 0, 1000, 0, 0, 0, 1000;
        hasInit = true;
    }
    if(isJustChangeLeg) {
        vOpOpswT_KF_x(2) = vOpOpswT_KF_x(1) - posFloatBase_world(0);
        vOpOpswT_KF_y(2) = vOpOpswT_KF_y(1) - posFloatBase_world(1);
        vOpOpswT_KF_z(2) = vOpOpswT_KF_z(1) - posFloatBase_world(2);

        sigma_vOpOpswT_KF_x(2, 2) = sigma_vOpOpswT_KF_x(1, 1) + Cov_LRT_r(0, 0);
        sigma_vOpOpswT_KF_y(2, 2) = sigma_vOpOpswT_KF_y(1, 1) + Cov_LRT_r(1, 1);
        sigma_vOpOpswT_KF_z(2, 2) = sigma_vOpOpswT_KF_z(1, 1) + Cov_LRT_r(2, 2);
    }

    // Prediction model (x(k) is defined as [v_O; p_O; p_swT])
    // 1) Process model
    //    X(k+1) = V(k)*delta_t
    //    V(k+1) = V(k) + dk*delta_t
    // 2) Measurement model (Take k_timestep p_O2stT_r as z(k))
    //    z(k) = x(k)(2) - x(k)(3)
    Eigen::Matrix<double, 3, 3> At = Eigen::Matrix<double,3 ,3>::Identity();
    At(1, 0) = deltaTime;
    Eigen::Matrix<double, 3, 1> Bt;
    Bt << deltaTime, 0, 0;
    Eigen::Matrix<double, 1, 3> Ct;
    Ct << 0, 1, -1;

    // Input
    double ut_x = a_world_(0);
    double ut_y = a_world_(1);
    double ut_z = a_world_(2)-9.81;



    Eigen::Matrix<double, 3 , 3> Qt_x, Qt_y, Qt_z;
    Qt_x.setZero();
    Qt_y.setZero();
    Qt_z.setZero();

    Qt_x(0, 0) = Cov_LinearAccelerator(0, 0)* pow((1+ abs(a_world_(0))), 2)* pow(deltaTime, 2);
    Qt_x(2, 2) = SlideToeCovariance(tau_*stepTime, 'x') * pow(deltaTime, 2);

    Qt_y(0, 0) = Cov_LinearAccelerator(1, 1)* pow((1+ abs(a_world_(1))), 2)* pow(deltaTime, 2);
    Qt_y(2, 2) = SlideToeCovariance(tau_*stepTime, 'y') * pow(deltaTime, 2);

    Qt_z(0, 0) = Cov_LinearAccelerator(2, 2)* pow((1+ abs(a_world_(2))), 2)* pow(deltaTime, 2);
    Qt_z(2, 2) = SlideToeCovariance(tau_*stepTime, 'z') * pow(deltaTime, 2);

    if (t_ < 0.1){
        Qt_x(2, 2) = 1000;
        Qt_y(2, 2) = 1000;
        Qt_z(2, 2) = 1000;
    }

    // Measurements
    double zt_x = posBase2ST_world_(0);
    double zt_y = posBase2ST_world_(1);
    double zt_z = posBase2ST_world_(2);

    // Rt (Measurement model noise) ????????
    double Rt_x = Cov_LRT_r(0, 0);
    double Rt_y = Cov_LRT_r(1, 1);
    double Rt_z = Cov_LRT_r(2, 2);

    /// KF algorithm (EKF)
    // x
    // Prediction
    Eigen::Matrix<double, 3, 1> vOpOpswT_KF_x_hat = At*vOpOpswT_KF_x + Bt*ut_x;
    Eigen::Matrix<double, 3, 3> sigma_vOpOpswT_KF_x_hat = At*sigma_vOpOpswT_KF_x*At.transpose() + Qt_x;
    // Correction
    vOpOpswT_KF_x = vOpOpswT_KF_x_hat + sigma_vOpOpswT_KF_x_hat*Ct.transpose()*((zt_x - Ct*vOpOpswT_KF_x_hat)/(Ct*sigma_vOpOpswT_KF_x_hat*Ct.transpose() + Rt_x));
    sigma_vOpOpswT_KF_x = (Eigen::Matrix<double, 3, 3>::Identity() - sigma_vOpOpswT_KF_x_hat*Ct.transpose()*(Ct/(Ct*sigma_vOpOpswT_KF_x_hat*Ct.transpose() + Rt_x)))*sigma_vOpOpswT_KF_x_hat;

    // y
    // Prediction
    Eigen::Matrix<double, 3, 1> vOpOpswT_KF_y_hat = At*vOpOpswT_KF_y + Bt*ut_y;
    Eigen::Matrix<double, 3, 3> sigma_vOpOpswT_KF_y_hat = At*sigma_vOpOpswT_KF_y*At.transpose() + Qt_y;
    // Correction
    vOpOpswT_KF_y = vOpOpswT_KF_y_hat + sigma_vOpOpswT_KF_y_hat*Ct.transpose()*((zt_y - Ct*vOpOpswT_KF_y_hat)/(Ct*sigma_vOpOpswT_KF_y_hat*Ct.transpose() + Rt_y));
    sigma_vOpOpswT_KF_y = (Eigen::Matrix<double, 3, 3>::Identity() - sigma_vOpOpswT_KF_y_hat*Ct.transpose()*(Ct/(Ct*sigma_vOpOpswT_KF_y_hat*Ct.transpose() + Rt_y)))*sigma_vOpOpswT_KF_y_hat;

    // x
    // Prediction
    Eigen::Matrix<double, 3, 1> vOpOpswT_KF_z_hat = At*vOpOpswT_KF_z + Bt*ut_z;
    Eigen::Matrix<double, 3, 3> sigma_vOpOpswT_KF_z_hat = At*sigma_vOpOpswT_KF_z*At.transpose() + Qt_z;
    // Correction
    vOpOpswT_KF_z = vOpOpswT_KF_z_hat + sigma_vOpOpswT_KF_z_hat*Ct.transpose()*((zt_z - Ct*vOpOpswT_KF_z_hat)/(Ct*sigma_vOpOpswT_KF_z_hat*Ct.transpose() + Rt_z));
    sigma_vOpOpswT_KF_z = (Eigen::Matrix<double, 3, 3>::Identity() - sigma_vOpOpswT_KF_z_hat*Ct.transpose()*(Ct/(Ct*sigma_vOpOpswT_KF_z_hat*Ct.transpose() + Rt_z)))*sigma_vOpOpswT_KF_z_hat;

    pos_world << vOpOpswT_KF_x(1), vOpOpswT_KF_y(1), vOpOpswT_KF_z(1);
    vel_world << vOpOpswT_KF_x(0), vOpOpswT_KF_y(0), vOpOpswT_KF_z(0);
}*/

void Robbie::controlStandingLeg() {
//    cout << "height: " << posBase2ST_RobBIE_(2) << endl;
    double dz = cmdHeight_ - posBase2ST_world_(2);

    Eigen::Matrix<double, 4, 1> forceFoot;
    forceFoot(0) = -(kp_Torso_Pitch * qFloatBase_(1) + kd_Torso_Pitch * dqFloatBase_(1)) / posBase2ST_RobBIE_(2);
    forceFoot(1) = (kp_Torso_Roll * qFloatBase_(0) + kd_Torso_Roll * dqFloatBase_(0)) / posBase2ST_RobBIE_(2);
    forceFoot(2) = -m_RobBIE*gravity - kp_Height*(dz)-kd_Height*(0-velFloatBaseFil_(2));
    if(standingLeg_ == RIGHTLEG)
        forceFoot(3) = kp_Yaw_Standing * (0-qMotorRight_(1)) + kd_Yaw_Standing*(0-dqMotorRight_(1));
    else
        forceFoot(3) = kp_Yaw_Standing * (0-qMotorLeft_(1)) + kd_Yaw_Standing*(0-dqMotorLeft_(1));

    if (standingLeg_ == LEFTLEG) {
        desTorqueLeftMotor_ = JacobiJ2F_left_.transpose()*forceFoot;
//        cout << "forceFoot: " << forceFoot.transpose() << endl;
//        cout << "desTorque: " << desTorqueLeftMotor_.transpose() << endl;
        cout << "forceFoot: " << forceFoot.transpose() << endl;
        cout << "desTorque: " << desTorqueLeftMotor_.transpose() << endl;
        cout << "qMotor:" << qMotorRight_.transpose() << endl;
        cout << "jacobi: \n" << JacobiJ2F_right_ << endl;
    } else {
        desTorqueRightMotor_ = JacobiJ2F_right_.transpose()*forceFoot;
        cout << "forceFoot: " << forceFoot.transpose() << endl;
        cout << "desTorque: " << desTorqueRightMotor_.transpose() << endl;
        cout << "qMotor:" << qMotorRight_.transpose() << endl;
        cout << "jacobi: \n" << JacobiJ2F_right_ << endl;
    }


}

void Robbie::calculateStateSSP(){
    averageStepTime = countStep_==0 ? stepTime : sumTime_ / countStep_;
    lambda_ = sqrt(gravity / posBase2ST_world_(2));
    sigma_ = lambda_ * sinh(averageStepTime / 2 * lambda_) / cosh(averageStepTime / 2 * lambda_);
    estimatedDesDyStateSSP_ = (cmdWs * sigma_) / 2;

    double t;
    t = averageStepTime- tau_ * stepTime;
    t = t<0?0:t;
    if(standingLeg_ == LEFTLEG){
        y0_ = -posBase2ST_world_(1);
        dy0_ = -velFloatBaseFil_(1);
    }else{
        y0_ = posBase2ST_world_(1);
        dy0_ = velFloatBaseFil_(1);
    }
    double c1 = 1./2*(y0_+ dy0_ / lambda_);
    double c2 = 1./2*(y0_- dy0_ / lambda_);
    double ySSP = c1* pow(M_E, lambda_ * t) + c2 * pow(M_E, -lambda_ * t);
    double dySSP = lambda_ * (c1 * pow(M_E, lambda_ * t) - c2 * pow(M_E, -lambda_ * t));
    estimatedStateSSP_(0) = ySSP;
    estimatedStateSSP_(1) = dySSP;

    return;
}

void Robbie::calculateDesPosSwingFoot() {
    static double startXFoot, startYFoot, startZFoot;
    double pxSWFoot, pySWFoot, pzSWFoot, newPxSWFoot, newPySWFoot, newPzSWFoot;
    static double nextX, nextY, nextZ, k_yx, k_zx;
    double vxBase, vxDes, pzBase;
    static Eigen::Matrix<double, 3, 1> nextDesiredSWFootPos;

    pxSWFoot = posSWingFoot_RobBIE_(0);
    pySWFoot = posSWingFoot_RobBIE_(1);
    pzSWFoot = posSWingFoot_RobBIE_(2);

    if(tau_<0.00){
        desPosSwingFoot_ << posSW2Base_world_;
    }else{
        Eigen::Matrix<double, 3, 1> beginPoint, endPoint, controlPoint1, controlPoint2;
        vxBase = velFloatBaseFil_(0);
        vxDes = cmdVelocity_;
        pzBase = posBase2ST_world_(2);

        if (vxBase * vxDes < 0)
            vxDes = 0;
        if (vxBase * vxBase < vxDes * vxDes)
            newPxSWFoot = vxDes > 0 ? -0.01 : 0.01;
        else if (vxBase < 0.)
            newPxSWFoot = -sqrt((vxBase * vxBase - vxDes * vxDes) * pzBase / gravity);
        else
            newPxSWFoot = sqrt((vxBase * vxBase - vxDes * vxDes) * pzBase / gravity);

        newPxSWFoot = velFloatBaseFil_(0)*stepTime/2 - kp_vx*(cmdVelocity_-velFloatBaseFil_(0));

        if (standingLeg_ == LEFTLEG) {
            float y_min = min(-estimatedStateSSP_(1)/sigma_ + kp_SWY*(estimatedStateSSP_(1)-estimatedDesDyStateSSP_), -estimatedStateSSP_(1)/lambda_);
            newPySWFoot = y_min;
//            newPySWFoot = -0.2;
        } else {
            float y_min = max(estimatedStateSSP_(1)/sigma_ - kp_SWY*(estimatedStateSSP_(1)-estimatedDesDyStateSSP_), estimatedStateSSP_(1)/lambda_);
            newPySWFoot = y_min;
//            newPySWFoot = 0.2;
        }
        newPzSWFoot = -cmdHeight_-0.02;

        beginPoint << posSW2Base_world_begein_;
        controlPoint1 << posSW2Base_world_begein_(0), posSW2Base_world_begein_(1), posSW2Base_world_begein_(2) + 0.1;
        controlPoint2 << newPxSWFoot, newPySWFoot, newPzSWFoot + 0.1;
        endPoint << newPxSWFoot, newPySWFoot, newPzSWFoot;
        desPosSwingFoot_ << bezierPoint3D(beginPoint, controlPoint1, controlPoint2, endPoint, tau_);
        desPosSwingFoot_(0) = limitRange(desPosSwingFoot_(0), posSWFootXLB, posSWFootXUB);

        if(standingLeg_==LEFTLEG)
            desPosSwingFoot_(1) = limitRange(desPosSwingFoot_(1), posSWFootYLB_Right, posSWFootYUB_Right);
        else
            desPosSwingFoot_(1) = limitRange(desPosSwingFoot_(1), posSWFootYLB_Left, posSWFootYUB_Left);

        /*if(standingLeg_==LEFTLEG){
            cout << "\nswing leg is right leg." << endl;
        }else{
            cout << "\nswing leg is left leg." << endl;
        }
        printf("start point: %f %f %f\n", beginPoint(0), beginPoint(1), beginPoint(2));
        printf("end point: %f %f %f\n", newPxSWFoot, newPySWFoot, newPzSWFoot);
        cout << tau_ << " tau, desired point: " << desPosSwingFoot_.transpose() << endl;*/
    }
    return;
}

/*void Robbie::calculateDesQSwingMotor() {

}*/
void Robbie::calculateDesQSwingMotor() {
/*/// v1.0 imu角度直接加上去算的， 不准
//gpt注释

        // v1.0 imu角度直接加上去算的， 不准
        // 定义局部变量用于存储腿部末端（foot）的期望位置，以及迭代过程中的中间变量
        double thetaX, thetaY, desiredXSW, desiredYSW, desiredZSW;
        desPosSwingFoot_RobBIE_ = rotationMatrix_.transpose() * desPosSwingFoot_;

        // 初始化迭代用的关节角度向量，位置向量，雅可比矩阵，以及期望位置向量
        Eigen::Matrix<double, 4, 1> qIter;
        Eigen::Matrix<double, 4, 1> pIter;
        Eigen::Matrix<double, 4, 4> JacobiIter;
        Eigen::Matrix<double, 4, 1> desP;

        // 容错和最大关节角度变化定义
        double permissionError = 0.01;
        double maxDq = 0.2;

        // 计算目标方向角度的调整值，限制在[-0.2, 0.2]范围内
        double dDirection = (cmdDirection - yawSTFoot_world_begin_) + (yawSWFoot_RobBIE_begin_ - yawSTFoot_RobBIE_begin_);
        if(dDirection < -0.2)
            dDirection = -0.2;
        if(dDirection > 0.2)
            dDirection = 0.2;

        // 根据时间系数tau计算实际使用的调整方向
        double s = tau_>1 ? 1 : tau_;
        desP << desPosSwingFoot_RobBIE_, yawSWFoot_RobBIE_begin_ + (dDirection) * s;
        cout << "desP: " << desP.transpose() << endl;

        // 根据当前摆动腿是左腿还是右腿，选择对应的电机数组并进行迭代求解
        if(standingLeg_ == LEFTLEG){
            Eigen::Matrix<double, 4, 1> dp, dq;
            int iterNum = 0;
            int maxIterNum = 10;
            double maxDpValue, maxDqValue;
            qIter = qMotorRight_;
            while (1){
                forward_kinematics_right(qIter, pIter); // 使用正运动学求解当前关节角度对应的足部位置
                dp = desP - pIter;
                maxDpValue = max(max(abs(dp(0)), abs(dp(1))), max(abs(dp(2)), abs(dp(3))));
                if(maxDpValue < permissionError || iterNum > maxIterNum)
                    break;
                else{
                    jacobianRight(qIter, JacobiIter); // 计算当前关节角度的雅可比矩阵
                    dq = JacobiIter.inverse() * dp;  // 使用雅可比逆矩阵计算角度增量
                    maxDqValue = max(max(abs(dq(0)), abs(dq(1))), max(abs(dq(2)), abs(dq(3))));
                    qIter += maxDq / maxDqValue * dq;  // 更新关节角度
                    ++iterNum;
                }
            }
        } else {
            // 与左腿类似，处理右腿
            Eigen::Matrix<double, 4, 1> dp, dq;
            int iterNum = 0;
            int maxIterNum = 100;
            double maxDpValue, maxDqValue;
            qIter = qMotorLeft_;
            while (1){
                forward_kinematics_left(qIter, pIter);
                cout << "forward_kinematics:" << pIter << endl;
                cout << "test:" << posSW2Base_RobBIE_ << endl;
                dp = desP - pIter;
                maxDpValue = max(max(abs(dp(0)), abs(dp(1))), max(abs(dp(2)), abs(dp(3))));
                if(maxDpValue < permissionError || iterNum > maxIterNum)
                    break;
                else{
                    jacobianLeft(qIter, JacobiIter);
                    dq = JacobiIter.inverse() * dp;
                    maxDqValue = max(max(abs(dq(0)), abs(dq(1))), max(abs(dq(2)), abs(dq(3))));
                    qIter += maxDq / maxDqValue * dq;
                    ++iterNum;
                }
            }
        }
        // 更新最终的期望关节角度
        desQSwingMotor_ << qIter;*/
    desQSwingMotor_ << 0,0,-0.5,1.1;
        return;
    }



void Robbie::controlSwingLeg() {

    calculateStateSSP();
    calculateDesPosSwingFoot();
    calculateDesQSwingMotor();

    if(standingLeg_ == LEFTLEG){
        if(isnan(desQSwingMotor_(1))||isnan(desQSwingMotor_(2))|| isnan(desQSwingMotor_(0))){
            desQSwingMotor_(0) = qMotorRight_(0);
            desQSwingMotor_(1) = qMotorRight_(1);
            desQSwingMotor_(2) = qMotorRight_(2);
        }

        desTorqueRightMotor_(0) = kp_Abd_Swing * (desQSwingMotor_(0) - qMotorRight_(0)) + (kd_Abd_Swing) * (0 - dqMotorRight_(0));
        desTorqueRightMotor_(1) = kp_Yaw_Swing * (desQSwingMotor_(1) - qMotorRight_(1)) + kd_Yaw_Standing * (0 - dqMotorRight_(1));
        desTorqueRightMotor_(2) = kp_Hip_Swing * (desQSwingMotor_(2) - qMotorRight_(2)) + kd_Hip_Swing * (0 - dqMotorRight_(2));
        desTorqueRightMotor_(3) = kp_Knee_Swing * (desQSwingMotor_(3) - qMotorRight_(3)) + kd_Knee_Swing * (0 - dqMotorRight_(3));

    }
    else{
        if(isnan(desQSwingMotor_(1))||isnan(desQSwingMotor_(2))|| isnan(desQSwingMotor_(0))){
            desQSwingMotor_(0) = qMotorLeft_(0);
            desQSwingMotor_(1) = qMotorLeft_(1);
            desQSwingMotor_(2) = qMotorLeft_(2);
        }
        desTorqueLeftMotor_(0) = kp_Abd_Swing * (desQSwingMotor_(0) - qMotorLeft_(0)) + kd_Abd_Swing * (0 - dqMotorLeft_(0));
        desTorqueLeftMotor_(1) = kp_Yaw_Swing * (desQSwingMotor_(1) - qMotorLeft_(1)) + kd_Yaw_Swing * (0 - dqMotorLeft_(1));
        desTorqueLeftMotor_(2) = kp_Hip_Swing * (desQSwingMotor_(2) - qMotorLeft_(2)) + kd_Hip_Swing * (0 - dqMotorLeft_(2));
        desTorqueLeftMotor_(3) = kp_Knee_Swing * (desQSwingMotor_(3) - qMotorLeft_(3)) + kd_Knee_Swing * (0 - dqMotorLeft_(3));
    }
}

vector<double> Robbie::controlMotor(vector<double>& control) {
    for (int i = 0; i < 4; ++i) {
        desTorqueLeftMotor_(i) = limitRange(desTorqueLeftMotor_(i), torqueLB, torqueUB);
        desTorqueRightMotor_(i) = limitRange(desTorqueRightMotor_(i), torqueLB, torqueUB);
    }
    if(control.size()!=8)
        control.resize(6);
    for (int i = 0; i < 4; ++i) {
        control[i] = desTorqueLeftMotor_(i);
        control[i+4] = desTorqueRightMotor_(i);
    }
    return control;
}
/*vector<double> Robbie::controlMotor(vector<double>& control) {

    if(control.size()!=8)
        control.resize(6);
    for (int i = 0; i < 4; ++i) {
        control[i] = 1;
        control[i+4] = 0;
    }
    return control;
}*/

void Robbie::logData(vector<double> &vectorData) {
    vectorData.clear();
    vectorData.push_back(t_);

    vectorData.push_back(qFloatBase_(0));
    vectorData.push_back(qFloatBase_(1));
    vectorData.push_back(qFloatBase_(2));
    vectorData.push_back(dqFloatBase_(0));
    vectorData.push_back(dqFloatBase_(1));
    vectorData.push_back(dqFloatBase_(2));

    vectorData.push_back(qMotorLeft_(0));
    vectorData.push_back(qMotorLeft_(1));
    vectorData.push_back(qMotorLeft_(2));
    vectorData.push_back(qMotorLeft_(3));

    vectorData.push_back(qMotorRight_(0));
    vectorData.push_back(qMotorRight_(1));
    vectorData.push_back(qMotorRight_(2));
    vectorData.push_back(qMotorRight_(3));

    vectorData.push_back(dqMotorLeft_(0));
    vectorData.push_back(dqMotorLeft_(1));
    vectorData.push_back(dqMotorLeft_(2));
    vectorData.push_back(dqMotorLeft_(3));

    vectorData.push_back(dqMotorRight_(0));
    vectorData.push_back(dqMotorRight_(1));
    vectorData.push_back(dqMotorRight_(2));
    vectorData.push_back(dqMotorRight_(3));

    vectorData.push_back(torqueMotorLeft_(0));
    vectorData.push_back(torqueMotorLeft_(1));
    vectorData.push_back(torqueMotorLeft_(2));
    vectorData.push_back(torqueMotorLeft_(3));

    vectorData.push_back(torqueMotorRight_(0));
    vectorData.push_back(torqueMotorRight_(1));
    vectorData.push_back(torqueMotorRight_(2));
    vectorData.push_back(torqueMotorRight_(3));

    vectorData.push_back(desQSwingMotor_(0));
    vectorData.push_back(desQSwingMotor_(1));
    vectorData.push_back(desQSwingMotor_(2));
    vectorData.push_back(desQSwingMotor_(3));

    vectorData.push_back(desTorqueLeftMotor_(0));
    vectorData.push_back(desTorqueLeftMotor_(1));
    vectorData.push_back(desTorqueLeftMotor_(2));
    vectorData.push_back(desTorqueLeftMotor_(3));

    vectorData.push_back(desTorqueRightMotor_(0));
    vectorData.push_back(desTorqueRightMotor_(1));
    vectorData.push_back(desTorqueRightMotor_(2));
    vectorData.push_back(desTorqueRightMotor_(3));

    vectorData.push_back(posSWingFoot_RobBIE_(0));
    vectorData.push_back(posSWingFoot_RobBIE_(1));
    vectorData.push_back(posSWingFoot_RobBIE_(2));

    vectorData.push_back(desPosSwingFoot_(0));
    vectorData.push_back(desPosSwingFoot_(1));
    vectorData.push_back(desPosSwingFoot_(2));

    vectorData.push_back(posBase2ST_RobBIE_(0));
    vectorData.push_back(posBase2ST_RobBIE_(1));
    vectorData.push_back(posBase2ST_RobBIE_(2));

    vectorData.push_back(velFloatBase_(0));
    vectorData.push_back(velFloatBase_(1));
    vectorData.push_back(velFloatBase_(2));

    vectorData.push_back(velFloatBaseFil_(0));
    vectorData.push_back(velFloatBaseFil_(1));
    vectorData.push_back(velFloatBaseFil_(2));

    vectorData.push_back(forceGroundLeft_(0));
    vectorData.push_back(forceGroundLeft_(1));
    vectorData.push_back(forceGroundLeft_(2));

    vectorData.push_back(forceGroundRight_(0));
    vectorData.push_back(forceGroundRight_(1));
    vectorData.push_back(forceGroundRight_(2));

    vectorData.push_back(forceGroundLeftFil_(0));
    vectorData.push_back(forceGroundLeftFil_(1));
    vectorData.push_back(forceGroundLeftFil_(2));

    vectorData.push_back(forceGroundRightFil_(0));
    vectorData.push_back(forceGroundRightFil_(1));
    vectorData.push_back(forceGroundRightFil_(2));

    vectorData.push_back(cmdHeight_);
    vectorData.push_back(cmdVelocity_);

    if(standingLeg_ == LEFTLEG){
        vectorData.push_back(-estimatedDesDyStateSSP_);
        vectorData.push_back(-estimatedStateSSP_(0));
        vectorData.push_back(-estimatedStateSSP_(1));
        vectorData.push_back(-estimatedStateSSP_(1)/sigma_);
    }else{
        vectorData.push_back(estimatedDesDyStateSSP_);
        vectorData.push_back(estimatedStateSSP_(0));
        vectorData.push_back(estimatedStateSSP_(1));
        vectorData.push_back(estimatedStateSSP_(1)/sigma_);
    }

    vectorData.push_back(pos_world(0));
    vectorData.push_back(pos_world(1));
    vectorData.push_back(pos_world(2));

    vectorData.push_back(vel_world(0));
    vectorData.push_back(vel_world(1));
    vectorData.push_back(vel_world(2));
    vectorData.push_back((double )standingLeg_);


}

void Robbie::getVelX(double* vel, double* velFil){
    *vel = velFloatBase_(0);
    *velFil = velFloatBaseFil_(0);
}