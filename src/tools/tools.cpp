//
// Created by yeyinong on 2023/3/30.
//
#include "tools.h"
vector<double> bezierPoint(double px1, double py1, double px2, double py2, double px3, double py3, double px4, double py4, double tau){
    if(tau>=1)
        tau = 1;
    if(tau <= 0)
        tau =  0;
    double px = px1*(pow(1-tau,3)) + 3*px2*tau*(pow(1-tau,2)) + 3*px3* pow(tau, 2)*(1-tau) + px4* pow(tau, 3);
    double py = py1*(pow(1-tau,3)) + 3*py2*tau*(pow(1-tau,2)) + 3*py3* pow(tau, 2)*(1-tau) + py4* pow(tau, 3);
//    res.x -= t*(src[3].x-src[0].x)/2;
    vector<double> res = vector<double>{px, py};
    return res;
}

Eigen::Matrix<double, 3, 1> bezierPoint3D(Eigen::Matrix<double, 3, 1> startPoint, Eigen::Matrix<double, 3, 1> controlPoint, Eigen::Matrix<double, 3, 1> endPoint,  double tau){
    if(tau>=1)
        tau = 1;
    if(tau <= 0)
        tau =  0;

    double px = startPoint(0)*(pow(1-tau,2)) + 2*controlPoint(0)*tau*(1-tau)  + endPoint(0) * pow(tau, 2);
    double py = startPoint(1)*(pow(1-tau,2)) + 2*controlPoint(1)*tau*(1-tau)  + endPoint(1) * pow(tau, 2);
    double pz = startPoint(2)*(pow(1-tau,2)) + 2*controlPoint(2)*tau*(1-tau)  + endPoint(2) * pow(tau, 2);
//    res.x -= t*(src[3].x-src[0].x)/2;
    Eigen::Matrix<double, 3, 1> point;
    point << px, py, pz;
    return point;
}

Eigen::Matrix<double, 3, 1> bezierPoint3D(Eigen::Matrix<double, 3, 1> startPoint, Eigen::Matrix<double, 3, 1> controlPoint1, Eigen::Matrix<double, 3, 1> controlPoint2, Eigen::Matrix<double, 3, 1> endPoint,  double tau){
    tau*=1.1;
    if(tau>=1)
        tau = 1;
    if(tau <= 0)
        tau =  0;

    double px = startPoint(0)*(pow(1-tau,3)) + 3*controlPoint1(0)*tau*(pow(1-tau,2)) + 3*controlPoint2(0)* pow(tau, 2)*(1-tau) + endPoint(0) * pow(tau, 3);
    double py = startPoint(1)*(pow(1-tau,3)) + 3*controlPoint1(1)*tau*(pow(1-tau,2)) + 3*controlPoint2(1)* pow(tau, 2)*(1-tau) + endPoint(1) * pow(tau, 3);
    double pz = startPoint(2)*(pow(1-tau,3)) + 3*controlPoint1(2)*tau*(pow(1-tau,2)) + 3*controlPoint2(2)* pow(tau, 2)*(1-tau) + endPoint(2) * pow(tau, 3);
//    res.x -= t*(src[3].x-src[0].x)/2;
    Eigen::Matrix<double, 3, 1> point;
    point << px, py, pz;
    return point;
}

void getSensorMessage(const mjData *d, vector<double>& floatMessage, vector<double>& motorMessage){
    int motorNum = 8;
    if(floatMessage.size()!=10)
        floatMessage.resize(10);
    if(motorMessage.size()!=motorNum*3)
        motorMessage.resize(motorNum*3);

    // 获取imu数据
    for (int i = 0; i < 10; ++i) {
        floatMessage[i] = d->sensordata[i];
    }

    // 获取电机的数据
    int sensorIndex_actuatorPos = 10;
    for (int i = 0; i < motorNum*3; ++i) {
        motorMessage[i] = d->sensordata[sensorIndex_actuatorPos+i];
    }

    int sensorPos = 34;
    Eigen::Matrix<double, 3, 1> posBase, posLeftFoot, posRightFoot;
    for (int i = 0; i < 3; ++i) {
        posBase(i) = d->sensordata[sensorPos + i];
        posLeftFoot(i) = d->sensordata[sensorPos + 3 + i];
        posRightFoot(i) = d->sensordata[sensorPos + 6 + i];
    }
    cout << "pos_leftFoot2Base_world_imu: " << (posLeftFoot-posBase).transpose() << endl;
    cout << "pos_rightFoot2Base_world_imu: " << (posRightFoot-posBase).transpose() << endl;

}
