//
// Created by yeyinong on 2023/3/30.
//
#include "controller.h"
#include "vector"
#include "robbie/Robbie.h"
#include "eigen3/Eigen/Geometry"
#include "thread"
#include "log/logFile.h"
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;


void init_controller(const mjModel* m, mjData* d){
    const int qNum = 18;

//    vector<double> qPos{
//            0,0,0,0,0,0,
//            0.,0.,-0.6,1.4,0.8,-1.4,
//            0.,0.,-0.6,1.4,0.8,-1.4
//    };
    vector<double> qPos{
            0,0,0,0,0,0,
            0,0,-0.5,1.1,-0.66,
            0,0,-0.5,1.1,-0.66
           /* 0,0,0,0,0,0,
            0.,-3./180*M_PI,-0.5,0.9,0.4,-0.9,
            0.,3./180*M_PI,-0.5,0.9,0.4,-0.9*/
    };
    vector<double> qVel{
            0,0,0,0,0,0,
            -0.,0.,-0.,0.,0,0,
            0.,0.,-0.,0.,0,0
    };

    for (int i = 0; i < qNum; ++i) {
        d->qpos[i] = qPos[i];
    }
    for (int i = 0; i < qNum; ++i) {
        d->qvel[i] = qVel[i];
    }
    //d->qvel[1] = -0.2;
}

void controllerHandler(const mjModel *m, mjData *d){
    static bool hasInited = false;
    vector<double> floatMessage;
    vector<double> motorMessage;
    vector<double> control;

    static double lastTime=0;
    static double time=0;

    time = d->time;
    if(time==0)
        return;

    getSensorMessage(d, floatMessage, motorMessage);

    if(!hasInited){
        robbie = new Robbie(RIGHTLEG, motorMessage, floatMessage);
        record = new thread(recordData, robbie, d);
        hasInited = true;
    }

    mutex_robbie.lock();
    robbie->update(motorMessage, floatMessage, time-lastTime);
    robbie->controlMotor(control);
    mutex_robbie.unlock();

    lastTime = time;
    for (int i = 0; i < 8; ++i) {
        d->ctrl[i] = control[i];
    }
    cout << "\n time: " << d->time << endl;

}
