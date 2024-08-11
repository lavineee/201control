//
// Created by leggedrobot on 23-7-9.
//

//
// Created by leggedrobot on 23-7-9.
//

#include "logFile.h"

#include "iostream"
#include "fstream"
#include "fstream"
#include <sys/stat.h>
#include "unistd.h"
#include "stdio.h"
#include "sys/time.h"
#include "vector"
#include "robbie/Robbie.h"

string listNames = "time,"
                   "rollTorsoAct,pitchTorsoAct,yawTorsoAct,wXAct,wYAct,wZAct,"
                   "posLYawAct,posLAbdAct,posLHipAct,posLKneeAct,posRYawAct,posRAbdAct,posRHipAct,posRKneeAct,"
                   "velLYawAct,velLAbdAct,velLHipAct,velLKneeAct,velRYawAct,velRAbdAct,velRHipAct,velRKneeAct,"
                   "torLYawAct,torLAbdAct,torLHipAct,torLKneeAct,torRYawAct,torRAbdAct,torRHipAct,torRKneeAct,"
                   ""
                   "posSWYawDes,posSWAbdDes,posSWHipDes,posSWKneeDes,"
                   "torLYawDes,torLAbdDes,torLHipDes,torLKneeDes,torRYawDes,torRAbdDes,torRHipDes,torRKneeDes,"
                   ""
                   "posSWFootXAct,posSWFootYAct,posSWFootZAct,"
                   "posSWFootXDes,posSWFootYDes,posSWFootZDes,"
                   ""
                   "posXAct,posYAct,posZAct,"
                   "velXAct,velYAct,velZAct,velXFil,velYFil,velZFil,"

                   "forceGroundLeftXAct,forceGroundLeftYAct,forceGroundLeftZAct,forceGroundRightXAct,forceGroundRightYAct,forceGroundRightZAct,"
                   "forceGroundLeftXFil,forceGroundLeftYFil,forceGroundLeftZFil,forceGroundRightXFil,forceGroundRightYFil,forceGroundRightZFil,"
                   ""
                   "posZCmd,velXCmd,"
                   ""
                   "dySSPDes,ySSPEst,dySSPEst,yFootSSPEst,"
                   ""
                   "pos_x_est,pos_y_est,pos_z_est,"
                   "vel_x_est,vel_y_est,vel_z_est,"
                   ""
                   "standingLeg,"
                   "\n";

void recordData(Robbie* robbie, mjData *d){
    vector<double> logData;
    const double duration = 1250;

    ofstream logFile_;
    char DirPath[64];
    time_t t = time(0);
    /*根据当前时间建立文件夹*/
    strftime(DirPath, sizeof(DirPath), "/home/yeyinong/RobBIE_Data/RobBIEv2/Data-%Y-%m-%d %H-%M-%S", localtime(&t));
    int isCreate = mkdir(DirPath, S_IRWXU  | S_IRWXG | S_IRWXO);
    if (!isCreate)
        printf("create path:%s\n", DirPath);
    else
        printf("create path failed! error\n");
    /*在文件夹下新建解析过的csv文件*/
    string filePath = string(DirPath) + "/Data_Csv_file.csv";
    logFile_.open(filePath);
    if(!logFile_.is_open()){
        cout << "cmd.txt open failed!" <<endl;
    }
    logFile_.tie(nullptr);

    struct timeval time1;
    struct timeval lastTime;
    gettimeofday(&lastTime, NULL);
    logFile_ << listNames;
    while(1){
        if (robbie == nullptr){
            usleep(1000);
            continue;
        }
        mutex_robbie.lock();
        robbie->logData(logData);
        mutex_robbie.unlock();
        for (double data : logData) {
            logFile_ << data << ",";
        }
        logFile_ << endl;

        gettimeofday(&time1, NULL);
        double deltaT = (time1.tv_sec-lastTime.tv_sec)*1000000 + (time1.tv_usec-lastTime.tv_usec);
        lastTime = time1;
        if (deltaT < duration)
            usleep(duration-deltaT);
    }
};