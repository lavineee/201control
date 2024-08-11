//
// Created by yeyinong on 2023/3/30.
//

#ifndef MUJOCOTEST3D_CONTROLLER_H
#define MUJOCOTEST3D_CONTROLLER_H

#include "robbie/Robbie.h"
#include "mujoco/mujoco.h"
#include "mutex"
#include "thread"
extern Robbie* robbie;
extern mutex mutex_robbie;
extern thread* record;
void init_controller(const mjModel* m, mjData* d);

void controllerHandler(const mjModel *m, mjData *d);

#endif //MUJOCOTEST3D_CONTROLLER_H
