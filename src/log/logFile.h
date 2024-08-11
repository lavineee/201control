//
// Created by leggedrobot on 23-7-9.
//

#ifndef BIPEDALNUCV3YM_LOGFILE_H
#define BIPEDALNUCV3YM_LOGFILE_H

#include "robbie/Robbie.h"
#include <mutex>

using namespace std;
extern mutex mutex_robbie;

void recordData(Robbie* robbie, mjData *d);

#endif //BIPEDALNUCV3YM_LOGFILE_H
