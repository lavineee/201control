//
// Created by yeyinong on 2023/11/8.
//
#include "robbie/RobBIE.h"
#include "iostream"
#include "math.h"

using namespace std;



int main(){
    for (int i = 0; i < 35; ++i) {
        cout << sin(2*M_PI*((double )i/35.))<<endl;
    }
/*    float Robot_State[22] = {
            1,7.109e-21,0,0,
            2.975e-18,2.525e-18,-4.138e-18,
            -2.545e-20,0.384,-0.6981,-9.222e-20,0.384,-0.6981,
            -1.018e-17,1.754e-17,-5.357e-17,-3.689e-17,-2.314e-19,-9.161e-18,
            -0.9195,-0.393,0.35
    };
    float Robot_State2[22] = {
            0.9998, 0.01856,    -0.003203,  5.947e-05,
            -0.01517,-0.01448,0.005035,
            0.01562,0.367,-0.6941,-0.01312,0.3963,-0.7129,
            -0.002073,0.004924,-0.0001753,-0.01729,0.01931,0.001503,
            -0.9195,-0.393,0.35
    };
    float *ptr;
    ptr=RobBIE_onnxruntime(Robot_State);
    for(int i=0;i<6;i++)
    {
        std::cout<<ptr[i]<<" ";
    }
    std::cout << endl;

    ptr=RobBIE_onnxruntime(Robot_State2);
    for(int i=0;i<6;i++)
    {
        std::cout<<ptr[i]<<" ";
    }
    std::cout << endl;*/

    std::cout<<"this is mian"<<std::endl;
    getchar();
    return 0;
}
