# Robbie: Lightweight Bipedal Walking Control

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Platform](https://img.shields.io/badge/platform-MuJoCo%20%7C%20C++-green)

## 📌 Project Overview
Robbie is a lightweight **bipedal walking controller** implemented in C++.  
The project demonstrates how **analytical models** (Capture Point, LIP dynamics) can be combined with **virtual force control** to achieve stable and robust humanoid walking **without complex trajectory optimization (e.g., HZD/MPC)**.

**Key Features:**
- Sagittal plane: **Capture Point (Raibert heuristic)** for velocity tracking
- Coronal plane: **Linear Inverted Pendulum (LIP)** with step-to-step velocity convergence
- Stance leg: **Virtual force control (JᵀF)** for body attitude and height stabilization
- Swing leg: **Bezier trajectory + geometric IK + joint PD**
- Robust engineering details: GRF-based contact detection, velocity filtering, torque/position saturation

---

## 🚀 Control Framework
The walking control loop is composed of four main modules:

1. **State Estimation**  
   - Floating base position/velocity estimation (using stance leg kinematics)  
   - Contact detection via ground reaction force (GRF)  

2. **Stance Leg Control (Virtual Force + JᵀF)**  
   - Stabilizes roll/pitch and height  
   - Maps virtual force to joint torques via stance leg Jacobian  

   $$
   \boldsymbol{\tau}_{stance} = J_p^\top \mathbf{F}
   $$

3. **Swing Leg Control (Bezier + IK + PD)**  
   - Smooth foot trajectory generation between step start and desired foot placement  
   - Geometric inverse kinematics + joint PD for torque output  

   $$
   \tau_{swing} = K_q (q^* - q) + D_q (\dot q^* - \dot q)
   $$

4. **Foot Placement Strategy**  
   - **Sagittal (x): Raibert heuristic**
     $$
     x_{foot} = \tfrac{T}{2}\dot x + b_x + k_v(\dot x - \dot x^*)
     $$
   - **Coronal (y): LIP + feedback**
     $$
     y_{foot} = -s \hat y_{td} + k_y(\hat{\dot y}_{td} - \dot y^*_{td}) + b_y^{(s)}
     $$

---

## 📂 Repository Structure
