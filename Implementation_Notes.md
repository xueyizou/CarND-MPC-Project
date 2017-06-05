# CarND-Controls-MPC Implementation Notes
Self-Driving Car Engineer Nanodegree Program.

by Xueyi Zou


## The Model
The model used is the same as that disscussed in the course.

The state include:
* `px`: x coordinate (unit: m);
* `py`: y coordinate (unit: m);
* `psi`: heading angle (unit: rad);
* `v`: speed (unit: m/s);
* `cte`: cross track error (unit m);
* `epsi`: heading angle error (unit rad);

Actuators are:
* `delta`: steering angle (unit: rad);
* `a`: acceleration (unit: m/s^2);

The update equations are as follows:
* px' = px + v \* cos(psi) \* dt;
* py' = px + v \* sin(psi) \* dt;
* psi' = psi + v / Lf\* delta  \* dt;
* v' = v + a \* dt;
* cte = coeffs[0]+coeffs[1] \* px + coeffs[2] \* px^2 + coeffs[3] \* px^3 - py;
* epsi = atan(coeffs[1] + 2\*coeffs[2] \* px + 3\*coeffs[3]\* px^3 - psi.


## Timestep Length and Elapsed Duration (N & dt)
The final timestep length (`N`) and elapsed duration (`dt`) I set is `N=30` and `dt=0.05`.

The reason is that in the MPC I want to predict a total of about **2 seconds** ahead, and to emulate the continuous vehicle motion with discrete steps, the resolution (i.e. `dt`) should be small. However, if `dt` is set to small, `N` should be large (because the total prediction time t = N*dt). This will result in a high computation cost and make the controller fail to be real-time. So, finally, I found `N=30` and `dt=0.05` is proper. That being said, I have also tried N=20, 25, 50 and dt=0.01, 0.02, 0.1. Their performances are not as good as `N=30` and `dt=0.05`.



## Polynomial Fitting and MPC Preprocessing

I first transformed the coordinates of the currently received waypoints to be relative to the current vehicle pose. Then I fitted a 3rd polynomial using the function `Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order)`. The code is located in lines 120~130 of `main.cpp`.




## Model Predictive Control with Latency

The idea to overcome control latency is to predict the state of the car after a delay. In this case, I predict the state after 100 millisecond, and treat that state as the initial state for the MPC. The core code to overcome latency is located in lines 114-117 of `main.cpp`. The model used to predict future state is the same as that disscussed above.
