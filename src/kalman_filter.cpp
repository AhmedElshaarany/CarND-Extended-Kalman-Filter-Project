#include "kalman_filter.h"
#include <math.h>
#include <iostream>

#define PI 3.14159265

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /*
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
  std::cout << "Start Lidar Update \n";
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */
  VectorXd z_pred(3);
  float px_pred = x_[0];
  float py_pred = x_[1];
  float vx_pred = x_[2];
  float vy_pred = x_[3];

  // calculate ro
  float ro = sqrt(px_pred*px_pred + py_pred*py_pred);

  // calcualte phi
  float phi = atan2(py_pred,px_pred);

  // calculate ro_dot
  float ro_dot = (px_pred*vx_pred + py_pred*vy_pred)/ro;
  
  // set the predicted state in polar coordinates
  z_pred << ro, phi, ro_dot;

  VectorXd y = z - z_pred;

  // check if phi is in the range -pi to pi, otherwise adjust
  if(y[1] < -PI){
    while(y[1] < -PI){
      y[1] += 2*PI;
    }
  }
  else if(y[1] > PI){
    while(y[1] > PI){
      y[1] -= 2*PI;
    }
  }

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}
