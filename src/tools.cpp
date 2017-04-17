#include <iostream>
#include "tools.h"
#include <math.h>

#define RMSE_SIZE 4

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // define rmse vector
  VectorXd rmse(RMSE_SIZE);

  // initialize rmse vector with zeros
  rmse << 0,0,0,0;

  // check correctness of estimation and ground truth vector sizes
  if( estimations.size() == 0){
    std::cout <<" CalculateRMSE() - Error - Estimation vector size is Zero";
  }
  else if(estimations.size() != ground_truth.size()){
    std::cout <<" CalculateRMSE() - Error - Estimation vector size not equal to ground truth vector size";
  }
  else{
    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
      VectorXd est = estimations[i];
      VectorXd gnd_tth = ground_truth[i];
      VectorXd vec_diff = est-gnd_tth;
      VectorXd diff_sqrd = vec_diff.array()*vec_diff.array();
      rmse = rmse + diff_sqrd;
    }
    
    //calculate the mean
    rmse = rmse/estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();
  }

  // return the resulting vector
  return rmse;
  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  // define the Jacobian Matrix
  MatrixXd Hj(3,4);

  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // calculate some values to avoid recalculation
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  // check division by zero
  if(fabs(c1) < 0.0001){
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
  }
  else{
    //compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
         -(py/c1), (px/c1), 0, 0,
          py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  }

  return Hj;  
  
}
