#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

float max_theta = 0.0;
float min_theta = 0.0;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) 
{
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
  // Predict the State
    x_ = F_ *x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_ ;
}

void KalmanFilter::Update(const VectorXd &z) {
  // Update Kalman Filter
      VectorXd y  = z - H_ * x_;
      MatrixXd Ht = H_.transpose();
      MatrixXd S  = H_ * P_ * Ht + R_;
      MatrixXd Si = S.inverse();
      MatrixXd K  = P_ * Ht * Si;

  // New State
      x_ = x_ + (K * y);
      int x_size = x_.size();
      MatrixXd I = MatrixXd::Identity(x_size, x_size);

      P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // Update the state by using the Extended Kalman Filter Equations
     
    // State Matrix x_
        float x  = x_(0);
        float y  = x_(1);
        float vx = x_(2);
        float vy = x_(3);

    // Cartesian to Polar
        float rho = sqrt(x*x + y*y);
        float theta = atan2(y, x);
        float rho_dot = (x*vx + y*vy) / rho;
      
    // Prediction
        VectorXd z_pred = VectorXd(3);
        z_pred << rho, theta, rho_dot;

        VectorXd Y = z - z_pred;

    // Normalize theta
        Y[1] = atan2( sin(Y[1]), cos(Y[1]));

  // Update Kalman Filter
      MatrixXd Ht = H_.transpose();
      MatrixXd S  = H_ * P_ * Ht + R_;
      MatrixXd Si = S.inverse();
      MatrixXd K  = P_ * Ht * Si;

  // New State
      int x_size = x_.size();
      MatrixXd I = MatrixXd::Identity(x_size, x_size);
     
      x_ = x_ + (K * Y);
      P_ = (I - K * H_) * P_;
}
