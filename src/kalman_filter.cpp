#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

const double EPS = 0.001;

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
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  
  // calculate matrices
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

double normalizeAngle_(double angle) {
  if (angle > M_PI) {
    while (angle < M_PI) {
      angle -= 2.0 * M_PI;
    }
  } else if (angle < -M_PI) {
    while (angle > -M_PI) {
      angle += 2.0 * M_PI;
    }
  }
  return angle;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  VectorXd z_pred(3);
  double rho_pred = std::sqrt(x_[0] * x_[0] + x_[1] * x_[1]);
  double phi_pred = 0.0;
  if (rho_pred < EPS) {
	  rho_pred = EPS;
  } else {
	  phi_pred = std::atan2(x_[1], x_[0]);
  }
  phi_pred = normalizeAngle_(phi_pred);
  double rhodot_pred = (x_[0] * x_[2] + x_[1] * x_[3]) / rho_pred;
  z_pred << rho_pred, phi_pred, rhodot_pred;

  // debugging
  z_ = z_pred;
  
  // calculate matrices
  VectorXd y = z - z_pred;

  // normalise the angle value of y
  y[1] = normalizeAngle_(y[1]);

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
