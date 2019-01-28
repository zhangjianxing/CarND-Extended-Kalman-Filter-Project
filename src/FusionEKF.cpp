#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // create a 4D state vector, we don't know yet the values of the x state
  ekf_.x_ = VectorXd(4);

  // state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;

  // the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  // set the acceleration noise components
  iterate = 0;
  noise_ax = 9;
  noise_ay = 9;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    cout << "Kalman Filter Initialization " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rho_dot = measurement_pack.raw_measurements_[2];
      double cos_phi = cos(phi);
      double sin_phi = sin(phi);

      double x = rho * cos_phi;
      double y = rho * sin_phi;
      double vx = rho_dot * cos_phi;
      double vy = rho_dot * sin_phi;
      if ( x < 0.0001 ) {
        x = 0.0001;
      }
      if ( y < 0.0001 ) {
        y = 0.0001;
      }
      ekf_.x_ << x, y, vx , vy;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_[0], 
                 measurement_pack.raw_measurements_[1], 
                 0, 
                 0;
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */
  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
              0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
              dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
              0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  // predict
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */
//  VectorXd measurement_position(2);
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    // measurement matrix
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);

    // measurement covariance
    ekf_.R_ = MatrixXd(3, 3);
    ekf_.R_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

//    double rho = raw_measurements_[0];
//    double phi = raw_measurements_[1];
//    double rho_dot = raw_measurements_[2];
//    double cos_phi = cos(phi);
//    double sin_phi = sin(phi);
//
//    double x = rho * cos_phi;
//    double y = rho * sin_phi;
//    measurement_position << x, y;
  } else {
    // TODO: Laser updates
      // measurement matrix
    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ << 1, 0, 0, 0,
               0, 1, 0, 0;

    // measurement covariance
    ekf_.R_ = MatrixXd(2, 2);
    ekf_.R_ << 0.0225, 0,
               0, 0.0225;
    ekf_.Update(measurement_pack.raw_measurements_);
//    measurement_position << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];
  }

  // print the output
//  iterate ++;
//  cout << endl << "iterate: " << iterate << endl;
//  cout << "dt = " << dt << endl;
//  cout << "measurement_pack.timestamp_ = " << measurement_pack.timestamp_ << endl;
//  cout << "measurement_pack.raw_measurements_ = " << endl << measurement_pack.raw_measurements_ << endl;
//  cout << "x_ = " << endl << ekf_.x_ << endl;
//  cout << "measurement_position = " << endl << measurement_position << endl;
//  cout << "P_ = " << endl << ekf_.P_ << endl;
}
