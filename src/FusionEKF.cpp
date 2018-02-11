#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  
  
		
  ekf_.Q_ = MatrixXd(4,4);
  noise_ax_ = 9;
  noise_ay_ = 9;
  
  H_laser_ << 1, 0, 0, 0,
      0, 1, 0, 0;
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.H_ = H_laser_;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    //ekf_.x_ << 1, 1, 1, 1;
	
	ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1000, 0,
		0, 0, 0, 1000;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	  float rho = measurement_pack.raw_measurements_[0];
	  float phi = measurement_pack.raw_measurements_[1];
	  float rhodot = measurement_pack.raw_measurements_[2];
	  ekf_.x_ << rho* cos(phi), rho*sin(phi), 0, 0;
	  
	  /*
	  if(fabs(rho* cos(phi)) < 0.0001){
                ekf_.x_[0] = 1;
                ekf_.P_(0,0) = 1000;
            }
	  if(fabs(rho*sin(phi)) < 0.0001){
				ekf_.x_[1] = 1;
				ekf_.P_(1,1) = 1000;
			}
	  */
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	  ekf_.x_ << measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1],0,0;
	  
    }
	previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
   previous_timestamp_ = measurement_pack.timestamp_;
   
   ekf_.F_ << 1, 0, dt, 0,
		  0, 1, 0, dt,
		  0, 0, 1, 0,
		  0, 0, 0, 1;
	
	float dt_4 = dt*dt*dt*dt;
    float dt_3 = dt*dt*dt;
    float dt_2 = dt*dt;

    ekf_.Q_ << dt_4/4*noise_ax_, 0, dt_3/2*noise_ax_, 0,
            0, dt_4/4*noise_ay_, 0, dt_3/2*noise_ay_,
            dt_3/2*noise_ax_, 0, dt_2*noise_ax_, 0,
            0, dt_3/2*noise_ay_, 0, dt_2*noise_ax_;

   ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	z_ = VectorXd(3);
	z_ << measurement_pack.raw_measurements_[0],
		measurement_pack.raw_measurements_[1],
		measurement_pack.raw_measurements_[2];
	ekf_.R_ = R_radar_;
    // Adding the RADAR measurement noise coefficient to our measurements
    //z_ += ekf_.R_.diagonal();
	
	ekf_.UpdateEKF(z_);
  } else {
    // Laser updates
	z_ = VectorXd(2);
	z_<< measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1];
	
	ekf_.R_ = R_laser_;
    // Adding the LASERmeasurement noise coefficient to our measurements
	ekf_.Update(z_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
