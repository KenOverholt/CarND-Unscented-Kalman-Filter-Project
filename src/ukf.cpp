#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  
  is_initialized_ = false;

  previous_timestamp_ = 0;
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;

  //set n_x_
  n_x_ = 5;  //KRO2 added
  
  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);  //KRO2 this works but should use n_x_, the state dimension

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2;	//KRO2: will need to modify this (original value: 30)

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;	//KRO2: will need to modify this (original value: 30)

  // Laser noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z_ = 3;
  
  //set augmented dimension
  n_aug_ = 7;	//KRO2 added

  //create matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //create matrix for cross correlation Tc_
  Tc_ = MatrixXd(n_x_, n_z_);

  //create matrix for z, incoming radar measurement
  //z_ = VectorXd(n_z_);
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  //define spreading parameter
  double lambda_ = 3 - n_aug_;

  //create augmented mean vector
  x_aug_ = VectorXd(n_aug_);	//KRO2 added

  //create augmented state covariance
  P_aug_ = MatrixXd(7, 7);	//KRO2 added

  //create sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);	//KRO2 added

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);  //KRO2 added
  
  //create vector for weights
  weights_ = VectorXd(2*n_aug_+1);  //KRO2 added
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  
  //KRO2: This entire method is from EKF so update/remove as necessary
  
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state x_ with the first measurement.
      * Create the covariance matrix.
    */
    // first measurement
    x_ << 1, 1, 0, 0, 0;  //KRO important for RMSE; first two will be overwritten but I should play with the last 3

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //measurement_pack.raw_measurements_(0) is rho
      //measurement_pack.raw_measurements_(1) is phi
      x_(0) = measurement_pack.raw_measurements_(0) * cos(measurement_pack.raw_measurements_(1)); 
      x_(1) = measurement_pack.raw_measurements_(0) * sin(measurement_pack.raw_measurements_(1));
      //set z_ from incoming radar measurement
      //KRO2 temporarily remove
	//      z_ <<  measurement_pack.raw_measurements_(0),
          //   measurement_pack.raw_measurements_(1),
            // measurement_pack.raw_measurements_(2);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //measurement_pack.raw_measurements_(0) is x
      //measurement_pack.raw_measurements_(1) is y

      x_(0) = measurement_pack.raw_measurements_(0);
      x_(1) = measurement_pack.raw_measurements_(1);
    }

    //KRO2 don't think this is needed.  Remove for now.
    //F_ << 1, 0, 0, 0,  //KRO 1 diagonal matrix
    //      0, 1, 0, 0,
    //      0, 0, 1, 0,
    //      0, 0, 0, 1;
    
    previous_timestamp_ = measurement_pack.timestamp_;
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for the Q matrix.
   */
  
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt expressed in seconds
	
  /*
//KRO2: Not sure if any of this is needed. It is from EKF.
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  
  //Modify the F matrix so that the time is integrated  L5, S8
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  // set the process covariance matrix Q as in L5, s9
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << noise_ax_*(dt_4/4), 0, noise_ax_*(dt_3/2), 0,
	     0, noise_ay_*(dt_4/4), 0, noise_ay_*(dt_3/2),
             noise_ax_*(dt_3/2), 0, noise_ax_*dt_2, 0,
             0, noise_ay_*(dt_3/2), 0, noise_ay_*dt_2;
  */
	
  Prediction(dt);
  previous_timestamp_ = measurement_pack.timestamp_;
	
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	  UpdateRadar(measurement_pack);
   /* KRO2: from EKF. don't think it is needed:
   // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_); // KRO set to Hj which should be set using the Jacobian function in tools.cpp
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    */
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    	UpdateLidar(measurement_pack);
	  /* KRO2: from EKF. don't think it is needed:
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
    */
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  ////////////////////////////////////////////////////////////////////////
  // Create Augmented Signa Points	//L7, sect 18
  ///////////////////////////////////////////////////////////////////////
  
  //create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug_;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)       = x_aug_ + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L.col(i);
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  // Predict Sigma Points   L7, sect 21
  ///////////////////////////////////////////////////////////////////////////////

    //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }  //end predict sigma points    for (int i = 0; i< 2*n_aug+1; i++)

  ////////////////////////////////////////////////////////////////////////////
  // Predict mean and covariance	//L7, sect 24
  ///////////////////////////////////////////////////////////////////////

    // set weights
  double weight_0 = lambda_/(lambda_ + n_aug_);
  weights_(0) = weight_0;
  double weight = 0.5/(n_aug_+lambda_);
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    weights_(i) = weight;
  }

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }  //end Predict mean and covariance

	
}  //end void UKF::Prediction(double delta_t) {

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  ////////////////////////////////////////////////////////////////////
  // Predict Radar Sigma Points   L7, sect 27
  ///////////////////////////////////////////////////////////
	
//cout << " n_aug_: " << n_aug_ << endl;
//cout << "Xsig_pred_: " << Xsig_pred_ << endl;
//cout << "Zsig_: " << Zsig_ << endl;
	
    //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
	  
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
	  
    // measurement model
    Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig_.col(i);  //KRO2 reset the weights_ var?
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();	  
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  //////////////////////////////////////////////////////////////////////////////
  // Update Radar  L7, sect 30
  //////////////////////////////////////////////////////////////////

    //calculate cross correlation matrix
//cout << "z_pred.size() " << z_pred.size() << endl;
//cout << "Tc_.size(): " << Tc_.size() << endl;
  Tc_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
//    cout << "i: " << i << endl;

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred;
//   cout << "z_diff.size(): " << z_diff.size() << endl;
	  
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
//   cout << "x_diff.size(): " << x_diff.size() << endl; 
   //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
//cout << "weights_.size(): " << weights_.size() << endl;
	  
    Tc_ = Tc_ + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc_ * S.inverse();

  //residual
//KRO2 temp comment out  VectorXd z_diff = z_ - z_pred;
VectorXd z_diff = z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

} // end void UKF::UpdateRadar(MeasurementPackage meas_package)
