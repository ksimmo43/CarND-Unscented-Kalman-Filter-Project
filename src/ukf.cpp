#include "ukf.h"
#include "tools.h"
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

  pi_;

  is_initialized_ = false;

  previous_timestamp_ = 0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.55;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Radar measurment dimension
  n_r_z_ = 3;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  ///* Augmented Sigma point spreading parameter
  lambda_aug_ = 3 - n_aug_;

  ///* the current NIS for radar
  NIS_radar_ = 0;

  ///* the current NIS for laser
  NIS_laser_ = 0;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 2, 0, 0,
        0, 0, 0, .02, 0,
        0, 0, 0, 0, .5;

  ///* Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_aug_/(lambda_aug_+n_aug_);
  weights_.fill(0.5/(lambda_aug_+n_aug_));
  weights_(0) = weight_0;

  //create sigma point matrix (NOT USED)
  Xsig_ = MatrixXd(n_x_, 2 * n_x_ + 1);

  //create augmented sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //create matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_r_z_, 2 * n_aug_ + 1);

  //mean predicted measurement
  z_pred_ = VectorXd(n_r_z_);

  //measurement covariance matrix S
  S_ = MatrixXd(n_r_z_,n_r_z_);

  //measurement noise covariance matrix - radar
  R_ = MatrixXd(n_r_z_,n_r_z_);
  R_ << std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  //measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;

  //measurement matrix - laser
  H_laser_ = MatrixXd(2, n_x_);
  H_laser_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho_i = measurement_pack.raw_measurements_(0);
      float phi_i = measurement_pack.raw_measurements_(1);
      float rho_dot_i = measurement_pack.raw_measurements_(2);
      // Make sure phi is  -pi < phi < pi
      while (phi_i> M_PI) phi_i-=2.*M_PI;
      while (phi_i<-M_PI) phi_i+=2.*M_PI;

      x_ << rho_i * cos(phi_i), rho_i * sin(phi_i), ((rho_dot_i > 0) - (rho_dot_i < 0))*sqrt(pow(rho_dot_i * cos(phi_i),2) + pow(rho_dot_i * sin(phi_i),2)), 0, 0;
      previous_timestamp_ = measurement_pack.timestamp_;
      // std::cout << " x_ = " << x_ << std::endl;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0, 0;

      previous_timestamp_ = measurement_pack.timestamp_;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  if((measurement_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_) || (measurement_pack.sensor_type_ == MeasurementPackage::LASER && use_laser_)) {
     //compute the time elapsed between the current and previous measurements
    float delta_t = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;
    // std::cout << " delta_t = " << delta_t << std::endl;
    Prediction(delta_t);

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Radar updates

      //measurement update
    	UpdateRadar(measurement_pack.raw_measurements_);
    } else {
      // Laser updates
      //measurement update
    	UpdateLidar(measurement_pack.raw_measurements_);
    }
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  AugmentedSigmaPoints();
  SigmaPointPrediction(delta_t);
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state coriance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const VectorXd &z) {
  //VectorXd z_pred_ = H_laser_ * x_;
	VectorXd y = z - H_laser_ * x_;
	//MatrixXd Ht = H_.transpose();
	MatrixXd S = H_laser_ * P_ * H_laser_.transpose() + R_laser_;
	//MatrixXd Si = S.inverse();
	//MatrixXd PHt = P_ * H_laser_.transpose();
	MatrixXd K = P_ * H_laser_.transpose() * S.inverse();

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_laser_) * P_;

  // Calculate NIS
  NIS_laser_ = y.transpose()*S.inverse()*y;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const VectorXd &z) {
  PredictRadarMeasurement();
  UpdateState(z);
}

void UKF::GenerateSigmaPoints() {

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //set first column of sigma point matrix
   Xsig_.col(0)  = x_;

   //set remaining sigma points
   for (int i = 0; i < n_x_; i++)
   {
     Xsig_.col(i+1)     = x_ + sqrt(lambda_+n_x_) * A.col(i);
     Xsig_.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
   }

}
void UKF::AugmentedSigmaPoints() {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  x_aug << x_, 0, 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)       = x_aug + sqrt(lambda_aug_+n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_aug_+n_aug_) * L.col(i);
  }

}
void UKF::SigmaPointPrediction(double delta_t) {

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
  }

}
void UKF::PredictMeanAndCovariance(){
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose() ;
  }
}
void UKF::PredictRadarMeasurement() {

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
    if (Zsig_(0,i) < 0.00001) Zsig_(0,i) = 0.001;
    Zsig_(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / Zsig_(0,i);   //r_dot
  }

  //mean predicted measurement
  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred_ += weights_(i) * Zsig_.col(i);
  }

  //measurement covariance matrix S
  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S_ += weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S_ += R_;

}
void UKF::UpdateState(const VectorXd &z) {
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_r_z_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_.inverse();

  //residual
  VectorXd z_diff = z - z_pred_;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K*S_*K.transpose();

  // Calculate NIS
  NIS_radar_ = z_diff.transpose()*S_.inverse()*z_diff;

}
