#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  const double pi_ = 2*acos(0.0);

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* time when the state is true, in us
  long long previous_timestamp_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Radar measurment dimension
  int n_r_z_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* Augmented Sigma point spreading parameter
  double lambda_aug_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* sigma points matrix
  MatrixXd Xsig_;

  // augmented sigma point matrix
  MatrixXd Xsig_aug_;

  // predicted sigma points as columns
  MatrixXd Xsig_pred_;

  // sigma points in measurement space
  MatrixXd Zsig_;

  //mean predicted measurement
  MatrixXd z_pred_;

  //measurement covariance matrix S
  MatrixXd S_;

  //measurement noise covariance matrix
  MatrixXd R_;

  //measurement covariance matrix - laser
  MatrixXd R_laser_;

  //measurement matrix - laser
  MatrixXd H_laser_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(const VectorXd &z);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const VectorXd &z);

  void GenerateSigmaPoints();
  void AugmentedSigmaPoints();
  void SigmaPointPrediction(double delta_t);
  void PredictMeanAndCovariance();
  void PredictRadarMeasurement();
  void UpdateState(const VectorXd &z);

};

#endif /* UKF_H */
