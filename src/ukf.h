#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {

private:

  ///* time when the state is true, in us
  long long time_us_;

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

  ///* Number of sigma points
  int n_sigma_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* matrix of predicted state sigma points
  MatrixXd Xsig_;

  ///* mean of predicted state sigma points
  MatrixXd Xsig_pred_;

  ///* mean and a state sigma point difference
  MatrixXd Xdiff_;

  ///* matrix of predicted measurement sigma points
  MatrixXd Zsig_;

  ///* mean of predicted measurement sigma points
  VectorXd Zsig_pred_;

  ///* mean and a measurement sigma point difference
  MatrixXd Zdiff_;

  ///* measurement covariance matrix
  MatrixXd S_;

public:

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x;

  ///* state covariance matrix
  MatrixXd P;

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar;

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
  void Prediction(double delta_t_);

  /**
   * Updates the state and the state covariance matrix using a sensor measurement
   * @param meas_package The measurement at k+1
   */
  void Update(MeasurementPackage meas_package);

  void InitiateMeasurement(MeasurementPackage meas_package);

  void AugmentedSigmaPoints();

  void SigmaPointPrediction(double delta_t_);

  void PredictMeanAndCovariance();

  void PredictMeasurement(MeasurementPackage meas_package);

  void UpdateState(MeasurementPackage meas_package) ;

};

#endif /* UKF_H */
