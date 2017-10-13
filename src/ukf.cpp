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

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(1);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5)

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Number of sigma points
  n_sigma_ = 2*n_aug_ + 1;

  // Weights of sigma points
  weights_ = VectorXd(n_sigma_);

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  //predicted sigma points state matrix
  Xsig_ = MatrixXd(n_aug_, n_sigma_);
  Xsig_pred_ = VectorXd(n_x_);

  // predicted sigma points radar measurement matrix
  Z_rad_sig_ = MatrixXd(3, n_sigma_);
  Z_rad_sig_pred_ = Vector(3);

  // predicted sigma points lidar measurement matrix
  Z_las_sig_ = MatrixXd(3, n_sigma_);
  Z_las_sig_pred_ = Vector(2);

  // Radar measurement covariance matrix
  S_rad = MatrixXd(3,3);

  // lidar measurement covariance matrix
  S_las = MatrixXd(2,2);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    /*****************************************************************************
   *  Initialization
   ****************************************************************************/
    if (!is_initialized_) {
    /**
      * Initialize the state x_ with the first measurement.
      */
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

            float ro     = measurement_pack.raw_measurements_(0);
            float phi    = measurement_pack.raw_measurements_(1);
            float ro_dot = measurement_pack.raw_measurements_(2);

            x_(0) = ro * cos(phi);
            x_(1) = ro * sin(phi);

            }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

            x_(0) = measurement_pack.raw_measurements_(0);
            x_(1) = measurement_pack.raw_measurements_(1);

            }
        previous_timestamp_ = measurement_pack.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
}

    float delta_t = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;

    AugmentedSigmaPoints();
    SigmaPointPrediction(delta_t);
    PredictMeanAndCovariance();

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

}

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

    MatrixXd Tc = MatrixXd(n_x_, 2);
    //calculate cross correlation matrix
    Tc.fill(0);
    MatrixXd diff_z = MatrixXd(2, 2);
    MatrixXd diff_x = MatrixXd(n_x_, n_x_);

    for (int i = 0; i < n_sigma_; i++) {
        diff_z = Z_las_sig_.col(i) - Z_las_sig_pred_;
        diff_x = Xsig_.col(i) - x_;
        Tc = Tc + weights(i) * (diff_x * diff_z.transpose());
    }

    //calculate Kalman gain K;
    MatrixXd K = MatrixXd(n_x, 2);
    K = Tc * S_las.inverse();

    //update state mean and covariance matrix
    VectorXd z = meas_package.raw_measurements_;
    x_ += K * (z - Z_las_sig_pred_);
    P_ -= K * S_las * K.transpose();
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
    MatrixXd Tc = MatrixXd(n_x_, 3);
    //calculate cross correlation matrix
    Tc.fill(0);
    MatrixXd diff_z = MatrixXd(3, 3);
    MatrixXd diff_x = MatrixXd(n_x_, n_x_);

    for (int i = 0; i < n_sigma_; i++) {
      diff_z = Z_rad_sig_.col(i) - Z_rad_sig_pred_;
      if (diff_z(1) < -M_PI) {
          diff_z(1) += 2*M_PI;
      }
      else if (diff_z(1) > M_PI) {
          diff_z(1) -= 2*M_PI;
      }
      diff_x = Xsig_.col(i) - x_;
      Tc = Tc + weights(i) * (diff_x * diff_z.transpose());
    }

    //calculate Kalman gain K;
    MatrixXd K = MatrixXd(n_x, 3);
    K = Tc * S_rad.inverse();

    //update state mean and covariance matrix
    VectorXd z = meas_package.raw_measurements_;
    x_ += K * (z - Z_rad_sig_pred_);
    P_ -= K * S_rad * K.transpose();
}

void UKF::AugmentedSigmaPoints() {
    matrixXd Xaug_ = matrixXd(n_aug_, n_aug_);
    Xaug_.fill(0);
    Xaug_.head(n_x_) = x_;

    Paug_.fill(0.0);
    Paug_.topLeftCorner(n_x_, n_x_) = P_;
    Paug_(5, 5) = std_a_;
    Paug_(6, 6) = std_yawdd_;

    MatrixXd A_ = P_.llt().matrixL();

    Xsig_.col(0) = Xaug_;
    for (int i = 0; i < n_aug; i++)
    {
        Xsig_.col(i+1)     = Xaug_ + sqrt(lambda+n_aug) * A.col(i);
        Xsig_.col(i+1+n_aug) = Xaug_ - sqrt(lambda+n_aug) * A.col(i);
    }
}

void UKF::SigmaPointPrediction(delta_t) {
    double px, py, v, yaw, yawd;
    double nu_a, nu_yawdd;

    double ppx, ppy, pv, pyaw, pyawd;

    for (int i = 0; i < n_sigma_; i++) {
        px   = Xsig_(0, i);
        py   = Xsig_(1, i);
        v    = Xsig_(2, i);
        yaw  = Xsig_(3, i);
        yawd = Xsig_(4, i);

        nu_a    = Xsig_(5, i);
        nu_yawdd = Xsig_(6, i);

        if (yawd<0.0001) {
            ppx = px + v*cos(yaw)*delta_t + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
            ppy = py + v*sin(yaw)*delta_t + 0.5*delta_t*delta_t*sin(yaw)*nu_a;
            yaw = 0.5 * delta_t*delta_t * nu_yawdd;
            }
        else {
            ppx = px + v/yawd * (sin(yaw+yawd*delta_t) - sin(yaw)) \
                + 0.5 *delta_t*delta_t*cos(yaw)*nu_a;
            ppy = py + v/yawd * (-cos(yaw+yawd*delta_t) + cos(yaw)) \
                + 0.5 *delta_t*delta_t*sin(yaw)*nu_a;
            yaw = yawd*delta_t + 0.5*delta_t*delta_t*nu_yawdd;
            }

        pv = v + delta_t * nu_a;
        pyawd = yawd + delta_t * nu_yawdd;

        Xsig_pred_.col(i) << ppx, ppy, pv, pyaw, pyawd;
}

void UKF::PredictMeanAndCovariance() {
    //set weights
    weights_(0) = lambda_ / (lambda_ + n_aug_);

    for (int i = 1; i < n_sigma_; i++) {
        weights_(i) = 1/(2*(lambda_ + n_aug_));
    }
    //predict state mean
    x_ = Xsig_pred_.col(0) * weights_(0);
    for (int i = 1; i < n_sigma_; i++) {
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }
    //predict state covariance matrix
    VectorXd diff = VectorXd(n_x_);
    diff = (Xsig_pred_.col(0) - x_.col(0));
    P_ = weights_(0) * (diff * diff.transpose());
    for (int i = 1; i < n_sigma_; i++) {
        diff = (Xsig_pred_.col(i) - x_);
        P_ = P_ + weights_(i) * (diff * diff.transpose());
    }
}

void UKF::PredictRadarMeasurement() {
    double px, py, v, yaw, yawd;
    double rho, phi, rhod;

    for (int i = 0; i < n_sigma_; i++) {
        double px = Xsig_pred_(0,i);
        double py = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);

        rho = sqrt(px*px + py*py);
        phi = atan2(py,px);
        rhod = (px*cos(yaw)*v + py*sin(yaw)*v) / rho;

        Z_rad_sig_.col(i) << rho, phi, rhod;
    }

    Z_rad_sig_pred_ = weights_(0) * Z_rad_sig_.col(0);
    for (int i = 1; i < n_sigma_; i++) {
        Z_rad_sig_pred_ = Z_rad_sig_pred_ + weights_(i) * Z_rad_sig_.col(i);
    }

    S_rad_.fill(0);
    MatrixXd diff = MatrixXd(3, 3);
    for (int i = 0; i < n_sigma_; i++) {
        diff = Z_rad_sig_.col(i) - Z_rad_sig_pred_;
        if (diff(1) < -M_PI) {
            diff(1) += 2*M_PI;
        }
        else if (diff(1) > M_PI) {
            diff(1) -= 2*M_PI;
        }
        S_rad_ = S_rad_ + weights_(i) * (diff * diff.transpose());
        }

    S_rad_(0, 0) += std_radr_*std_radr_;
    S_rad_(1 ,1) += std_radphi_*std_radphi_;
    S_rad_(2, 2) += std_radrd_*std_radrd_;
}

void UKF::PredictLidarMeasurement() {
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        Z_las_sig_(0) = Xsig_pred_(0,i);
        Z_las_sig_(1) = Xsig_pred_(1,i);
    }

    Z_las_sig_pred = weights_(0) * Z_las_sig_.col(0);
    for (int i = 1; i < n_sigma_; i++) {
        Z_las_sig_pred_ = Z_las_sig_pred_ + weights_(i) * Z_las_sig_.col(i);
    }

    S_las_.fill(0);
    MatrixXd diff = MatrixXd(3, 3);
    for (int i = 0; i < n_sigma_; i++) {
        diff = Z_las_sig_.col(i) - Z_las_sig_pred_;
        S_las_ = S_las_ + weights_(i) * (diff * diff.transpose());
        }

    S_las_(0, 0) += std_laspx_*std_laspx_;
    S_las_(1 ,1) += std_laspy_*std_laspy_;
}

