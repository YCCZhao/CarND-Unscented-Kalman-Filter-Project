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

  // initial state vector
  x = VectorXd(5);
  x.fill(1);

  // initial covariance matrix
  P = MatrixXd::Identity(5, 5)


  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

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

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // Weights of sigma points
  weights_ = VectorXd(n_sigma_);

  // matrix of predicted state sigma points
  Xsig_ = MatrixXd(n_aug_, n_sigma_);

  // mean of predicted state sigma points
  Xsig_pred_ = VectorXd(n_x_);

  // mean and a state sigma point difference
  Xdiff_ = MatrixXd(n_x_, n_aug_);


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
    if (!is_initialized_) {
    /**
      * Initialize the state x with the first measurement.
      */

        InitiateMeasurement()

        previous_timestamp_ = measurement_pack.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

  /*****************************************************************************
   *  Turn on/off Radar/Laser Measurement
   ****************************************************************************/
    if ((measurement_pack.sensor_type_ == MeasurementPackage::RADAR & !use_radar_) || \
    (measurement_pack.sensor_type_ == MeasurementPackage::LASER & !use_laser_)) {
        return;
    }

  /*****************************************************************************
   *  Predict
   ****************************************************************************/

    // time step
    float delta_t_ = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;

    // predict k+1
    Prediction(delta_t_);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

    Update();

}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t_) {
  /**
  Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    //
    AugmentedSigmaPoints();

    //
    SigmaPointPrediction(delta_t_);

    //
    PredictMeanAndCovariance();

}
/**
 * Updates the state and the state covariance matrix using a sensor measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::Update(MeasurementPackage meas_package) {
  //
  PredictMeasurement(meas_package);
  UpdateState(meas_package);
}

void UKF::InitiateMeasurement(MeasurementPackage meas_package) {
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
}

void UKF::AugmentedSigmaPoints() {

    VectorXd Xaug_ = VectorXd(n_aug_);
    MatrixXd Paug_ = MatrixXd(n_aug_, n_aug_);

    Xaug_.fill(0);
    Xaug_.head(n_x_) = x;

    Paug_.fill(0.0);
    Paug_.topLeftCorner(n_x_, n_x_) = P;
    Paug_(5, 5) = std_a_;
    Paug_(6, 6) = std_yawdd_;

    MatrixXd Psqrt_ = P_.llt().matrixL();

    Xsig_.col(0) = Xaug_;
    for (int i = 0; i < n_aug_; i++)
    {
        Xsig_.col(i + 1)          = Xaug_ + sqrt(lambda_ + n_aug_) * Psqrt_.col(i);
        Xsig_.col(i + 1 + n_aug_) = Xaug_ - sqrt(lambda_ + n_aug_) * Psqrt_.col(i);
    }
}

void UKF::SigmaPointPrediction(delta_t_) {

    double px_, py_, v_, yaw_, yawd_;
    double nu_a_, nu_yawdd_;

    double ppx_, ppy_, pv_, pyaw_, pyawd_;

    for (int i = 0; i < n_sigma_; i++) {

        px_   = Xsig_(0, i);
        py_   = Xsig_(1, i);
        v_    = Xsig_(2, i);
        yaw_  = Xsig_(3, i);
        yawd_ = Xsig_(4, i);

        nu_a_    = Xsig_(5, i);
        nu_yawdd_ = Xsig_(6, i);

        if (yawd_ < 0.0001) {
            ppx_ = px_ + v_*cos(yaw_)*delta_t_ + 0.5*delta_t_*delta_t_*cos(yaw_)*nu_a_;
            ppy_ = py_ + v_*sin(yaw_)*delta_t_ + 0.5*delta_t_*delta_t_*sin(yaw_)*nu_a_;
            yaw_ = 0.5 * delta_t_*delta_t_ * nu_yawdd_;
            }
        else {
            ppx_ = px_ + v_/yawd_ * (sin(yaw_+yawd_*delta_t_) - sin(yaw_)) \
                 + 0.5 * delta_t_*delta_t_*cos(yaw_)*nu_a_;
            ppy_ = py_ + v_/yawd_ * (-cos(yaw_+yawd_*delta_t_) + cos(yaw_)) \
                 + 0.5 * delta_t_*delta_t_*sin(yaw_)*nu_a_;
            yaw_ = yawd_ * delta_t_ + 0.5 * delta_t_*delta_t_ * nu_yawdd_;
            }

        pv_    = v_ + delta_t_ * nu_a_;
        pyawd_ = yawd_ + delta_t_ * nu_yawdd_;

        Xsig_pred_.col(i) << ppx_, ppy_, pv_, pyaw_, pyawd_;
}
}

void UKF::PredictMeanAndCovariance() {

    VectorXd diff = VectorXd(n_x_);

    // set weights
    weights_(0) = lambda_ / (lambda_ + n_aug_);

    for (int i = 1; i < n_sigma_; i++) {
        weights_(i) = 1 / (2 * (lambda_ + n_aug_));
    }

    // predict state mean
    x = Xsig_pred_.col(0) * weights_(0);
    for (int i = 1; i < n_sigma_; i++) {
        x = x + weights_(i) * Xsig_pred_.col(i);
    }

    // mean and sigma point different
    for (int i = 0; i< n_sigma_; i++) {
        Xdiff_.col(i) = Xsig_pred_.col(i) - x;
    }

    // predict covariance matrix
    P = weights_(0) * (diff * diff.transpose());
    for (int i = 1; i < n_sigma_; i++) {
        P = P + weights_(i) * (Xdiff_.col(i) * Xdiff_.col(i).transpose());
    }

}

void UKF::PredictMeasurement(MeasurementPackage meas_package) {

    int n_z_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

        double px_, py_, v_, yaw_, yawd_;
        double rho_, phi_, rhod_;

        n_z_ = 3;

        for (int i = 0; i < n_sigma_; i++) {

            px_  = Xsig_pred_(0, i);
            py_  = Xsig_pred_(1, i);
            v_   = Xsig_pred_(2, i);
            yaw_ = Xsig_pred_(3, i);

            rho_  = sqrt(px*px + py*py);
            phi_  = atan2(py, px);

            if (rho < 0.0001) {
                rhod_ = 0.0001;
            }
            else {
            rhod_ = (px*cos(yaw)*v + py*sin(yaw)*v) / rho;
            }

            Zsig_.col(i) << rho, phi, rhod;
        }
    }

    else {

        n_z_ = 2;

        for (int i = 0; i < n_sigma_; i++) {
        Z_las_sig_(0) = Xsig_pred_(0,i);
        Z_las_sig_(1) = Xsig_pred_(1,i);
        }
    }

    Zsig_      = MatrixXd(n_z_, n_sigma_);
    Zsig_pred_ = VectorXd(n_z_);
    Zdiff_     = MatrixXd(n_z_, n_sigma_);
    S_         = MatrixXd(n_z_, n_z_);

    Zsig_pred_ = weights_(0) * Zsig_.col(0);
    for (int i = 1; i < n_sigma_; i++) {
        Zsig_pred_ = Zsig_pred_ + weights_(i) * Zsig_.col(i);
    }

    for (int i = 0; i < n_sigma_; i++) {
        Zdiff_.col(i) = Zsig_.col(i) - Zsig_pred_;
    }

    S_.fill(0);
    for (int i = 0; i < n_sigma_; i++) {
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            if (Zdiff_(1, i) < -M_PI) {
                Zdiff_(1, i) += 2*M_PI;
            }
            else if (Zdiff_(1, i) > M_PI) {
                Zdiff_(1, i) -= 2*M_PI;
            }
        }

        S_ = S_ + weights_(i) * (Zdiff_.col(i) * Zdiff_.col(i).transpose());
    }

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

        S_(0, 0) += std_radr_*std_radr_;
        S_(1 ,1) += std_radphi_*std_radphi_;
        S_(2, 2) += std_radrd_*std_radrd_;
    }
    else {
        S_(0, 0) += std_laspx_*std_laspx_;
        S_(1 ,1) += std_laspy_*std_laspy_;
    }
}

void UKF::UpdateState(MeasurementPackage meas_package) {
  /**

  Complete this function! Use sensor data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  calculate the lidar NIS.
  */

    MatrixXd Tc_ = MatrixXd(n_x_, 2);
    MatrixXd K_  = MatrixXd(n_x, 2);
    VectorXd z_  = meas_package.raw_measurements_;


    //calculate cross correlation matrix
    Tc.fill(0);

    for (int i = 0; i < n_sigma_; i++) {
        Tc_ = Tc_ + weights_(i) * (Xdiff_.col(i) * Zdiff_.col(i).transpose());
    }

    //calculate Kalman gain K;
    K_ = Tc_ * S_.inverse();

    //update state mean and covariance matrix
    x += K_ * (z_ - Zsig_pred_);
    P -= K_ * S_ * K_.transpose();
}


