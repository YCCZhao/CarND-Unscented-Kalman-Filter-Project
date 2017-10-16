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

  // initial time
  time_us_ = 0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(n_sigma_);

  // matrix of predicted state sigma points
  Xsig_ = MatrixXd(n_aug_, n_sigma_);

  // mean of predicted state sigma points
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);

  // mean and a state sigma point difference
  Xdiff_ = MatrixXd(n_x_, n_sigma_);


  // initial state vector
  x = VectorXd(n_x_);
  x.fill(1);
  x(2) = 5.2;
  //x(3) = 0.00005;
  x(4) = 0.01;

  // initial covariance matrix
  P = MatrixXd::Identity(5, 5);
  //P(0, 0) = 5;
  //P(1, 1) = 0.15;
  P(2, 2) = 0.5;
  P(4, 4) = 0.01;
  P(3, 3) = 0.01;

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
	  
        InitiateMeasurement(meas_package);
        time_us_ = meas_package.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

  /*****************************************************************************
   *  Turn on/off Radar/Laser Measurement
   ****************************************************************************/
    if ((meas_package.sensor_type_ == MeasurementPackage::RADAR & !use_radar_) || \
    (meas_package.sensor_type_ == MeasurementPackage::LASER & !use_laser_)) {
        return;
    }

  /*****************************************************************************
   *  Predict
   ****************************************************************************/

    // time step
    float delta_t_ = (meas_package.timestamp_ - time_us_)/1000000.0;

    // predict k+1
    Prediction(delta_t_);
  /*****************************************************************************
   *  Update
   ****************************************************************************/

    Update(meas_package);

    time_us_ = meas_package.timestamp_;
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t_) {
  /**
  Estimate the object's location. Modify the state
  vector, x. Predict sigma points, the state, and the state covariance matrix.
  */

    // Generate augmented points
    AugmentedSigmaPoints();
	//  cout << "augmentation" << endl;


    // Predict k+1 state of generated augmented points
    SigmaPointPrediction(delta_t_);
	//  cout << "sigma predict" << endl;

	
    //Calculated mean and covariance of augmented points
    PredictMeanAndCovariance();
	 // cout << "predict mean" << endl;


}
/**
 * Updates the state and the state covariance matrix using a sensor measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::Update(MeasurementPackage meas_package) {
  /**
  Combine predicted object's location measurement and actual sensor measurement,
  to update the estimated object's location.
  Modify the state vector and the state covariance matrix.
  */

  // Convert predicted k+1 location to sensor measurement
  PredictMeasurement(meas_package);
  //cout << "predictmeasurement" << endl;

  // update state using k+1 sensor measurement
  UpdateState(meas_package);
   // cout << "update"<< endl;

}

void UKF::InitiateMeasurement(MeasurementPackage meas_package) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

            float rho_     = meas_package.raw_measurements_(0);
            float phi_     = meas_package.raw_measurements_(1);
            float rho_dot_ = meas_package.raw_measurements_(2);
			
			float vx_ = rho_dot_ * cos(phi_);
			float vy_ = rho_dot_ * sin(phi_);
			float v_  = sqrt(vx_ * vx_ + vy_ * vy_);
			
            x(0) = rho_ * cos(phi_);
            x(1) = rho_ * sin(phi_);
			x(2) = v_;

            }

        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

            x(0) = meas_package.raw_measurements_(0);
            x(1) = meas_package.raw_measurements_(1);

            }
}

void UKF::AugmentedSigmaPoints() {

    VectorXd Xaug_ = VectorXd(n_aug_);
    MatrixXd Paug_ = MatrixXd(n_aug_, n_aug_);

    Xaug_.fill(0);
    Xaug_.head(n_x_) = x;

    Paug_.fill(0.0);
    Paug_.topLeftCorner(n_x_, n_x_) = P;
    Paug_(5, 5) = std_a_*std_a_;
    Paug_(6, 6) = std_yawdd_*std_yawdd_;

    MatrixXd Psqrt_ = Paug_.llt().matrixL();

    Xsig_.col(0) = Xaug_;
    for (int i = 1; i < n_aug_; i++)
    {
        Xsig_.col(i + 1)          = Xaug_ + sqrt(lambda_ + n_aug_) * Psqrt_.col(i);
        Xsig_.col(i + 1 + n_aug_) = Xaug_ - sqrt(lambda_ + n_aug_) * Psqrt_.col(i);
    }
}

void UKF::SigmaPointPrediction(double delta_t_) {

    double px_, py_, v_, yaw_, yawd_;
    double nu_a_, nu_yawdd_;

    double ppx_, ppy_, pv_, pyaw_, pyawd_;

    for (int i = 0; i < n_sigma_; i++) {

        px_   = Xsig_(0, i);
        py_   = Xsig_(1, i);
        v_    = Xsig_(2, i);
        yaw_  = Xsig_(3, i);
        yawd_ = Xsig_(4, i);

        nu_a_     = Xsig_(5, i);
        nu_yawdd_ = Xsig_(6, i);

        if (fabs(yawd_) < 0.00001) {
            ppx_  = px_ + v_*cos(yaw_)*delta_t_ + 0.5*delta_t_*delta_t_*cos(yaw_)*nu_a_;
            ppy_  = py_ + v_*sin(yaw_)*delta_t_ + 0.5*delta_t_*delta_t_*sin(yaw_)*nu_a_;
            }

        else {
            ppx_  = px_ + v_/yawd_ * ( sin(yaw_+yawd_*delta_t_) - sin(yaw_)) \
                  + 0.5 * delta_t_*delta_t_*cos(yaw_)*nu_a_;
            ppy_  = py_ + v_/yawd_ * (-cos(yaw_+yawd_*delta_t_) + cos(yaw_)) \
                  + 0.5 * delta_t_*delta_t_*sin(yaw_)*nu_a_;
            }

        pv_    = v_ + delta_t_ * nu_a_;
	    pyaw_ = yaw_ + yawd_ * delta_t_ + 0.5 * delta_t_*delta_t_ * nu_yawdd_;
        pyawd_ = yawd_ + delta_t_ * nu_yawdd_;
cout<<"here"<<endl;
cout<<pyaw_<<endl;
        Xsig_pred_.col(i) << ppx_, ppy_, pv_, pyaw_, pyawd_;
}
}

void UKF::PredictMeanAndCovariance() {

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
	//cout<<"predicted"<<endl;
	//cout<<Xsig_pred_<<endl;
    // mean and sigma point different
    for (int i = 0; i< n_sigma_; i++) {
        Xdiff_.col(i) = Xsig_pred_.col(i) - x;
		
		if (Xdiff_(3, i) < -M_PI) {
			Xdiff_(3, i) += 2*M_PI;
            }
        else if (Xdiff_(3, i) > M_PI) {
			Xdiff_(3, i) -= 2*M_PI;
            }
			
		if (Xdiff_(4, i) < -M_PI) {
			Xdiff_(4, i) += 2*M_PI;
            }
        else if (Xdiff_(4, i) > M_PI) {
			Xdiff_(4, i) -= 2*M_PI;
            }
			
    }

    // predict covariance matrix
    P = weights_(0) * (Xdiff_.col(0) * Xdiff_.col(0).transpose());
    for (int i = 1; i < n_sigma_; i++) {
        P = P + weights_(i) * (Xdiff_.col(i) * Xdiff_.col(i).transpose());
    }
}

void UKF::PredictMeasurement(MeasurementPackage meas_package) {

    int n_z_   = meas_package.raw_measurements_.size() ;
    Zsig_      = MatrixXd(n_z_, n_sigma_);
    Zsig_pred_ = VectorXd(n_z_);
    Zdiff_     = MatrixXd(n_z_, n_sigma_);
    S_         = MatrixXd(n_z_, n_z_);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) { //radar

        double px_, py_, v_, yaw_;
        double v1, v2, rho_, phi_, rhod_;

        for (int i = 0; i < n_sigma_; i++) {

            px_  = Xsig_pred_(0, i);
            py_  = Xsig_pred_(1, i);
            v_   = Xsig_pred_(2, i);
            yaw_ = Xsig_pred_(3, i);
			
			v1 = cos(yaw_) * v_;
			v2 = sin(yaw_) * v_;

            rho_  = sqrt(px_*px_ + py_*py_);
            phi_  = atan2(py_, px_);

            if (fabs(rho_) < 0.0001) {
                rhod_ = 0;
            }
            else {
            rhod_ = (px_*v1 + py_*v2) / rho_;
            }

            Zsig_.col(i) << rho_, phi_, rhod_;
        }
    }

    else {//laser

        for (int i = 0; i < n_sigma_; i++) {
			Zsig_(0, i) = Xsig_pred_(0, i);
			Zsig_(1, i) = Xsig_pred_(1, i);
        }
    }


    Zsig_pred_ = weights_(0) * Zsig_.col(0);
    for (int i = 1; i < n_sigma_; i++) {
        Zsig_pred_ = Zsig_pred_ + weights_(i) * Zsig_.col(i);
    }

    for (int i = 0; i < n_sigma_; i++) {
        Zdiff_.col(i) = Zsig_.col(i) - Zsig_pred_;
    }

    S_.fill(0);
    for (int i = 0; i < n_sigma_; i++) {
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            if (Zdiff_(1, i) < -M_PI) {
                Zdiff_(1, i) += 2*M_PI;
            }
            else if (Zdiff_(1, i) > M_PI) {
                Zdiff_(1, i) -= 2*M_PI;
            }
        }

        S_ = S_ + weights_(i) * (Zdiff_.col(i) * Zdiff_.col(i).transpose());
    }

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

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
  position. Modify the state vector, x, and covariance, P.

  Calculate the lidar NIS.
  */
    int n_z_   = meas_package.raw_measurements_.size() ;

    MatrixXd Tc_ = MatrixXd(n_x_, n_z_);
    MatrixXd K_  = MatrixXd(n_x_, n_z_);
    VectorXd z_  = meas_package.raw_measurements_;


    //calculate cross correlation matrix
    Tc_.fill(0);
    for (int i = 0; i < n_sigma_; i++) {
        Tc_ = Tc_ + weights_(i) * (Xdiff_.col(i) * Zdiff_.col(i).transpose());
    }

    //calculate Kalman gain K;
    K_ = Tc_ * S_.inverse();

    //update state mean and covariance matrix
	VectorXd z_cor = (z_ - Zsig_pred_);
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            if (z_cor(1) < -M_PI) {
                z_cor(1) += 2*M_PI;
            }
            else if (z_cor(1) > M_PI) {
                z_cor(1) -= 2*M_PI;
            }
        }
		
    x += K_ * z_cor;
    P -= K_ * S_ * K_.transpose();
}


