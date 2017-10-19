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
  is_initialized = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar = true;

  // initial time
  time_us_ = 0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.9;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.9;

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
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);

  // mean and a state sigma point difference
  Xdiff_ = MatrixXd(n_x_, n_sigma_);


  // initial state vector
  x = VectorXd(n_x_);
  x.fill(0.0);

  // initial covariance matrix
  P = MatrixXd::Identity(5, 5);
  //P(0, 0) = 0.1;
  //P(1, 1) = 0.1;
  //P(2, 2) = 0.1;
  //P(3, 3) = 0.1;
  //P(4, 4) = 0.1;



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
    if (!is_initialized) {
    /**
      * Initialize the state x with the first measurement.
      */
        InitiateMeasurement(meas_package);
        time_us_ = meas_package.timestamp_;

        // done initializing, no need to predict or update
        is_initialized = true;
        return;
    }

  /*****************************************************************************
   *  Turn on/off Radar/Laser Measurement
   ****************************************************************************/
    if ((meas_package.sensor_type_ == MeasurementPackage::RADAR & !use_radar) || \
    (meas_package.sensor_type_ == MeasurementPackage::LASER & !use_laser)) {
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

    // Predict k+1 state of generated augmented points
    SigmaPointPrediction(delta_t_);

    //Calculated mean and covariance of augmented points
    PredictMeanAndCovariance();


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

  // update state using k+1 sensor measurement
  UpdateState(meas_package);

}

/**
 * Initiates state vector using first sensor measurement
 * @param {MeasurementPackage} meas_package The first measurement data of
 * either radar or laser.
 */
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
			
	    // set weights
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < n_sigma_; i++) {
        weights_(i) = 1 / (2 * (lambda_ + n_aug_));
    }
}

/**
 * Creates Augmented state vector from state vector.
 * then generates sigma points.
 */
void UKF::AugmentedSigmaPoints() {

    VectorXd Xaug_ = VectorXd(n_aug_);
    MatrixXd Paug_ = MatrixXd(n_aug_, n_aug_);

    Xaug_.fill(0.0);
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

/**
 * Predicts k+1 states of augmented sigma points.
 */
void UKF::SigmaPointPrediction(double delta_t_) {

    for (int i = 0; i < n_sigma_; i++)
	{
		//extract values for better readability
		double p_x = Xsig_(0, i);
		double p_y = Xsig_(1, i);
		double v = Xsig_(2, i);
		double yaw = Xsig_(3, i);
		double yawd = Xsig_(4, i);
		double nu_a = Xsig_(5, i);
		double nu_yawdd = Xsig_(6, i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd*delta_t_) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_t_));
		}
		else {
			px_p = p_x + v*delta_t_*cos(yaw);
			py_p = p_y + v*delta_t_*sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd*delta_t_;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a*delta_t_*delta_t_ * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t_*delta_t_ * sin(yaw);
		v_p = v_p + nu_a*delta_t_;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t_*delta_t_;
		yawd_p = yawd_p + nu_yawdd*delta_t_;

		//write predicted sigma point into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
}
}

/**
 * calculates the mean and the covariance
 * for augmented sigma points at k+1.
 */
void UKF::PredictMeanAndCovariance() {



    // predict state mean
	x.fill(0.0);
	for (int i = 0; i < n_sigma_; i++) {  //iterate over sigma points
		x = x + weights_(i) * Xsig_pred_.col(i);
    }
    // mean and sigma point different
    for (int i = 0; i < n_sigma_; i++) {
        Xdiff_.col(i) = Xsig_pred_.col(i) - x;

		while (Xdiff_(3, i) < -M_PI) Xdiff_(3, i) += 2*M_PI;
        while (Xdiff_(3, i) > M_PI) Xdiff_(3, i) -= 2*M_PI;
		
		//while (Xdiff_(4, i) < -M_PI) Xdiff_(4, i) += 2*M_PI;
        //while (Xdiff_(4, i) > M_PI) Xdiff_(4, i) -= 2*M_PI;



    }

    // predict covariance matrix
	P.fill(0.0);
    for (int i = 0; i < n_sigma_; i++) {
        P = P + weights_(i) * (Xdiff_.col(i) * Xdiff_.col(i).transpose());
    }
}

/**
 * converts k+1 state vectors to measurement vector
 * calculates the mean and the covariance of the
 * sigma point measurement prediction.
 */
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
                rhod_ = 0.0001;
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

	Zsig_pred_.fill(0.0);
    for (int i = 0; i < n_sigma_; i++) {
        Zsig_pred_ = Zsig_pred_ + weights_(i) * Zsig_.col(i);
    }

    for (int i = 0; i < n_sigma_; i++) {
        Zdiff_.col(i) = Zsig_.col(i) - Zsig_pred_;
    }

    S_.fill(0.0);
    for (int i = 0; i < n_sigma_; i++) {
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			while (Zdiff_(1, i) < -M_PI) Zdiff_(1, i) += 2*M_PI;
			while (Zdiff_(1, i) > M_PI) Zdiff_(1, i) -= 2*M_PI;

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

/**
 * Update state vectors with the knowledge of predicted measurement and
 * the actual measurement difference.
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
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
    Tc_.fill(0.0);
    for (int i = 0; i < n_sigma_; i++) {
        Tc_ = Tc_ + weights_(i) * (Xdiff_.col(i) * Zdiff_.col(i).transpose());
    }

    //calculate Kalman gain K;
    K_ = Tc_ * S_.inverse();

    //update state mean and covariance matrix
	VectorXd z_cor_ = (z_ - Zsig_pred_);
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		
			while (z_cor_(1) < -M_PI) z_cor_(1) += 2*M_PI;
			while (z_cor_(1) > M_PI) z_cor_(1) -= 2*M_PI;

        }

	//calculate NIS
	
    double NIS = z_cor_.transpose() * S_.inverse() * z_cor_;
	cout<<"NIS:"<<endl;
	cout<<NIS<<endl;
	
    x += K_ * z_cor_;
    P -= K_ * S_ * K_.transpose();
	while (x(3)> M_PI) x(3)-=2.*M_PI;
    while (x(3)<-M_PI) x(3)+=2.*M_PI;
	
	cout<<x<<endl;
	cout<<P<<endl;
}


