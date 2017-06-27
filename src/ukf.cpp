#include <iostream>
#include "ukf.h"
using namespace std;
/**
* Initializes Unscented Kalman filter
*/
UKF::UKF() {
  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  //lambda_ = 3 - n_aug_;
  lambda_ = 3 - n_aug_;
  NIS_radar_ = 1;
  NIS_laser_ = 1;
  previous_timestamp_ = 0;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // augmented sigma point state vector
  x_aug_ = VectorXd(n_aug_);

  // augmented sigma points  covariance matrix
  P_aug_ = MatrixXd(n_aug_, n_aug_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

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

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //create matrix with predicted augmented sigma points as columns
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // weights defined
  weights_ = VectorXd(2 * n_aug_ + 1);

  //create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z_ = 3;

  //create matrix for sigma points in measurement space for radar
  Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

  z_pred_ = VectorXd(n_z_);

  Tc_ = MatrixXd(n_x_, n_z_);

  //measurement covariance matrix S
  S_ = MatrixXd(n_z_, n_z_);

  z_radar_ = VectorXd(n_z_);

  // lidar incoming measurement vector
  z_lidar_ = VectorXd(2);

  // noise for laser sensors px and py
  R_laser_ = MatrixXd(2, 2);  // laser sensor measurement covariance matrix
  R_laser_ << std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
}



UKF::~UKF() {}
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /////////////  Initialization ////////////////////////
  // First measurement
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      if (((fabs(meas_package.raw_measurements_[0])) < 0.0001) && ((fabs(meas_package.raw_measurements_[1])) < 0.0001))
      {
        meas_package.raw_measurements_[0] = 0.001;
        meas_package.raw_measurements_[1] = 0.001;
      }
      //Converting radar measurements: from polar to cartesian
      double x_cart = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
      double y_cart = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
      //Initialize the state x_ with the first measurement
      x_ << x_cart, y_cart, 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      if (((fabs(meas_package.raw_measurements_[0])) < 0.0001) && ((fabs(meas_package.raw_measurements_[1])) < 0.0001))
      {
        meas_package.raw_measurements_[0] = 0.001;
        meas_package.raw_measurements_[1] = 0.001;
      }
      /////// Init states ////////
      double x = meas_package.raw_measurements_[0];
      double y = meas_package.raw_measurements_[1];
      //Initialize the state x_ with the first measurement
      x_ << x, y, 0, 0, 0;
    }
    // The first measurement is reliable so use the
    // identity matrix as the initial covariance matrix

    P_ << MatrixXd::Identity(n_x_, n_x_);

    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = meas_package.timestamp_;
    return;
  }


  ////////////////  Prediction  ///////////////////////////
  // 3 steps:
  // (1) AugmentedSigmaPoints
  // (2) SigmaPointPrediction
  // (3) PredictMeanAndCovariance
  //////////////////////////////////////////////////////////

  delta_t_ = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  AugmentedSigmaPoints();
  SigmaPointPrediction(delta_t_);
  PredictMeanAndCovariance();

  //////////  Measurement Update  ////////////////////////
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
  else {
    PredictRadarMeasurement();
    UpdateRadar(meas_package);
  }
  ///////////////////////////////////////////////////////

}


void UKF::AugmentedSigmaPoints() {

  //create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5, 5) = P_;
  P_aug_(5, 5) = std_a_*std_a_;
  P_aug_(6, 6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();
  //create augmented sigma points
  Xsig_aug_.col(0) = x_aug_;
  for (int i = 0; i< n_aug_; i++)
  {

    Xsig_aug_.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
  }
}

void UKF::SigmaPointPrediction(const double delta_t) {

  //predict sigma points
  for (int i = 0; i< 2 * n_aug_ + 1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0, i);
    double p_y = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double yaw = Xsig_aug_(3, i);
    double yawd = Xsig_aug_(4, i);
    double nu_a = Xsig_aug_(5, i);
    double nu_yawdd = Xsig_aug_(6, i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v / yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
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
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
}


void UKF::PredictMeanAndCovariance() {
  //Next x_ and P_ are calculated here
  // with the weight added

  //  Step 3 in Lecture Diagram
  // set weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i<2 * n_aug_ + 1; i++) {  //2n+1 weights
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }

  //predicted state mean -
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  //cout << "predicted state mean" << endl;
  //cout << x_ << endl;
  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //	VectorXd x_diff = Xsig_pred_.col(i) -Xsig_pred_.col(0);

    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}


void UKF::PredictRadarMeasurement() {
  //transform sigma points into measurement space
  //calculate mean predicted measurement
  //calculate measurement covariance matrix S

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    //////////////////////////////////////////////
    if (fabs(p_x)<0.001 && fabs(p_y)<0.001) {
      p_x = 0.001;
      p_y = 0.001;
    }
    //////////////////////////////////////////////

    Zsig_(0, i) = sqrt(p_x*p_x + p_y*p_y);  // r
    Zsig_(1, i) = atan2(p_y, p_x);			// phi
    Zsig_(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);  // r_dot


  }

  //mean predicted measurement

  z_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    //cout << "z_diff = Zsig_.col(i) - z_pred_" << z_diff << endl;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_, n_z_);
  R << std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0, std_radrd_*std_radrd_;
  S_ = S_ + R;
}


/**
* Updates the state and the state covariance matrix using a laser measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  // z is the measured value of the sensors
  z_lidar_ <<
           meas_package.raw_measurements_[0],
          meas_package.raw_measurements_[1];

  H_ = MatrixXd(2, 5);
  H_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0;
  VectorXd z_pred = H_ * x_;
  VectorXd y = z_lidar_ - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);

  long long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

  //You'll also need to calculate the lidar NIS.
  NIS_laser_ = (z_pred - z_lidar_).transpose()*
               Si*(z_pred - z_lidar_);
  cout << "NIS_laser: " << NIS_laser_ << endl;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {


  double rho = sqrt((x_[0] * x_[0] + x_[1] * x_[1]));
  double phi = atan2(x_[1], x_[0]);
  double rho_dot = (x_[0] * x_[2] + x_[1] * x_[3]) / sqrt((x_[0] * x_[0] + x_[1] * x_[1]));

  //z_pred_<< rho, phi, rho_dot;

  z_radar_ <<
           meas_package.raw_measurements_[0],
          meas_package.raw_measurements_[1],
          meas_package.raw_measurements_[2];

  //calculate cross correlation matrix
  Tc_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

    Tc_ = Tc_ + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc_ * S_.inverse();

  //residual
  VectorXd z_diff = z_radar_ - z_pred_;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_*K.transpose();
  // calculate the radar NIS.
  //

  // update NIS
  NIS_radar_ = (z_pred_ - z_radar_).transpose()*
               S_.inverse()*(z_pred_ - z_radar_);
  cout << "NIS_radar: " << NIS_radar_ << endl;
}

