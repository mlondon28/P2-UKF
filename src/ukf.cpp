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
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI / 13;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3 - n_x_;

  n_sig_ = 2 * n_aug_ + 1;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);

  P_ <<   1,0,0,    0,0,
          0,1,0,    0,0,
          0,0,1000, 0,0,
          0,0,0,    1,0,
          0,0,0,    0,1;

  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  //create vector for weights
  weights_ = VectorXd(n_sig_);
  for(int i = 0; i < n_sig_; i++){
    if( i == 0){
      weights_(i) = lambda_/(lambda_ + n_aug_);
    } else {
      weights_(i) = 1 / (2 * (lambda_ + n_aug_));
    }
  }

}

UKF::~UKF() {}

VectorXd UKF::Polar2Cartesian(const VectorXd &polar){
  float rho = polar(0);
  float phi = polar(1);
  float rho_dot = polar(2);

  float px = rho * cos(phi);
  float py = rho * sin(phi);
  float v = rho_dot;

  VectorXd Cartesian(3);
  Cartesian << px, py, v;

  return Cartesian;
}

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

  long long timestamp = meas_package.timestamp_;

  if(!is_initialized_) {
    previous_timestamp_ = timestamp;

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      VectorXd cartesian = Polar2Cartesian(meas_package.raw_measurements_);

      float px = cartesian(0);
      float py = cartesian(1);
      float v  = cartesian(2);

      x_ << px, py, v, 0, 0;
      cout << "x_: " << x_ << endl;

    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      /**
      Initialize state.
      */
      x_ <<   meas_package.raw_measurements_(0),
              meas_package.raw_measurements_(1),
              0,
              0,
              0;
    }
    cout << "x_: " << x_ << endl;

    is_initialized_ = true;
    return;
  }
  double dt = (timestamp - previous_timestamp_) / 1000000.;
  previous_timestamp_ = timestamp;

  cout << "dt: " << dt << endl;

  Prediction(dt);

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
    UpdateRadar(meas_package);
      cout << "radar update" << endl;
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
    UpdateLidar(meas_package);
      cout << "laser update" << endl;
  }

}

MatrixXd UKF::GenerateSigmaPoints(){
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(5) = x_;

    ///*initialize x_aug. top 5 rows == x_
    cout << "x_aug in gen sigma points function: " << x_aug << endl;

  MatrixXd Q = MatrixXd(2,2);
  Q <<  pow(std_a_, 2), 0,
        0,              pow(std_yawdd_, 2);

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2,2) = Q;

  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

  MatrixXd L = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++){
    Xsig_aug.col(i + 1)         = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_)= x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  return Xsig_aug;
}


VectorXd UKF::ApplyMotionModel(VectorXd &x_aug, double dt){
  ///*extract values for better readability
  double px     = x_aug(0);
  double py     = x_aug(1);
  double v      = x_aug(2);
  double psi    = x_aug(3);
  double psid   = x_aug(4);
  double nu_a   = x_aug(5);
  double nu_psidd = x_aug(6);

  VectorXd state_trans_vec = VectorXd(n_x_);
  VectorXd covariance_trans_vec = VectorXd(n_x_);


  if(fabs(psid) < 0.00001){
      state_trans_vec <<  v * cos(psi) * dt,
                          v * sin(psi) * dt,
                          0,
                          psid * dt,
                          0;

  } else {
      state_trans_vec <<  (v / psid) * ( sin(psi + psid * dt) - sin(psi)),
                          (v / psid) * (-cos(psi + psid * dt) + cos(psi)),
                          0,
                          psid * dt,
                          0;
  }

  covariance_trans_vec <<     0.5 * dt * dt * cos(psi) * nu_a,
                              0.5 * dt * dt * sin(psi) * nu_a,
                              dt * nu_a,
                              0.5 * dt * dt * nu_psidd,
                              dt * nu_psidd;

  cout << "Xsig_aug + state + covariance models: " << x_aug.head(5) + state_trans_vec + covariance_trans_vec << endl;
  return x_aug.head(5) + state_trans_vec + covariance_trans_vec;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  MatrixXd XSig_aug = GenerateSigmaPoints();

  for(int i = 0; i < n_sig_; i++){
    VectorXd xsig = XSig_aug.col(i);
    Xsig_pred_.col(i) = ApplyMotionModel(xsig, dt);
  }

  x_.fill(0.0);
  for (int i = 0; i <n_sig_; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
//    x_diff(3) -= (2 * M_PI) * floor((x_diff(3) + M_PI) / (2 * M_PI));
    //angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

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

  /**
   * Use standard Kalman filter equations to update lidar
   */

  MatrixXd z = meas_package.raw_measurements_;

  int n_z = 2;

  MatrixXd H = MatrixXd(n_z, n_x_);
  MatrixXd R = MatrixXd(n_z, n_z);

  ///* init measurement matrix
  H <<  1,0,0,0,0,
        0,1,0,0,0;

  ///* init measurement covariance matrix
  R <<  pow(std_laspx_,2),  0,
        0,                  std_laspy_;

  VectorXd y = z - H * x_;
  MatrixXd S = H * P_ * H.transpose() + R;
  MatrixXd K = P_ * H.transpose() * S.inverse();

  ///* Update state (x_) and covariance (P_)
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;
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

  int n_z = 3;

  ///* radar sigma point matrix
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  ///* mean predicted measurement matrix
  VectorXd z_pred = VectorXd(n_z);

  ///* measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  z_pred.fill(0.0);

  ///* transform sigma points into measurement space
  for(int i = 0; i < n_sig_; i++){
    VectorXd xsig = Xsig_pred_.col(i);
    double px   = xsig(0);
    double py   = xsig(1);
    double v    = xsig(2);
    double psi  = xsig(3);
    double psid = xsig(4);

    ///* catch division by 0
    float div_zero_check = 0.00001;
    if (fabs(px) < div_zero_check){
      px = div_zero_check;
    }
    if (fabs(py) < div_zero_check) {
      py = div_zero_check;
    }

    ///* convert into symbolically correct representation
    float rho = sqrt(pow(px, 2) + pow(py, 2));
    float phi = atan2(py, px);
    float rhod = v * (px * cos(psi) + py * sin(psi)) / rho;

    ///* fill measurement vector column with measurements
    Zsig.col(i) << rho, phi, rhod;

    ///* calculate predicted mean measurement
    z_pred += weights_(i) * Zsig.col(i);
  }

  ///* Calculate predicted covariance matrix S
  S.fill(0.0);

  for(int i = 0; i < 2 *n_aug_ + 1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;

//    //angle normalization
//    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
//    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;


    S += weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z, n_z);

  R <<    pow(std_radr_,2),  0,                  0,
          0,                pow(std_radphi_,2),  0,
          0,                0,                  pow(std_radrd_,2);

  S += R;

  /**
   * Update with measured values
   */

  MatrixXd Tc = MatrixXd(n_x_,n_z);

  //calculate cross correlation matrix
  Tc.fill(0);
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points

    //residual
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;

////    //angle normalization
//    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
//    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
//    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
//    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    Tc +=  weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  MatrixXd z = meas_package.raw_measurements_;
  //update state mean and covariance matrix
  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K * S * K.transpose();

}