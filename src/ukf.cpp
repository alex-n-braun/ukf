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
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector; dim = 5
  n_x_ = 5;
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);
    // initialize to identity matrix;
  //HINT: to be optimized

  //define spreading parameter
  lambda_ = 3 - n_x_;

  // dimension of augmented state vector
  n_aug_ = n_x_ + 2;

  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;
    // according to lecture, this is a sensible value for cars;
    // however, we are dealing with bicycles here.
  //HINT: to be optimized.

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI;
    // starting value; needs to be adapted. Probably too high
  //HINT: to be optimized.

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
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage & measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    Initialize(measurement_pack);
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // delta time (in seconds)
  double dt((measurement_pack.timestamp_ - time_us_)*1.0e-6);

  Prediction(dt);

  // update time stamp
  time_us_ = measurement_pack.timestamp_;


  /*****************************************************************************
   *  Measurement Update
   ****************************************************************************/
  // check if radar or lidar measurement
  if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    UpdateLidar(measurement_pack);
  else
    UpdateRadar(measurement_pack);

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
/*    // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
//  auto px(x_(0)); not needed!
//  auto py(x_(1));
  auto v(x_(2));
  auto psi(x_(3));
  auto psi_d(x_(4));

  // compute delta for state vector
  double d_v(0.0);
  double d_psi(psi*delta_t);
  double d_psi_d(0.0);

  // computation of new position; distinguish psi == 0 and psi != 0
  double d_px;
  double d_py;
  if (abs(psi)<1e-15) { // psi is (close to) zero
    d_px = v * cos(psi) * delta_t;
    d_py = v * sin(psi) * delta_t;
  } else {
    double vpsi(v/psi);
    double psi_psi_d_t(psi + psi_d * delta_t);
    d_px = vpsi * (sin(psi_psi_d_t)-sin(psi));
    d_py = vpsi * (cos(psi)-cos(psi_psi_d_t));
  }

  VectorXd delta_x(5);
  delta_x << d_px, d_py, d_v, d_psi, d_psi_d;
  x_ = x_ + delta_x; */

/*  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  GenerateSigmaPoints(Xsig);*/

  {
    //create augmented sigma point matrix
    MatrixXd Xsig_aug(n_aug_, 2 * n_aug_ + 1);
    AugmentedSigmaPoints(Xsig_aug);

    // predict sigma points
    SigmaPointPrediction(Xsig_aug, delta_t);
  }

  // predict mean and covariance matrix
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage & measurement_pack) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  const int n_z(2);
  //create matrix for sigma points in measurement space
  MatrixXd Z_sig(n_z, 2 * n_aug_ + 1);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // predict measurement from state vector / sigma matrix
  PredictLaserMeasurement(Z_sig, z_pred, S);

  // update state and covariance, store NIS
  UpdateLaserMeasurement(Z_sig, z_pred, S, measurement_pack.raw_measurements_);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage & measurement_pack) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS == Normalized Innovation Squared.
  */

  const int n_z(3);
  //create matrix for sigma points in measurement space
  MatrixXd Z_sig(n_z, 2 * n_aug_ + 1);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // predict measurement from state vector / sigma matrix
  PredictRadarMeasurement(Z_sig, z_pred, S);

  // update state and covariance, store NIS
  UpdateRadarMeasurement(Z_sig, z_pred, S, measurement_pack.raw_measurements_);
}


void UKF::Initialize(const MeasurementPackage & measurement_pack) {
  cout << "UKF: " << endl;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    /**
    Convert radar from polar to cartesian coordinates and initialize state.
    */
    float rho(measurement_pack.raw_measurements_[0]);
    float phi(measurement_pack.raw_measurements_[1]);
    float c_phi=cos(phi); float s_phi=sin(phi);
    float rho_dot(measurement_pack.raw_measurements_[2]);
    x_<<rho*c_phi, rho*s_phi, rho_dot, 0.0, 0.0;
  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    /**
    Initialize state.
    */
    x_<<measurement_pack.raw_measurements_[0],
        measurement_pack.raw_measurements_[1], 0.0, 0.0, 0.0;
  }
  time_us_ = measurement_pack.timestamp_;

  // done initializing, no need to predict or update
  is_initialized_ = true;
}

/*void UKF::GenerateSigmaPoints(MatrixXd & Xsig) {
  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //calculate sigma points ...
  //set sigma points as columns of matrix Xsig

  double factor(sqrt(lambda_+n_x_));
  A *= factor;

  Xsig.col(0)=x_;
  for (int i(0); i<n_x_; ++i) {
    Xsig.col(i+1)=x_+A.col(i);
    Xsig.col(i+1+n_x_)=x_-A.col(i);
  }
}*/


void UKF::AugmentedSigmaPoints(MatrixXd & Xsig_aug) const {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  double factor(sqrt(lambda_+n_aug_));

  L *= factor;

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - L.col(i);
  }
}

void UKF::SigmaPointPrediction(const MatrixXd & Xsig_aug, double delta_t) {
  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 1e-3) {
        px_p = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw) );
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw + yawd*delta_t) );
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


void UKF::PredictMeanAndCovariance() {
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //write result
  x_ = x;
  P_ = P;
}

void UKF::PredictRadarMeasurement(MatrixXd & Z_sig, VectorXd & z_pred, MatrixXd & S) const {

  //set measurement dimension, radar can measure r, phi, and r_dot
  const int n_z(3);

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
    Z_sig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Z_sig(1,i) = atan2(p_y,p_x);                                 //phi
    Z_sig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement

  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Z_sig.col(i);
  }

  //measurement covariance matrix S

  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Z_sig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;
}


void UKF::UpdateRadarMeasurement(const MatrixXd & Z_sig, const VectorXd & z_pred, const MatrixXd & S,
                      const VectorXd & z) {
  //set measurement dimension, radar can measure r, phi, and r_dot
  const int n_z(3);

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff(Z_sig.col(i) - z_pred);
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd Si(S.inverse());
  MatrixXd K = Tc * Si;

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix & write result
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // compute NIS
  NIS_radar_ = z_diff.transpose() * Si * z_diff;
}

void UKF::PredictLaserMeasurement(MatrixXd & Z_sig, VectorXd & z_pred, MatrixXd & S) const {

  //set measurement dimension, radar can measure r, phi, and r_dot
  const int n_z(3);

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
    Z_sig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Z_sig(1,i) = atan2(p_y,p_x);                                 //phi
    Z_sig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement

  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Z_sig.col(i);
  }

  //measurement covariance matrix S

  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Z_sig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;
}


void UKF::UpdateLaserMeasurement(const MatrixXd & Z_sig, const VectorXd & z_pred, const MatrixXd & S,
                      const VectorXd & z) {
  //set measurement dimension, radar can measure r, phi, and r_dot
  const int n_z(3);

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff(Z_sig.col(i) - z_pred);
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd Si(S.inverse());
  MatrixXd K = Tc * Si;

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix & write result
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // compute NIS
  NIS_radar_ = z_diff.transpose() * Si * z_diff;
}

