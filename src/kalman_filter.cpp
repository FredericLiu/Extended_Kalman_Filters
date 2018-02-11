#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ =  F_*x_;
  MatrixXd Ft = F_.transpose();
  P_ = (F_*P_*Ft)+ Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  
  
	VectorXd y_;
	y_ = z - (H_*x_);
	
	MatrixXd Ht = H_.transpose();
	MatrixXd S_ = (H_*P_*Ht) + R_;
	MatrixXd Si = S_.inverse();
	
	MatrixXd K_;
	K_ = P_*Ht*Si;
	// new state
	MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
	x_ = x_+ (K_*y_);
	P_ = (I-K_*H_)*P_;
  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  
  
    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);
    
	float rho = sqrt(px * px + py * py);
	if (fabs(rho)<0.0001){
		rho = 0.0001;
	}
	
    float phi = atan2(py, px);
    float rhoDot = (px * vx + py * vy) / rho;
    VectorXd h = VectorXd(3);
    h << rho, phi, rhoDot;
	
  	VectorXd y_= z - h;
	
	// adjust the phi angle so that it is between -pi and pi
    if (y_(1) > M_PI) {
		y_(1) = y_(1) - M_PI * 2.0;
    } else if (y_(1) < -M_PI) {
		y_(1) += M_PI * 2.0;
    }
	
	MatrixXd Hj = tools_.CalculateJacobian(x_);
	MatrixXd Hjt = Hj.transpose();
	MatrixXd S_ = (Hj*P_*Hjt) + R_;
	MatrixXd Si = S_.inverse();
	MatrixXd PHt = P_ * Hjt;
	MatrixXd K_ = PHt*Si;
	
	// new state
	x_ = x_+ (K_*y_);
	MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
	P_ = (I-K_*Hj)*P_;
  
}
