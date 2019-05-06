#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;

#define EPS 0.0001 // A very small number

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 6;       

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 7;  
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */
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
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;
  // time when the state is true, in us
  time_us_ = 0.0;
   //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //define spreading parameter
  lambda_ = 3 - n_x_;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //create vector for weights
  weights_ = VectorXd(2 * n_aug_ + 1);

  // the current NIS for radar
  NIS_radar_ = 0.0;

  // the current NIS for laser
  NIS_laser_ = 0.0;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) 
{
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   * total 3 steps  �ܹ�3������
   */
	// ǰ�����ж�������������ȷ
	if (((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_) ||
		((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_))
	{
		// 1. initialization ��ʼ��
		//*************************************************************************
		if (is_initialized_ != true)
		{
			// Э������� ��ʼ��
			P_ << 
				0.15, 0, 0, 0, 0,
				0, 0.15, 0, 0, 0,
				0, 0, 1, 0, 0,
				0, 0, 0, 1, 0,
				0, 0, 0, 0, 1;

			// ʱ��� ��ʼ�� 
			time_us_ = meas_package.timestamp_;

			// �������� ��ʼ��
			x_ << 1, 1, 1, 1, 0.1;
			if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_)
			{
				float x = meas_package.raw_measurements_(0);
				float y = meas_package.raw_measurements_(1);;
				x_(0) = x;
				x_(1) = y;
			}
			else if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_)
			{
				float ro = meas_package.raw_measurements_(0);
				float theta = meas_package.raw_measurements_(1);
				float ro_dot = meas_package.raw_measurements_(2);
				x_(0) = ro * cos(theta);
				x_(1) = ro * sin(theta);
			}

			// ��ʼ�����
			is_initialized_ = true;
		}
		else
		{
			// 2. Prediction Ԥ��
			//*************************************************************************
			// Calculate the timestep between measurements in seconds �������ι۲�֮���ʱ���
			float dt = (meas_package.timestamp_ - time_us_);
			time_us_ = meas_package.timestamp_;
			cout << "dt" << dt << endl;
			if ((dt > 0) && (dt < 10000000))
			{
				dt /= 1000000.0;  // ΢��ת������
			}
			else 
			{
				cout << "error time data!" << endl;
				dt = 0;
				is_initialized_ = false;
			}// if (dt > 0) ����

			Prediction(dt);

			// 3. Update ����
			//*************************************************************************
			if (dt > 0)
			{
				if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_)
				{
					UpdateLidar(meas_package);
				}
				else if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_)
				{
					UpdateRadar(meas_package);
				}// ���½���
			}
			
		}// if (is_initialized_ != true) ����
	}// �ж�������������ȷ ����

}

void UKF::Prediction(double delta_t) 
{
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   * total 3 steps  �ܹ�4������
   */

	// 1. Augmentation Assignment
	//*************************************************************************
	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	//create augmented mean state
	x_aug.head(n_x_) = x_;
	x_aug(n_aug_ - 2) = 0;
	x_aug(n_aug_ - 1) = 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(n_aug_ - 2, n_aug_ - 2) = std_a_ * std_a_;
	P_aug(n_aug_ - 1, n_aug_ - 1) = std_yawdd_ * std_yawdd_;

	//create square root matrix
	MatrixXd Root_aug = P_aug.llt().matrixL();   // �����ƽ����

	//create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i < n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * Root_aug.col(i);    // ��Ҫ��дΪn_x, 
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * Root_aug.col(i);
	}// 1.end Augmentation Assignment

	// 2. Sigma Point Prediction Assignment
	//*************************************************************************

	//predict sigma points
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		double p_x = Xsig_aug(0, i);
		double p_y = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);          // ���ٶȵ�����
		double nu_yawd = Xsig_aug(6, i);       // �Ǽ��ٶȵ�����

		// Ԥ�⿪ʼ  
		double EXS = 0.001;					   // һ���ǳ�С����
		double p_x_p = p_x;					   // p_xԤ��     
		double p_y_p = p_y;					   // p_yԤ�� 
		double v_p = v;						   // vԤ�� 
		double yaw_p = yaw + yawd * delta_t;   // yawԤ��
		double yawd_p = 0;					   // yawdԤ��

		//avoid division by zero
		if (fabs(yawd) > EXS)
		{
			p_x_p = p_x_p + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
			p_y_p = p_y_p + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t) );
		}
		else
		{
			p_x_p = p_x_p + v * delta_t * cos(yaw);
			p_y_p = p_y_p + v * delta_t * sin(yaw);
		}

		// �������
		p_x_p = p_x_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
		p_y_p = p_y_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
		v_p = v_p + nu_a * delta_t;
		yaw_p = yaw_p + 0.5 * nu_yawd * delta_t * delta_t;
		yawd_p = yawd_p + delta_t * nu_yawd;

		//write predicted sigma points into right column  д��
		Xsig_pred_(0, i) = p_x_p;
		Xsig_pred_(1, i) = p_y_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}// 2. end Sigma Point Prediction Assignment

	// 3. Predicted Mean and Covariance Assignment
	//*************************************************************************

	//set weights
	weights_(0) = lambda_ / (lambda_ + n_aug_);
	for (int i = 1; i < 2 * n_aug_ + 1; i++)
	{
		weights_(i) = 0.5 / (lambda_ + n_aug_);
	}

	//predict state mean  ����Ԥ��ƽ��ֵ
	x_.fill(0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}

	//predict state covariance matrix
	P_.fill(0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		// ����Ԥ��ƽ��ֵ �� sigmaԤ��� ֮��Ĳ��
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		// �Ƕȹ淶��
		while (x_diff(3) > M_PI)
		{
			x_diff(3) = x_diff(3) - 2.0 * M_PI;
		}

		while (x_diff(3) < -M_PI)
		{
			x_diff(3) = x_diff(3) + 2.0 * M_PI;
		}
		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
	}// end 3. Predicted Meanand Covariance Assignment
}

void UKF::UpdateLidar(MeasurementPackage meas_package) 
{
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

   // 1. ׼������
   //*************************************************************************
	// �������� 
	VectorXd z = meas_package.raw_measurements_;

	// lidar ���Բ���p_x��p_y������������ά��Ϊ2
	int n_z = 2;

	// �˾����sigma��Ͳ���ֵ��ϵ����
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	// ��sigama������������
	for (int i = 0; i < 2 * n_aug_ + 1; i++)  // 2 * n_aug_ + 1 ��sigma��
	{
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);

		Zsig(0, i) = p_x;
		Zsig(1, i) = p_y;
	}
    
	// �Ƶ� ������ƽ��ֵ
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	// �Ƶ� ������Э�������
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);

	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S = S + weights_(i) * z_diff * z_diff.transpose();
	}	 

	// ��Ӳ�������
	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_laspx_ * std_laspx_, 0,
		0, std_laspy_* std_laspy_;
	S = S + R;

	// ���¾��� Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	// 2. ��ʼ����
	//*************************************************************************
	Tc.fill(0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	// kalman gain K;
	MatrixXd K = Tc * S.inverse();

	// 
	VectorXd z_diff = z - z_pred;

	// ���� NIS
	NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

	// ����ƽ��ֵ Э�������
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S * K.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) 
{
  /**`
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

	// 1. ׼������
	//*************************************************************************
	
	// �������� 
	VectorXd z = meas_package.raw_measurements_;

	//set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 3;
	
	// �˾����sigma��Ͳ���ֵ��ϵ����
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	// ��sigama������������
	for (int i = 0; i < 2 * n_aug_ + 1; i++)  // 2 * n_aug_ + 1 ��sigma��
	{
		// ����ת��
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);
     
		double v1 = v * cos(yaw);
		double v2 = v * sin(yaw);

		Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);   // r  �뾶
		Zsig(1, i) = atan2(p_y, p_x);   // phi �Ƕ�
		Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x*p_x + p_y*p_y);   // r_dot �����ٶ�
	} 

	// �Ƶ� ������ƽ��ֵ
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	// �Ƶ� ������Э�������
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);

	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;
		// �Ƕȹ淶��
		while (z_diff(1) > M_PI)
		{
			z_diff(1) = z_diff(1) - 2.0 * M_PI;
		}
		while (z_diff(1) < -M_PI)
		{
			z_diff(1) = z_diff(1) + 2.0 * M_PI;
		}

		S = S + weights_(i) * z_diff * z_diff.transpose();		
	}

	// ��Ӳ�������
	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_radr_* std_radr_, 0, 0,
		0, std_radphi_* std_radphi_, 0,
		0, 0, std_radrd_* std_radrd_;
	S = S + R;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	// 1. ��ʼ����
	//*************************************************************************

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd z_diff = Zsig.col(i) - z_pred;
		// �Ƕȹ淶��
		while (z_diff(1) > M_PI)
		{
			z_diff(1) = z_diff(1) - 2.0 * M_PI;
		}
		while (z_diff(1) < -M_PI)
		{
			z_diff(1) = z_diff(1) + 2.0 * M_PI;
		}

		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		// �Ƕȹ淶��
		while (x_diff(3) > M_PI)
		{
			x_diff(3) = x_diff(3) - 2.0 * M_PI;
		}
		while (x_diff(3) < -M_PI)
		{
			x_diff(3) = x_diff(3) + 2.0 * M_PI;
		}
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//update state mean and covariance matrix
	VectorXd z_diff = z - z_pred;   // z����������Ԥ��Ĳ��
	  // �Ƕȹ淶��
	while (z_diff(1) > M_PI)
	{
		z_diff(1) = z_diff(1) - 2.0 * M_PI;
	}
	while (z_diff(1) < -M_PI)
	{
		z_diff(1) = z_diff(1) + 2.0 * M_PI;
	}

	//calculate NIS
	NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

	x_ = x_ + K * z_diff;
	P_ = P_ - K * S * K.transpose();
}