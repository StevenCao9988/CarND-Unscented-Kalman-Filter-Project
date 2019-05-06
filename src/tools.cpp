#include "tools.h"
#include <iostream>
using namespace std;

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;
	// check the validity of input  �������ݵĺϷ��Լ�⣬���ӳ���³����
	// The estimation vector size should not be zero �������ݵĳ��Ȳ���Ϊ0
	if (estimations.size() == 0) {
		cout << "Inout is empty!" << endl;
		return rmse;
	}
	// The estimation vector size should equal ground truth vector size �������ݳ���Ӧ������ʵֵ�ĳ���һ��
	if (estimations.size() != ground_truth.size()) {
		cout << "estimations and ground_truth should be the same size!" << endl;
		return rmse;
	}
	// Accumulate squared residuals �ۻ�����ƽ��
	for (unsigned int i = 0; i < estimations.size(); i++) {
		VectorXd residual = estimations[i] - ground_truth[i];
		residual = residual.array() * residual.array();
		rmse += residual;
	}

	// calculate the mean ����ƽ��ֵ
	rmse = rmse / estimations.size();
	rmse = rmse.array().sqrt();
	return rmse;
}