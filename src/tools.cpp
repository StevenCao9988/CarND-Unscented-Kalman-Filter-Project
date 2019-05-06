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
	// check the validity of input  输入数据的合法性检测，增加程序鲁棒性
	// The estimation vector size should not be zero 测量数据的长度不能为0
	if (estimations.size() == 0) {
		cout << "Inout is empty!" << endl;
		return rmse;
	}
	// The estimation vector size should equal ground truth vector size 测量数据长度应该与真实值的长度一致
	if (estimations.size() != ground_truth.size()) {
		cout << "estimations and ground_truth should be the same size!" << endl;
		return rmse;
	}
	// Accumulate squared residuals 累积误差的平方
	for (unsigned int i = 0; i < estimations.size(); i++) {
		VectorXd residual = estimations[i] - ground_truth[i];
		residual = residual.array() * residual.array();
		rmse += residual;
	}

	// calculate the mean 计算平均值
	rmse = rmse / estimations.size();
	rmse = rmse.array().sqrt();
	return rmse;
}