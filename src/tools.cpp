#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    VectorXd result(4);
    result<<0,0,0,0;
    {
        auto it_e(estimations.cbegin());
        auto it_g(ground_truth.cbegin());
        for (; it_e!=estimations.cend() & it_g!=ground_truth.cend(); ++it_e, ++it_g) {
            VectorXd diff(*it_e-*it_g);
            diff = diff.array()*diff.array();
            result = result + diff;
        }
    }
    result /= estimations.size();
    result = result.array().sqrt();
    return result;
}
