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
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check if vectors are the same size
    if((estimations.size() == 0) || (estimations.size() != ground_truth.size())){
        cout << "Error! Estimations vector incorrect size" << endl;
        cout << "estimations.size(): " << estimations.size() << endl;
        cout << "ground_truth.size(): " << ground_truth.size() << endl;
        return rmse;
    }


    // accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
        VectorXd residual = estimations[i] - ground_truth[i];
//        cout << "estimations: " << estimations[i] << endl;
        // coeff-wise mult
        residual = residual.array() * residual.array();
        rmse += residual;
    }

    // calculate the mean
    rmse = rmse / estimations.size();

    // calc sq root
    rmse = rmse.array().sqrt();

    cout << "RMSE: " << rmse << endl;
    return rmse;
}