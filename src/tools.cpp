#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{
	// Initialize Variables
		VectorXd rmse(4);
		rmse << 0, 0, 0, 0;

	// Check the validity of inputs
		if( ( estimations.size() != ground_truth.size() ) || (estimations.size() == 0) )
		{
			cout << "Invalid estimation or grount_truth data" << endl;
			return rmse;
		} 

	// Accumulate squared residuals
		for(unsigned int i=0; i < estimations.size(); ++i)
		{
			VectorXd residual = estimations[i] - ground_truth[i];

			// Coefficient-wise multiplication
				residual = residual.array()*residual.array();
				rmse += residual;
		}

	// Calculate the Mean
		rmse = rmse/estimations.size();

	// Square Root
		rmse = rmse.array().sqrt();

	// // Print Sum
	// 	cout << "Sum = " << rmse.sum() << endl;

	// Return
		return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

	MatrixXd Hj(3,4);

	// Recover state parameters
		float px = x_state(0);
		float py = x_state(1);
		float vx = x_state(2);
		float vy = x_state(3);

	// Pre-compute a set of terms to avoid repeated calculation
		float c1 = px*px + py*py;
		float c2 = sqrt(c1);
		float c3 = (c1 * c2);

	// Check devision by zero
		if(fabs(c1) < 0.001)
		{
			cout << "\n-------------------" << endl;
			cout << "		Error" << endl;
			cout << "---------------------" << endl;
			cout << "CalculateJacobian() -- divide by zero" << endl;
			return Hj;
		}

	// Compute Jacobian
		Hj << (px/c2), (py/c2), 0, 0,
			  -(py/c1), (px/c1), 0, 0,
			  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}
