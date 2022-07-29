#pragma once
#include <iostream>
#include <ctime>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

using namespace std;
using namespace Eigen;

#define Inf 1e308

class Sparse_Optimization_Solver
{
public:
	Sparse_Optimization_Solver();
	Sparse_Optimization_Solver(int n, int m, int T, double mu, double rou, double gamma);
	~Sparse_Optimization_Solver();

	void set_Matrix(MatrixXd& matrix, int row, int col);
	void set_Matrix(MatrixXcd& matrix, int row, int col);

	void setAllMatrix();

	void set_param(MatrixXcd& x);

	double norm_21(MatrixXcd & matrix);
	double norm_11(MatrixXcd & matrix);

	void update_x();
	void update_y();
	void update_z();
	void update_u();

	void solve(MatrixXcd& x);

public:
	int n_;
	int m_;
	int T_;
	double mu_;
	double rou_;
	double gamma_;

	MatrixXcd z_;
	MatrixXcd y_;

	MatrixXcd u_y_;
	MatrixXcd u_z_;

	MatrixXcd x_;
	MatrixXcd A;
	MatrixXcd B;
};

