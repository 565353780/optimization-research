#pragma once
#include <iostream>
#include <ctime>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

#include <opencv2/core/core.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include <opencv2/opencv.hpp>

using namespace std;
using namespace Eigen;
using namespace cv;

#define Inf 1.0e308

class Interior_Point_Solver
{
public:
	Interior_Point_Solver();
	Interior_Point_Solver(int n, int m, int p, double gamma, double epsilon);
	~Interior_Point_Solver();

	void set_Matrix(MatrixXd& matrix, int row, int col);

	void setAllMatrix();

	void set_param(MatrixXd& x, double t);

	double f_0(MatrixXd& x);

	double f_i(MatrixXd& x, int i);

	double fai(MatrixXd& x);

	double f(MatrixXd& x, double t);

	MatrixXd d_f_0(MatrixXd& x);

	MatrixXd d_f_i(MatrixXd& x, int i);

	MatrixXd d_fai(MatrixXd& x);

	MatrixXd d_f(MatrixXd& x, double t);

	MatrixXd d2_f_0(MatrixXd& x);

	MatrixXd d2_f_i(MatrixXd& x, int i);

	MatrixXd d2_fai(MatrixXd& x);

	MatrixXd d2_f(MatrixXd& x, double t);

	bool check_feasible(MatrixXd& x, double t);

	Point2d map_point(double x, double y);

	MatrixXd compute_x_star();

	bool is_finished();

	void update_t();

	void solve(MatrixXd& x, double t);

	void draw_f_value();

	void draw_norm_value();

public:
	int n_;
	int m_;
	int p_;

	double gamma_;
	double epsilon_;

	MatrixXd P;
	MatrixXd Q;
	MatrixXd r;
	MatrixXd G;
	MatrixXd h;
	MatrixXd A;
	MatrixXd b;

	double t_;
	MatrixXd x_;
	MatrixXd v_;

	vector<Point2d> x_list;
	Mat frame;
	Mat frame_show;
	Point2d start_p;
	double bigger_times;

	clock_t clk_start;
};

