#include "Interior_Point_Solver.h"

Interior_Point_Solver::Interior_Point_Solver()
{

}

Interior_Point_Solver::Interior_Point_Solver(int n, int m, int p, double gamma, double epsilon)
{
	n_ = n;
	m_ = m;
	p_ = p;
	gamma_ = gamma;
	epsilon_ = epsilon;

	setAllMatrix();

	frame = Mat(500, 500, CV_8UC3, Scalar(255, 255, 255));
	start_p = Point(frame.rows / 2, frame.cols / 2);
	bigger_times = 10;
}

Interior_Point_Solver::~Interior_Point_Solver()
{

}

void Interior_Point_Solver::set_Matrix(MatrixXd& matrix, int row, int col)
{
	matrix.resize(row, col);
	matrix.setZero();

	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
		{
			matrix(i, j) = 1.0 * (rand() % 10 - 5);
		}
	}
}

void Interior_Point_Solver::setAllMatrix()
{
	set_Matrix(P, n_, n_);
	P = P * P.transpose();
	if (n_ == 2)
	{
		P.setIdentity();
	}

	set_Matrix(Q, n_, 1);

	set_Matrix(r, 1, 1);

	set_Matrix(G, m_, n_);

	set_Matrix(h, m_, 1);

	set_Matrix(A, p_, n_);
	while (true)
	{
		FullPivHouseholderQR<MatrixXd> qr(A);
		if (qr.rank() == p_)
		{
			break;
		}
		set_Matrix(A, p_, n_);
	}

	set_Matrix(b, p_, 1);

	v_.resize(p_, 1);

	v_.setZero();

	if (n_ < 100)
	{
		cout << P << endl << endl;
		cout << Q << endl << endl;
		cout << r << endl << endl;
		cout << G << endl << endl;
		cout << h << endl << endl;
		cout << A << endl << endl;
		cout << b << endl << endl;
	}
}

void Interior_Point_Solver::set_param(MatrixXd& x, double t)
{
	x_ = x;
	t_ = t;

	set_Matrix(x_, n_, 1);

	bool init_ok = false;

	while (!init_ok)
	{
		init_ok = true;

		if (fai(x_) == Inf)
		{
			init_ok = false;

			set_Matrix(x_, n_, 1);
		}
		else
		{
			MatrixXd temp_check = G * x_ - h;

			for (int i = 0; i < m_; ++i)
			{
				if (temp_check(i, 0) > 0)
				{
					init_ok = false;

					set_Matrix(x_, n_, 1);

					break;
				}
			}
		}
	}

	cout << "init x and t done." << endl;
}

double Interior_Point_Solver::f_0(MatrixXd& x)
{
	return (1.0 / 2.0 * (x.transpose() * P * x) + Q.transpose() * x + r)(0, 0);
}

double Interior_Point_Solver::f_i(MatrixXd& x, int i)
{
	return (G.row(i) * x - h.row(i))(0, 0);
}

double Interior_Point_Solver::fai(MatrixXd& x)
{
	double sum = 0;

	for (int i = 0; i < m_; ++i)
	{
		double current_f_i = f_i(x, i);

		if (current_f_i < 0)
		{
			sum -= log(-current_f_i);
		}
		else
		{
			return Inf;
		}
	}

	return sum;
}

double Interior_Point_Solver::f(MatrixXd& x, double t)
{
	double fai_ = fai(x);

	if (fai_ != Inf)
	{
		return t * f_0(x) + fai_;
	}

	return Inf;
}

MatrixXd Interior_Point_Solver::d_f_0(MatrixXd& x)
{
	return P * x + Q;
}

MatrixXd Interior_Point_Solver::d_f_i(MatrixXd& x, int i)
{
	return G.row(i).transpose();
}

MatrixXd Interior_Point_Solver::d_fai(MatrixXd& x)
{
	MatrixXd d_sum(n_, 1);

	d_sum.setZero();

	for (int i = 0; i < m_; ++i)
	{
		double current_f_i = f_i(x, i);

		if (current_f_i < 0)
		{
			d_sum -= 1.0 / current_f_i * d_f_i(x, i);
		}
		//else if (current_f_i > 0)
		else
		{
			for (int j = 0; j < n_; ++j)
			{
				d_sum(j, 0) = Inf;
			}

			return d_sum;
		}
	}

	return d_sum;
}

MatrixXd Interior_Point_Solver::d_f(MatrixXd& x, double t)
{
	MatrixXd d_fai_ = d_fai(x);

	if (d_fai_(0, 0) != Inf)
	{
		return t * d_f_0(x) + d_fai_;
	}

	return d_fai_;
}

MatrixXd Interior_Point_Solver::d2_f_0(MatrixXd& x)
{
	return P.transpose();
}

MatrixXd Interior_Point_Solver::d2_f_i(MatrixXd& x, int i)
{
	MatrixXd d2(n_, n_);
	d2.setZero();

	return d2;
}

MatrixXd Interior_Point_Solver::d2_fai(MatrixXd& x)
{
	MatrixXd d2_sum(n_, n_);

	d2_sum.setZero();

	for (int i = 0; i < m_; ++i)
	{
		double current_f_i = f_i(x, i);

		if (current_f_i < 0)
		{
			d2_sum += 1.0 / current_f_i / current_f_i * d_f_i(x, i) * d_f_i(x, i).transpose();

			d2_sum -= 1.0 / current_f_i * d2_f_i(x, i);
		}
		//else if (current_f_i > 0)
		else
		{
			for (int j = 0; j < n_; ++j)
			{
				for (int k = 0; k < n_; ++k)
				{
					d2_sum(j, k) = Inf;
				}
			}

			return d2_sum;
		}
	}

	return d2_sum;
}

MatrixXd Interior_Point_Solver::d2_f(MatrixXd& x, double t)
{
	MatrixXd d2_fai_ = d2_fai(x);

	if (d2_fai_(0, 0) != Inf)
	{
		return t * d2_f_0(x) + d2_fai_;
	}

	return d2_fai_;
}

bool Interior_Point_Solver::check_feasible(MatrixXd& x, double t)
{
	double f_ = f(x, t);

	if (f_ == Inf)
	{
		return false;
	}

	MatrixXd d_f_ = d_f(x, t);

	if (d_f_(0, 0) == Inf)
	{
		return false;
	}

	MatrixXd d2_f_ = d2_f(x, t);

	if (d2_f_(0, 0) == Inf)
	{
		return false;
	}

	return true;
}

Point2d Interior_Point_Solver::map_point(double x, double y)
{
	Point2d p(x, y);

	return bigger_times * p + start_p;
}

MatrixXd Interior_Point_Solver::compute_x_star()
{
	clock_t clk_start = clock();

	if (n_ == 2)
	{
		x_list.clear();

		x_list.push_back(map_point(x_(0, 0), x_(1, 0)));

		cout << x_.transpose() << endl;
	}

	MatrixXd Left(n_ + p_, n_ + p_);

	Left.setZero();

	Left.block(0, n_, n_, p_) = A.transpose();

	Left.block(n_, 0, p_, n_) = A;

	double tao = 0.49;

	double gamma = 0.5;

	MatrixXd old_x = x_;

	while (true)
	{
		Left.block(0, 0, n_, n_) = d2_f(x_, t_);

		MatrixXd Right(n_ + p_, 1);

		Right.setZero();

		Right.block(0, 0, n_, 1) -= d_f(x_, t_);

		Right.block(0, 0, n_, 1) -= A.transpose() * v_;

		Right.block(n_, 0, p_, 1) -= A * x_ - b;

		MatrixXd result = Left.fullPivHouseholderQr().solve(Right);

		MatrixXd delta = result.block(0, 0, n_, 1);
		//delta = -d_f(x_, t_);

		MatrixXd delta_v = result.block(n_, 0, p_, 1);

		double alpha = 1.0;

		MatrixXd r_dual = d_f(x_, t_);

		r_dual += A.transpose() * v_;

		MatrixXd r_pri = A * x_ - b;

		double norm = sqrt(r_dual.norm() * r_dual.norm() + r_pri.norm() * r_pri.norm());
		//norm = r_dual.norm();

		while (true)
		{
			MatrixXd current_x = x_ + alpha * delta;

			MatrixXd current_v = v_ + alpha * delta_v;

			MatrixXd current_r_dual = d_f(current_x, t_);

			current_r_dual += A.transpose() * current_v;

			MatrixXd current_r_pri = A * current_x - b;
			
			if (!check_feasible(current_x, t_))
			{
				alpha *= gamma;

				continue;
			}
			
			double current_norm = sqrt(current_r_dual.norm() * current_r_dual.norm() + current_r_pri.norm() * current_r_pri.norm());
			//current_norm = current_r_dual.norm();

			double current_coeff = 1.0 - tao * alpha;

			if (current_norm <= current_coeff * norm || alpha == 0)
			{
				break;
			}

			alpha *= gamma;
		}

		if (clock() - clk_start > 5000)
		{
			t_ = m_ / epsilon_ + 1;

			return x_;
		}

		x_ += alpha * delta;

		v_ += alpha * delta_v;

		if (n_ == 2)
		{
			x_list.push_back(map_point(x_(0, 0), x_(1, 0)));

			line(frame, x_list[x_list.size() - 2], x_list[x_list.size() - 1], Scalar(0, 255, 0), 1);

			flip(frame, frame_show, 0);

			imshow("test", frame_show);

			waitKey(1);
		}

		/*if (clock() - clk_start > 10000)
		{
			t_ = m_ / epsilon_;

			break;
		}*/

		if (norm <= epsilon_)
		{
		cout << "=========================================" << endl;
			cout << "feasible : " << check_feasible(x_, t_) << endl;
			cout << "f_0 = " << f_0(x_) << endl;
			cout << "fai = " << fai(x_) << endl;
			cout << "f = " << f(x_, t_) << endl;
			cout << "d_f = " << d_f(x_, t_).transpose() << endl;
			cout << "Ax -b = " << (A * x_ - b).transpose() << endl;
			cout << "norm = " << norm << endl;
			cout << "alpha = " << alpha << endl;
			cout << "delta = " << delta.transpose() << endl;
			cout << "x = " << x_.transpose() << endl;
			cout << "t = " << t_ << endl;
			cout << "=========================================" << endl;

			break;
		}
	}

	return x_;
}

bool Interior_Point_Solver::is_finished()
{
	if (m_ / t_ < epsilon_)
	{
		return true;
	}

	return false;
}

void Interior_Point_Solver::update_t()
{
	t_ *= gamma_;
}

void Interior_Point_Solver::solve(MatrixXd& x, double t)
{
	set_param(x, t);

	if (n_ == 2)
	{
		draw_f_value();
	}

	while (!is_finished())
	{
		compute_x_star();

		update_t();
	}

	waitKey(0);

	x = x_;
}

void Interior_Point_Solver::draw_f_value()
{
	line(frame, Point2d(0, frame.cols / 2), Point2d(frame.rows, frame.cols / 2), Scalar(0, 0, 0), 1);
	line(frame, Point2d(frame.rows / 2, 0), Point2d(frame.rows / 2, frame.cols), Scalar(0, 0, 0), 1);

	line(frame, map_point(-4, 1), map_point(-4, 1), Scalar(0, 0, 255), 1);

	double cut_num = 100;

	MatrixXd f_value(100, 100);

	f_value.setZero();

	double value_min = 0;
	double value_max = 0;

	for (int i = -50; i < 50; ++i)
	{
		for (int j = -50; j < 50; ++j)
		{
			MatrixXd cur_x(n_, 1);
			cur_x(0, 0) = i;
			cur_x(1, 0) = j;

			if (f(cur_x, t_) == Inf)
			{
				f_value(i + 50, j + 50) = Inf;
				continue;
			}

			f_value(i + 50, j + 50) = f(cur_x, t_) / t_;

			if (value_min == 0 && value_max == 0)
			{
				value_min = f_value(i + 50, j + 50);
				value_max = f_value(i + 50, j + 50);
			}
			else if (f_value(i + 50, j + 50) < value_min)
			{
				value_min = f_value(i + 50, j + 50);
			}
			else if (f_value(i + 50, j + 50) > value_max)
			{
				value_max = f_value(i + 50, j + 50);
			}
		}
	}

	value_max -= value_min;

	if (cut_num > value_max)
	{
		cut_num = value_max;
	}

	for (int i = 0; i < 100; ++i)
	{
		for (int j = 0; j < 100; ++j)
		{
			if (f_value(i, j) > cut_num && f_value(i, j) != Inf)
			{
				f_value(i, j) = 255;
			}
			if (f_value(i, j) != Inf)
			{
				f_value(i, j) = 254 * (f_value(i, j) - value_min) / cut_num;
			}
		}
	}

	for (int i = -50; i < 50; ++i)
	{
		for (int j = -50; j < 50; ++j)
		{
			if (f_value(i + 50, j + 50) != Inf)
			{
				int cur_num = max(0, int(f_value(i + 50, j + 50)));
				cur_num = min(255, cur_num);

				line(frame, map_point(i, j), map_point(i, j), Scalar(255 - cur_num, 0, cur_num), 3);
			}
		}
	}
}

void Interior_Point_Solver::draw_norm_value()
{
	line(frame, Point2d(0, frame.cols / 2), Point2d(frame.rows, frame.cols / 2), Scalar(0, 0, 0), 1);
	line(frame, Point2d(frame.rows / 2, 0), Point2d(frame.rows / 2, frame.cols), Scalar(0, 0, 0), 1);

	line(frame, map_point(-4, 1), map_point(-4, 1), Scalar(0, 0, 255), 1);

	double cut_num = 10;

	MatrixXd f_value(100, 100);

	f_value.setZero();

	double value_min = 0;
	double value_max = 0;

	for (int i = -50; i < 50; ++i)
	{
		for (int j = -50; j < 50; ++j)
		{
			MatrixXd cur_x(n_, 1);
			cur_x(0, 0) = i;
			cur_x(1, 0) = j;

			MatrixXd r_dual = d_f(cur_x, t_);

			r_dual += A.transpose() * v_;

			MatrixXd r_pri = A * cur_x - b;

			double norm = sqrt(r_dual.norm() * r_dual.norm() + r_pri.norm() * r_pri.norm());

			if (f(cur_x, t_) == Inf)
			{
				f_value(i + 50, j + 50) = Inf;
				continue;
			}

			f_value(i + 50, j + 50) = norm / t_;

			if (value_min == 0 && value_max == 0)
			{
				value_min = f_value(i + 50, j + 50);
				value_max = f_value(i + 50, j + 50);
			}
			else if (f_value(i + 50, j + 50) < value_min)
			{
				value_min = f_value(i + 50, j + 50);
			}
			else if (f_value(i + 50, j + 50) > value_max)
			{
				value_max = f_value(i + 50, j + 50);
			}
		}
	}

	value_max -= value_min;

	if (cut_num > value_max)
	{
		cut_num = value_max;
	}

	for (int i = 0; i < 100; ++i)
	{
		for (int j = 0; j < 100; ++j)
		{
			if (f_value(i, j) > cut_num && f_value(i, j) != Inf)
			{
				f_value(i, j) = 255;
			}
			if (f_value(i, j) != Inf)
			{
				f_value(i, j) = 254 * (f_value(i, j) - value_min) / cut_num;
			}
		}
	}

	for (int i = -50; i < 50; ++i)
	{
		for (int j = -50; j < 50; ++j)
		{
			if (f_value(i + 50, j + 50) != Inf)
			{
				int cur_num = max(0, int(f_value(i + 50, j + 50)));
				cur_num = min(255, cur_num);

				line(frame, map_point(i, j), map_point(i, j), Scalar(255 - cur_num, 0, cur_num), 3);
			}
		}
	}
}