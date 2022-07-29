#include "Sparse_Optimization_Solver.h"

Sparse_Optimization_Solver::Sparse_Optimization_Solver()
{

}

Sparse_Optimization_Solver::Sparse_Optimization_Solver(int n, int m, int T, double mu, double rou, double gamma)
{
	n_ = n;
	m_ = m;
	T_ = T;
	mu_ = mu;
	rou_ = rou;
	gamma_ = gamma;

	setAllMatrix();
}

Sparse_Optimization_Solver::~Sparse_Optimization_Solver()
{


}

void Sparse_Optimization_Solver::set_Matrix(MatrixXd& matrix, int row, int col)
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

void Sparse_Optimization_Solver::set_Matrix(MatrixXcd& matrix, int row, int col)
{
	matrix.resize(row, col);
	matrix.setZero();

	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
		{
			matrix(i, j) = std::complex<double>(1.0 * (rand() % 10 - 5), 1.0 * (rand() % 10 - 5));
		}
	}
}

void Sparse_Optimization_Solver::setAllMatrix()
{
	set_Matrix(A, n_, m_);

	set_Matrix(B, n_, 1);

	y_.resize(n_, 1);
	y_.setZero();

	z_.resize(m_, 1);
	z_.setZero();

	u_y_.resize(n_, 1);
	u_y_.setZero();

	u_z_.resize(m_, 1);
	u_z_.setZero();

	cout << A << endl << endl;
	cout << B << endl << endl;
}

void Sparse_Optimization_Solver::set_param(MatrixXcd& x)
{
	set_Matrix(x, m_, 1);

	x_ = x;
}

double Sparse_Optimization_Solver::norm_21(MatrixXcd & matrix)
{
	double norm_ = 0;

	for (int i = 0; i < matrix.rows(); ++i)
	{
		double current_norm = 0;

		for (int j = 0; j < matrix.cols(); ++j)
		{
			current_norm += matrix(i, j).real() * matrix(i, j).real();
			current_norm += matrix(i, j).imag() * matrix(i, j).imag();
		}

		norm_ += sqrt(current_norm);
	}

	return norm_;
}

double Sparse_Optimization_Solver::norm_11(MatrixXcd & matrix)
{
	double norm_ = 0;

	for (int i = 0; i < matrix.rows(); ++i)
	{
		for (int j = 0; j < matrix.cols(); ++j)
		{
			double current_norm = 0;

			current_norm += matrix(i, j).real() * matrix(i, j).real();
			current_norm += matrix(i, j).imag() * matrix(i, j).imag();

			norm_ += sqrt(current_norm);
		}
	}

	return norm_;
}

void Sparse_Optimization_Solver::update_x()
{
	MatrixXcd Left(m_, m_);
	Left.setZero();

	MatrixXcd Right(m_, 1);
	Right.setZero();

	for (int i = 0; i < m_; ++i)
	{
		Left(i, i) += 1.0;

		Right(i, 0) += z_(i, 0) - u_z_(i, 0);

		for (int j = 0; j < n_; ++j)
		{
			for (int k = 0; k < m_; ++k)
			{
				Left(i, k) += A(j, i) * A(j, k);
			}

			Right(i, 0) += A(j, i) * (y_(j, 0) + B(j, 0) - u_y_(j, 0));
		}
	}

	x_ = Left.fullPivHouseholderQr().solve(Right);
}

void Sparse_Optimization_Solver::update_y()
{
	MatrixXcd Left(n_, n_);
	Left.setZero();

	MatrixXcd Right(n_, 1);
	Right.setZero();

	for (int i = 0; i < n_; ++i)
	{
		Left(i, i) += 2.0 + rou_;
	}

	Right += rou_ * (A * x_ - B + u_y_);

	y_ = Left.fullPivHouseholderQr().solve(Right);
}

void Sparse_Optimization_Solver::update_z()
{
	MatrixXcd Left(m_, m_);
	Left.setZero();

	MatrixXcd Right(m_, 1);
	Right.setZero();

	for (int i = 0; i < m_; ++i)
	{
		Left(i, i) += 2.0 * mu_ + rou_;
	}

	Right += rou_ * (x_ + u_z_);

	z_ = Left.fullPivHouseholderQr().solve(Right);
}

void Sparse_Optimization_Solver::update_u()
{
	u_y_ += gamma_ * (A * x_ - y_ - B);

	u_z_ += gamma_ * (x_ - z_);
}

void Sparse_Optimization_Solver::solve(MatrixXcd& x)
{
	MatrixXcd AX_B = A * x_ - B;

	double current_norm = mu_ * norm_21(x_) + norm_11(AX_B);

	double old_norm = Inf;

	cout << "====================================================" << endl;
	cout << "x = " << x_.transpose() << endl;
	cout << "y = " << y_.transpose() << endl;
	cout << "z = " << z_.transpose() << endl;
	cout << "u y = " << u_y_.transpose() << endl;
	cout << "u z = " << u_z_.transpose() << endl;
	cout << "Ax-b = " << (A * x_ - B).transpose() << endl;
	cout << "norm = " << current_norm << endl;
	cout << "====================================================" << endl;

	while(true)
	{
		old_norm = current_norm;

		update_x();

		update_y();

		update_z();

		update_u();

		AX_B = A * x_ - B;

		current_norm = mu_ * norm_21(x_) + norm_11(AX_B);

		if ((old_norm - current_norm) / current_norm < 1e-10)
		{
			break;
		}

		cout << "====================================================" << endl;
		cout << "x = " << x_.transpose() << endl;
		cout << "y = " << y_.transpose() << endl;
		cout << "z = " << z_.transpose() << endl;
		cout << "u y = " << u_y_.transpose() << endl;
		cout << "u z = " << u_z_.transpose() << endl;
		cout << "Ax-b = " << (A * x_ - B).transpose() << endl;
		cout << "norm = " << current_norm << endl;
		cout << "====================================================" << endl;
	}
}