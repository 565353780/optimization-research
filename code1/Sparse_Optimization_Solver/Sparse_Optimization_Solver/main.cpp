#include "Sparse_Optimization_Solver.h"

int main()
{
	//Sparse_Optimization_Solver solver(10, 8, 1, 1.0, 0.1, 0.1);

	//Sparse_Optimization_Solver solver(50, 1000, 1, 1.0, 0.1, 0.1);

	Sparse_Optimization_Solver solver(2, 2, 1, 1.0, 0.1, 0.1);

	MatrixXcd x;

	solver.set_param(x);

	solver.solve(x);

	cout << "Solve finished! Results are :" << endl;

	cout << "x = " << solver.x_.transpose() << endl << endl;
	cout << "Ax-b = " << (solver.A * solver.x_ - solver.B).transpose() << endl << endl;

	/*Sparse_Optimization_Solver solver(2, 2, 1, 1.0, 0.1, 0.1);

	MatrixXcd A;

	MatrixXcd B;

	MatrixXcd x;

	solver.A = A;

	solver.B = B;

	solver.x_ = x;

	solver.solve(x);*/

	return 1;
}