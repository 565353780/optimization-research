#include "Interior_Point_Solver.h"

int main()
{
	//Interior_Point_Solver solver(10, 6, 4, 1.1, 1e-10);

	Interior_Point_Solver solver(2, 1, 1, 1.1, 1e-10);

	//Interior_Point_Solver solver(500, 7, 3, 1.1, 1e-10);

	MatrixXd x;

	solver.solve(x, 1.0);

	cout << "Solve finished! Results are : " << endl << endl;

	cout << "x = " << solver.x_.transpose() << endl << endl;

	cout << "Gx-h = " << (solver.G * solver.x_ - solver.h).transpose() << endl << endl;

	cout << "Ax-b = " << (solver.A * solver.x_ - solver.b).transpose() << endl << endl;

	/*Interior_Point_Solver solver(2, 1, 1, 1.1, 1e-10);

	MatrixXd P;
	MatrixXd Q;
	MatrixXd r;
	MatrixXd G;
	MatrixXd h;
	MatrixXd A;
	MatrixXd b;

	MatrixXd x;

	solver.P = P;
	solver.Q = Q;
	solver.r = r;
	solver.G = G;
	solver.h = h;
	solver.A = A;
	solver.b = b;

	solver.x_ = x;

	solver.solve(x, 1.0);*/

	return 1;
}