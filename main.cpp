#include "headers.h"
#include "functions.h"

int main()
{
	setlocale(LC_CTYPE, "RUS");
	cout.precision(8);

	vector <Point> points;
	ImportBenchmarks(points);

	vector <Meas> sections;
	ImportMeaseurement(sections);

	MatrixXd W(8, 1);
	GetResiduals(sections,points, W);

	MatrixXd B(8, 13);
	setMatrixB(B);

	MatrixXd P(13, 13);
	setMatrixP(P, sections);

	MatrixXd W_h(8, 1); 
	setMatrixW_h(W_h, sections);

	MatrixXd A(5, 13);
	setMatrixA(A);

	MatrixXd R = B * P.inverse() * B.transpose();

	MatrixXd K = -R.inverse() * W;

	MatrixXd Qy = P.inverse() - P.inverse() * B.transpose() * R.inverse() * B * P.inverse();	//матрица превышений	cout << "Matrix_excess Qy" << "\n";

	MatrixXd Qh = A * Qy * A.transpose();

	MatrixXd V = P.inverse() * B.transpose() * K;

	double u = setu(u, V, P);

	MatrixXd mh(13, 1);
	setMatrixMh(mh, Qy, u);

	MatrixXd mH(5, 1);
	setMatrixMH(mH, Qh, u);	

	getAllValues(u, V, mh, mH, P, Qy);
	return 0;
}