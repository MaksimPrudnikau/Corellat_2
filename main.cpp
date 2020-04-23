#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include "Eigen/Dense"
#include "boost/math/distributions.hpp"

using namespace std;
using namespace Eigen;
using namespace boost::math;
struct  Point
{
	string benchmark_name;
	double height;
};

void ImportBenchmarks(vector<Point>& points)
{
	ifstream heights_file("Var02_heights.txt");
	string benchmark_name;
	double height;
	vector <string> benchmark_names;
	vector<double> heights;
	while (heights_file >> benchmark_name >> height)
	{
		//cout << benchmark_name << "" << height << "\n";
		Point p;
		p.benchmark_name = benchmark_name;
		p.height = height;
		points.push_back(p);

	}
}
struct Meas
{
	string from_station;
	string to_station;
	double elevation;
	double distance;
	string leveling_class;
};
void ImportMeaseurement(vector<Meas>& sections)
{
	ifstream measurements_file("Var02_measurements.txt");
	string from_station;
	string to_station;
	double elevation;
	double distance;
	string leveling_class;


	while (measurements_file >> from_station >> to_station >> elevation
		>> distance >> leveling_class)
	{
		Meas m;
		m.from_station = from_station;
		m.to_station = to_station;
		m.elevation = elevation;
		m.distance = distance;
		m.leveling_class = leveling_class;
		sections.push_back(m);
	}
}

MatrixXd GetResiduals(vector <Meas>& sections,MatrixXd& W)
{
	W(0, 0) = -sections[0].elevation + sections[2].elevation - sections[9].elevation;
	W(1, 0) = sections[1].elevation + sections[11].elevation - sections[2].elevation;
	W(2, 0) = sections[3].elevation - sections[12].elevation - sections[1].elevation;
	W(3, 0) = -sections[6].elevation + sections[5].elevation + sections[12].elevation;
	W(4, 0) = sections[0].elevation + sections[10].elevation - sections[3].elevation;
	W(5, 0) = sections[8].elevation - sections[10].elevation - sections[7].elevation;
	W(6, 0) = sections[6].elevation - sections[8].elevation - sections[4].elevation;
	return W;
}

void setMatrixB(MatrixXd& B) {
	B(0, 0) = -1;	 B(1, 0) = 0;	 B(2, 0) = 0;    B(3, 0) = 0;	 B(4, 0) = 1;     B(5, 0) = 0;     B(6, 0) = 0;
	B(0, 1) = 0;	 B(1, 1) = 1;	 B(2, 1) = -1;   B(3, 1) = 0;	 B(4, 1) = 0;     B(5, 1) = 0;     B(6, 1) = 0; 
	B(0, 2) = 1;	 B(1, 2) = -1;	 B(2, 2) = 0;    B(3, 2) = 0;	 B(4, 2) = 0;     B(5, 2) = 0;     B(6, 2) = 0; 
	B(0, 3) = 0;	 B(1, 3) = 0;	 B(2, 3) = 1;    B(3, 3) = 0;	 B(4, 3) = -1;    B(5, 3) = 0;     B(6, 3) = 0; 
	B(0, 4) = 0;	 B(1, 4) = 0;	 B(2, 4) = 0;    B(3, 4) = 0;    B(4, 4) = 0;     B(5, 4) = 0;     B(6, 4) = -1; 
	B(0, 5) = 0;	 B(1, 5) = 0;	 B(2, 5) = 0;    B(3, 5) = 1;	 B(4, 5) = 0;     B(5, 5) = 0;     B(6, 5) = 0;
	B(0, 6) = 0;	 B(1, 6) = 0;	 B(2, 6) = 0;    B(3, 6) = -1;	 B(4, 6) = 0;     B(5, 6) = 0;     B(6, 6) = 1; 
	B(0, 7) = 0;	 B(1, 7) = 0;	 B(2, 7) = 0;    B(3, 7) = 0;	 B(4, 7) = 0;     B(5, 7) = -1;    B(6, 7) = 0; 
	B(0, 8) = 0;	 B(1, 8) = 0;	 B(2, 8) = 0;    B(3, 8) = 0;	 B(4, 8) = 0;     B(5, 8) = 1;     B(6, 8) = -1; 
	B(0, 9) = -1;	 B(1, 9) = 0;	 B(2, 9) = 0;    B(3, 9) = 0;	 B(4, 9) = 0;     B(5, 9) = 0;     B(6, 9) = 0; 
	B(0, 10) = 0;	 B(1, 10) = 0;	 B(2, 10) = 0;   B(3, 10) = 0;   B(4, 10) = 1;    B(5, 10) = -1;   B(6, 10) = 0;
	B(0, 11) = 0;	 B(1, 11) = 1;   B(2, 11) = 0;   B(3, 11) = 0;   B(4, 11) = 0;    B(5, 11) = 0;    B(6, 11) = 0;
	B(0, 12) = 0;	 B(1, 12) = 0;   B(2, 12) = -1;  B(3, 12) = 1;   B(4, 12) = 0;    B(5, 12) = 0;    B(6, 12) = 0;
}

void setMatrixP(MatrixXd& P, vector<Meas> sections) {
	for (size_t i = 0; i < P.rows(); ++i)
	{
		for (int j = 0; j < P.cols(); j++)
		{
			if (i == j)
			{
				P(i, j) = pow((1.0 / (10.0 * sqrt(sections[i].distance) / 1000.0)), 2);			}
			else 
			{
				P(i, j) = 0;
			}
		}
	}
		
}

//допустимые невязки
void setMatrixW_h(MatrixXd& W_h, vector<Meas> sections) {
	W_h(0, 0) = (10.0 / 1000.0) * sqrt(+sections[0].distance + sections[2].distance + sections[9].distance);
	W_h(1, 0) = (10.0 / 1000.0) * sqrt(sections[1].distance + sections[11].distance + sections[2].distance);
	W_h(2, 0) = (10.0 / 1000.0) * sqrt(sections[3].distance + sections[12].distance + sections[1].distance);
	W_h(3, 0) = (10.0 / 1000.0) * sqrt(+sections[6].distance + sections[5].distance + sections[12].distance);
	W_h(4, 0) = (10.0 / 1000.0) * sqrt(sections[0].distance + sections[10].distance + sections[3].distance);
	W_h(5, 0) = (10.0 / 1000.0) * sqrt(sections[8].distance + sections[10].distance + sections[7].distance);
	W_h(6, 0) = (10.0 / 1000.0) * sqrt(sections[6].distance + sections[8].distance + sections[4].distance);
}

void setMatrixA(MatrixXd& A) {
	//Rp 4 = 1 -------
	//Rp 5 = 2
	//Rp 6 = 3
	//Rp 7 = 4
	//Rp 3 = 5
	A(0, 0) = 1;  A(1, 0) = 0;   A(2, 0) = 0;	A(3, 0) = 0;  A(4, 0) = 0; 
	A(0, 1) = 0;  A(1, 1) = 1;   A(2, 1) = 0;	A(3, 1) = 0;  A(4, 1) = 1;
	A(0, 2) = 0;  A(1, 2) = 0;   A(2, 2) = 1;	A(3, 2) = 0;  A(4, 2) = 0;
	A(0, 3) = 1;  A(1, 3) = 0;   A(2, 3) = 0;   A(3, 3) = 1;  A(4, 3) = 1;
	A(0, 4) = 0;  A(1, 4) = 0;   A(2, 4) = 0;	A(3, 4) = 0;  A(4, 4) = 0;
	A(0, 5) = 0;  A(1, 5) = 0;   A(2, 5) = 0;	A(3, 5) = 0;  A(4, 5) = 0;
	A(0, 6) = 0;  A(1, 6) = 0;   A(2, 6) = 0;	A(3, 6) = 0;  A(4, 6) = 0;
	A(0, 7) = 0;  A(1, 7) = 0;   A(2, 7) = 0;	A(3, 7) = 0;  A(4, 7) = 0;
	A(0, 8) = 0;  A(1, 8) = 0;   A(2, 8) = 0;	A(3, 8) = 0;  A(4, 8) = -1;
	A(0, 9) = 0;  A(1, 9) = 0;   A(2, 9) = 0;	A(3, 9) = 0;  A(4, 9) = 0;
	A(0, 10) = 0; A(1, 10) = 0;  A(2, 10) = 0;	A(3, 10) = 0; A(4, 10) = 0;
	A(0, 11) = 0; A(1, 11) = 0;  A(2, 11) = 0;	A(3, 11) = 0; A(4, 11) = 0;
	A(0, 12) = 0; A(1, 12) = 0;  A(2, 12) = 0;	A(3, 12) = 0; A(4, 12) = 0;
}

double setSRT(double& srt, MatrixXd V, MatrixXd P) {
	srt = sqrt((V.transpose()*P*V)(0,0) / (N-k)); 
	return srt;
}

void setMatrixMh(MatrixXd& mh, MatrixXd Qy, double srt) {
	for (size_t i = 0; i < Qy.rows(); i++) 
	{
		for (int j = 0; j < Qy.cols(); j++) 
		{
			if (i == j)
			{
				mh(i, 0) = srt * sqrt(Qy(i, j));
			}
		}
	}	
}

void setMatrixMH(MatrixXd& mH, MatrixXd Qh, double srt) {
	for (size_t i = 0; i < Qh.rows(); i++)
	{
		mH(i, 0) = srt * sqrt(Qh(i, i));
	}
}

void getAllValues(double srt, MatrixXd V, MatrixXd mh, MatrixXd mH, MatrixXd P) {
	cout << "u = " << srt << endl;
	for (size_t i = 0; i < V.rows(); ++i)
		cout << "V" + to_string(i + 1) << " = " << V(i, 0) << endl;
	cout << "----------------------------------------------------------------------------\n";
	for (size_t i = 0; i < mh.rows(); ++i)
		cout << "mh" + to_string(i + 1) << " = " << mh(i, 0) << endl;
	cout << "----------------------------------------------------------------------------\n";
	for (size_t i = 0; i < mH.rows(); ++i)
		cout << "mH" + to_string(i + 1) << " = " << mh(i, 0) << endl;
	cout << "VT*P*V = " << (V.transpose() * P * V)(0,0) << endl;

	cout << "χ21,r,α/2 (a=5%) = " << setX21() << endl;
	cout << "χ22,r,1-α/2 (a=5%) = " << setX22() << endl;

	int N = 13;
	int r = 2;
	double p = 0.95;
	double p_a = pow(p, 1.0 / N);
	students_t myStud(r - 1);
	double t = quantile(myStud, p);

	setMatrixS();
	for (size_t i = 0; i < mH.rows(); ++i)
		cout << "S" + to_string(i + 1) << " = " << S(i, 0) << endl;
}

double setChi() {
	int r = 2;
	students_t myStud(r - 1);
	double	Pa;
	return 0;
}

int main()
{
	setlocale(LC_CTYPE, "RUS");
	cout.precision(8);


	vector <Point> points;
	ImportBenchmarks(points);

	vector <Meas> sections;
	ImportMeaseurement(sections);
	vector<double> W_;

	MatrixXd W(7, 1);
	GetResiduals(sections, W);

	MatrixXd B(7, 13);
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

	double srt = setSRT(srt, V, P);

	MatrixXd mh(13, 1);
	setMatrixMh(mh, Qy, srt);

	MatrixXd mH(5, 1);
	setMatrixMH(mh, Qh, srt);	

	getAllValues(srt, V, mh, mH, P);
	return 0;
}