#pragma once
#include "headers.h"

void ImportBenchmarks(vector<Point>&);

void ImportMeaseurement(vector<Meas>&);

MatrixXd GetResiduals(vector <Meas>&, MatrixXd&);

void setMatrixB(MatrixXd&);

void setMatrixP(MatrixXd&, vector<Meas>);

void setMatrixW_h(MatrixXd&, vector<Meas>);

void setMatrixA(MatrixXd&);

double setu(double&, MatrixXd, MatrixXd);

void setMatrixMh(MatrixXd&, MatrixXd, double);

void setMatrixMH(MatrixXd&, MatrixXd, double);

//���� ������� ��� getAllValues
void getAnyVerticalMatrix(MatrixXd, string);

void getMatrixS(MatrixXd, MatrixXd, double);

double setXi21(double);

double setXi22(double);

void getXi();

double set_t(int, double);

void getTao();

void getAllValues(double, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd);