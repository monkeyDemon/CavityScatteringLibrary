#pragma once
#include "InhomogeneousSingalRectangleCavity.h"

class  _declspec(dllexport) InhomogeneousSingalRectangleCavityAAMM2 : public InhomogeneousSingalRectangleCavity
{
public:
	InhomogeneousSingalRectangleCavityAAMM2(unsigned int cavityType);

protected:
	complex<double> f(double x, double y);
	complex<double> f_(double x, double y);
	complex<double> u(double x, double y);
	complex<double> u_(double x, double y);

	VectorXcd compute_g(MatrixXcd &G, vector<vector<double>> &nbound);

	//新增一个具体算例时，需要考虑以下方法是否需要覆写
	/*
	Matrix2d beta(double x, double y);
	Matrix2d beta_(double x, double y);
	complex<double> q(double x, double y);
	complex<double> q_(double x, double y);
	complex<double> f(double x, double y);
	complex<double> f_(double x, double y);
	complex<double> u(double x, double y);
	complex<double> u_(double x, double y);
	double a(double x, double y);
	Vector2d b(double x, double y);
	*/
};