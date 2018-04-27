#include "InhomogeneousSingalRectangleCavityAAMM1.h"

InhomogeneousSingalRectangleCavityAAMM1::InhomogeneousSingalRectangleCavityAAMM1(unsigned int cavityType) :InhomogeneousSingalRectangleCavity(cavityType)
{
	int test = 0;
}

complex<double> InhomogeneousSingalRectangleCavityAAMM1::f(double x, double y)
{
	//源函数取值（计算后需要乘-1）
	complex<double> f = (2.0 + (k2 - 2.0*M_PI*M_PI) * (y - 1 / 2.0)*(y - 1 / 2.0))*sin(M_PI*x)*cos(M_PI*y) - M_PI*4.0 *(y - 1 / 2.0)*sin(M_PI*x)* sin(M_PI*y);
	f = -f;
	return f;
}

complex<double> InhomogeneousSingalRectangleCavityAAMM1::f_(double x, double y)
{
	//非均匀介质源函数取值（计算后需要乘-1）
	complex<double> f_ = (2.0 + (k2_int - 2.0*M_PI*M_PI) * (y - 1 / 2.0)*(y - 1 / 2.0))*sin(M_PI*x)*sin(M_PI*y) + M_PI*4.0 *(y - 1 / 2.0)*sin(M_PI*x)* cos(M_PI*y);
	f_ = -f_;
	return f_;
}

complex<double> InhomogeneousSingalRectangleCavityAAMM1::u(double x, double y)
{
	//腔体内精确解的值，用于获取边界条件或考察算法精确度
	
	double real = (y - 1 / 2.0)*(y - 1 / 2.0)* sin(M_PI* x) * cos(M_PI*y);
	complex<double> u(real, 0);
	return u;
}

complex<double> InhomogeneousSingalRectangleCavityAAMM1::u_(double x, double y)
{
	//腔体内非均匀介质精确解的值，用于获取边界条件或考察算法精确度

	double real = (y - 1 / 2.0)*(y - 1 / 2.0)* sin(M_PI* x) * sin(M_PI*y);
	complex<double> u_(real, 0);
	return u_;
}


double InhomogeneousSingalRectangleCavityAAMM1::a(double x, double y)
{
	double u = (y - 1 / 2.0)*(y - 1 / 2.0)* sin(M_PI* x) * cos(M_PI*y);
	double u_int = (y - 1 / 2.0)*(y - 1 / 2.0)* sin(M_PI* x) * sin(M_PI*y);
	double out = u - u_int;
	return out;
}

Vector2d InhomogeneousSingalRectangleCavityAAMM1::b(double x, double y)
{
	double b1 = 0; 
	double b2 = 0; 

	double uy = 2 * (y - 0.5)*sin(M_PI*x)*cos(M_PI*y) - M_PI*(y - 0.5)*(y - 0.5)*sin(M_PI*x)*sin(M_PI*y);
	double uy_int = 2 * (y - 0.5)*sin(M_PI*x)*sin(M_PI*y) + M_PI*(y - 0.5)*(y - 0.5)*sin(M_PI*x)*cos(M_PI*y);
	b2 = uy - uy_int;

	Vector2d out(b1, b2);

	return out;
}


VectorXcd InhomogeneousSingalRectangleCavityAAMM1::compute_g(MatrixXcd &G, vector<vector<double>> &nbound)
{
	int apertureNum = nbound[0].size();

	VectorXd uu(apertureNum);
	VectorXd UyN1(apertureNum);
	for (int k = 0; k < apertureNum; k++)
	{
		double x = nbound[1][k];
		// AAMM 1
		uu(k) = 1 / 4.0 * sin(x*M_PI);
		UyN1(k) = -sin(M_PI*x);

		// AAMM 2
		// uu(k, 1) = exp(x)*sin(0.5*k0*x)*sin(0.5*k0 + 0.25*pi);
		// UyN1(k, 1) = x*uu(k, 1) + (0.5*k0 + 0.25*pi)*exp(x)*sin(0.5*k0*x)*cos(0.5*k0 + 0.25*pi);
	}
	VectorXcd g(apertureNum);
	g = UyN1 - G * uu;
	return g;
}