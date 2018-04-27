#pragma once

#include<Eigen/Dense>
//#include <Eigen/Sparse>

#define _USE_MATH_DEFINES //若不使用此宏，math.h中无法找到π的定义M_PI
#include<math.h>
#include<vector>

using namespace std;
using namespace Eigen;

class ApertureIntegral
{
public:
	static MatrixXcd computeG(int n, double wavek, double a, double b);
	
private:
	static double r[6];
	static double s[7];
	static double p[5];
	static double q[5];

	static double Jr[6];   //use in function Bessel_J1
	static double Js[6];
	static double Jp[5];
	static double Jq[5];

	static void SubMatrixG(int n, double wavek, double a, double b, vector<double> &BMRE2, vector<double> &BMIM);
	static MatrixXcd toeplitz(vector<double> BMRE2, vector<double> BMIM);

	static double phi(int m, int n, vector<double> p, double x);
	static double Bessel_J1(double x);
	static double Bessel_mdf1Y1(double x);
	static double Bessel_mdf2Y1(double x);
	static vector<double> BI1_Algorithm_III(int n, vector<double> p, vector<double> col, vector<double> BMRE);
};