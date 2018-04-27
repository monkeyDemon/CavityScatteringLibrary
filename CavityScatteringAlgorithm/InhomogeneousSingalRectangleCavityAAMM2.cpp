#include "InhomogeneousSingalRectangleCavityAAMM2.h"

InhomogeneousSingalRectangleCavityAAMM2::InhomogeneousSingalRectangleCavityAAMM2(unsigned int cavityType) :InhomogeneousSingalRectangleCavity(cavityType)
{
	int test = 0;
}

complex<double> InhomogeneousSingalRectangleCavityAAMM2::f(double x, double y)
{
	//Դ����ȡֵ���������Ҫ��-1��
	complex<double> u = exp(x*y)*sin(0.5*k0*x)*sin((0.5*k0 + 0.25*M_PI)*y);
	complex<double> f = (x*x + y*y - 0.5*k0*k0 - M_PI*M_PI / 16.0 - 0.25*k0*M_PI + k0*k0*epr)*u + k0*y*exp(x*y)*cos(0.5*k0*x)*sin((0.5*k0 + 0.25*M_PI)*y) + (k0 + 0.5*M_PI)*x*exp(x*y)*sin(0.5*k0*x)*cos((0.5*k0 + 0.25*M_PI)*y);
	f = -f;
	return f;
}

complex<double> InhomogeneousSingalRectangleCavityAAMM2::f_(double x, double y)
{
	//�Ǿ��Ƚ���Դ����ȡֵ���������Ҫ��-1��
	complex<double> u = exp(x*y)*sin(0.5*k0*x)*sin((0.5*k0 + 0.25*M_PI)*y);
	complex<double> f_ = (x*x + y*y - 0.5*k0*k0 - M_PI*M_PI / 16.0 - 0.25*k0*M_PI + k0*k0*epr_int)*u + k0*y*exp(x*y)*cos(0.5*k0*x)*sin((0.5*k0 + 0.25*M_PI)*y) + (k0 + 0.5*M_PI)*x*exp(x*y)*sin(0.5*k0*x)*cos((0.5*k0 + 0.25*M_PI)*y);
	f_ = -f_;
	return f_;
}

complex<double> InhomogeneousSingalRectangleCavityAAMM2::u(double x, double y)
{
	//ǻ���ھ�ȷ���ֵ�����ڻ�ȡ�߽������򿼲��㷨��ȷ��
	double real = exp(x*y)*sin(0.5*k0*x)*sin((0.5*k0 + 0.25*M_PI)*y);
	complex<double> u(real, 0);
	return u;
}

complex<double> InhomogeneousSingalRectangleCavityAAMM2::u_(double x, double y)
{
	//ǻ���ڷǾ��Ƚ��ʾ�ȷ���ֵ�����ڻ�ȡ�߽������򿼲��㷨��ȷ��
	double real = exp(x*y)*sin(0.5*k0*x)*sin((0.5*k0 + 0.25*M_PI)*y);
	complex<double> u_(real, 0);
	return u_;
}


VectorXcd InhomogeneousSingalRectangleCavityAAMM2::compute_g(MatrixXcd &G, vector<vector<double>> &nbound)
{
	int apertureNum = nbound[0].size();

	VectorXd uu(apertureNum);
	VectorXd UyN1(apertureNum);
	for (int k = 0; k < apertureNum; k++)
	{
		double x = nbound[1][k];

		// AAMM 2
		double temp_uuk = exp(x)*sin(0.5*k0*x)*sin(0.5*k0 + 0.25*M_PI);
		uu(k) = temp_uuk;
		UyN1(k) = x*temp_uuk + (0.5*k0 + 0.25*M_PI)*exp(x)*sin(0.5*k0*x)*cos(0.5*k0 + 0.25*M_PI);
	}
	VectorXcd g(apertureNum);
	g = UyN1 - G * uu;
	return g;
}