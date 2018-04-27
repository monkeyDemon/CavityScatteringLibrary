#include"ApertureIntegral.h"


double ApertureIntegral::r[6] = { -0.4900604943e13, 0.1275274390e13, -0.5153438139e11, 0.7349264551e9, -0.4237922726e7, 0.8511937935e4 };
double ApertureIntegral::s[7] = { 0.2499580570e14, 0.4244419664e12, 0.3733650367e10, 0.2245904002e8, 0.1020426050e6, 0.3549632885e3, 1.0 };
double ApertureIntegral::p[5] = { 1.0, 0.183105e-2, -0.3516396496e-4, 0.2457520174e-5, -0.240337019e-6 };
double ApertureIntegral::q[5] = { 0.04687499995, -0.2002690873e-3, 0.8449199096e-5, -0.88228987e-6, 0.105787412e-6 };


double ApertureIntegral::Jr[6] = { 72362614232.0, -7895059235.0, 242396853.1, -2972611.439, 15704.48260, -30.16036606 };
double ApertureIntegral::Js[6] = { 144725228442.0, 2300535178.0, 18583304.74, 99447.43394, 376.9991397, 1.0 };
double ApertureIntegral::Jp[5] = { 1.0, 0.183105e-2, -0.3516396496e-4, 0.2457520174e-5, -0.240337019e-6 };
double ApertureIntegral::Jq[5] = { 0.04687499995, -0.2002690873e-3, 0.8449199096e-5, -0.88228987e-6, 0.105787412e-6 };


MatrixXcd ApertureIntegral::computeG(int n, double wavek, double a, double b)
{
	vector<double> BMRE2;
	vector<double> BMIM;
	SubMatrixG(n, wavek, a, b, BMRE2, BMIM);

	int Gsize = BMRE2.size();
	MatrixXcd G = toeplitz(BMRE2, BMIM);
	return G;
}


void ApertureIntegral::SubMatrixG(int n, double wavek, double a, double b, vector<double> &BMRE2, vector<double> &BMIM)
{
	vector<double> BMRE;
	BMRE.resize(n);
	BMRE2.resize(n);
	BMIM.resize(n);

	double h = (b - a) / (n + 1);

	vector<double> P;
	P.resize(n + 2);
	for (int ind = 0; ind < P.size(); ind++)
		P[ind] = a + ind*h;

	vector<double> col;
	col.resize(n);
	for (int ind = 0; ind < col.size(); ind++)
		col[ind] = P[ind + 1];

	for (int j = 0; j < n; j++)
	{
		double m1 = 1;
		double subh = 2 * h;
		double sum2 = 0;
		m1 = m1 * 2;
		subh = subh * 0.5;
		double sum1 = sum2;
		sum2 = 0;
		for (int k = 1; k <= m1 - 1; k += 2)
		{
			double xx = P[j] + 2 * k * h / m1;
			double fi = phi(j, n, P, xx);
			double fenmu = abs(xx - P[1]);
			double wjmtemp;
			if (fenmu > 1e-12)
			{
				double xxx = wavek * fenmu;
				double bj11;
				if (j == 0)
				{
					bj11 = Bessel_mdf1Y1(xxx);
				}
				else
				{
					bj11 = Bessel_mdf2Y1(xxx);
				}
				wjmtemp = bj11 / fenmu;
			}
			else
			{
				if (j == 0)
					wjmtemp = wavek * (-0.196057111305026);
				else
					wjmtemp = 0;
			}
			sum2 = sum2 + wjmtemp * fi;
		}
		sum2 = 0.5 * sum1 + subh * sum2;
		while (abs(sum1 - sum2) > pow(h, 3))
		{
			m1 = m1 * 2;
			subh = subh * 0.5;
			sum1 = sum2;
			sum2 = 0;
			for (int k = 1; k <= m1 - 1; k += 2)
			{
				double xx = P[j] + 2 * k * h / m1;
				double fi = phi(j, n, P, xx);
				double fenmu = abs(xx - P[1]);
				double wjmtemp;
				if (fenmu > 1e-12)
				{
					double xxx = wavek * fenmu;
					double bj11;
					if (j == 0)
						bj11 = Bessel_mdf1Y1(xxx);
					else
						bj11 = Bessel_mdf2Y1(xxx);
					wjmtemp = bj11 / fenmu;
				}
				else
				{
					if (j == 0)
						wjmtemp = wavek * (-0.196057111305026);
					else
						wjmtemp = 0;
				}
				sum2 = sum2 + wjmtemp * fi;
			}
			sum2 = 0.5 * sum1 + subh * sum2;
		}
		BMRE[j] = sum2 * 0.5 * wavek;
	}
	BMRE[0] += pow(wavek, 2) * h * (log(wavek * h) - 1.5) / (2 * M_PI);

	//       2.3 Compute BMIM(1, j)
	for (int j = 0; j < n; j++)
	{
		double m1 = 1;
		double subh = 2 * h;
		double sum2 = 0;
		m1 = m1 * 2;
		subh = subh * 0.5;
		double sum1 = sum2;
		sum2 = 0;
		for (int k = 1; k <= m1 - 1; k += 2)
		{
			double xx = P[j] + 2 * k * h / m1;
			double fi = phi(j, n, P, xx);
			double fenmu = abs(xx - P[1]);
			double wjmtemp;
			if (fenmu > 1e-12)
			{
				double xxx = wavek * fenmu;
				double bj11 = Bessel_J1(xxx);
				wjmtemp = bj11 / fenmu;
			}
			else
			{
				wjmtemp = 0.5 * wavek;
			}
			sum2 = sum2 + wjmtemp * fi;
		}
		sum2 = 0.5 * sum1 + subh * sum2;
		while (abs(sum1 - sum2) > pow(h, 3))
		{
			m1 = m1 * 2;
			subh = subh * 0.5;
			sum1 = sum2;
			sum2 = 0;
			for (int k = 1; k <= m1 - 1; k += 2)
			{
				double xx = P[j] + 2 * k * h / m1;
				double fi = phi(j, n, P, xx);
				double fenmu = abs(xx - P[1]);
				double wjmtemp;
				if (fenmu > 1e-12)
				{
					double xxx = wavek * fenmu;
					double bj11 = Bessel_J1(xxx);
					wjmtemp = bj11 / fenmu;
				}
				else
				{
					wjmtemp = 0.5 * wavek;
				}
				sum2 = sum2 + wjmtemp * fi;
			}
			sum2 = 0.5 * sum1 + subh * sum2;
		}
		BMIM[j] = -sum2 * 0.5 * wavek;
	}

	//  2.5 Compute the hypersingular part
	BMRE2 = BI1_Algorithm_III(n, P, col, BMRE);
}


MatrixXcd ApertureIntegral::toeplitz(vector<double> BMRE2, vector<double> BMIM)
{
	int Gsize = BMRE2.size();
	MatrixXcd G(Gsize, Gsize);

	for (int i = 0; i < Gsize; i++)
	{
		// fill BMRE2[i] and BMIM[i]
		double curReal = -BMRE2[i];
		double curImg = -BMIM[i];

		int iterNum = Gsize - i;
		int rowIndex = 0;
		int colIndex = i;
		while (iterNum>0)
		{
			G(rowIndex, colIndex) = complex<double>(curReal, curImg);
			iterNum--;
			rowIndex++;
			colIndex++;
		}

		iterNum = Gsize - i;
		rowIndex = i;
		colIndex = 0;
		while (iterNum>0)
		{
			G(rowIndex, colIndex) = complex<double>(curReal, curImg);
			iterNum--;
			rowIndex++;
			colIndex++;
		}
	}

	return G;
}

double ApertureIntegral::phi(int m, int n, vector<double> P, double x)
{
	double fi;
	if (x > P[m] && x < P[m + 1])
	{
		fi = (x - P[m]) / (P[m + 1] - P[m]);
	}
	else
	{
		if (x >= P[m + 1] && x < P[m + 2])
			fi = (P[m + 2] - x) / (P[m + 2] - P[m + 1]);
		else
			fi = 0;
	}
	return fi;
}

double ApertureIntegral::Bessel_J1(double x)
{
	double bj1;

	double ax = abs(x);
	if (ax < 8)
	{
		double y = pow(x, 2);
		double t = x*(Jr[0] + y*(Jr[1] + y*(Jr[2] + y*(Jr[3] + y*(Jr[4] + y*Jr[5])))));
		bj1 = t / (Js[0] + y*(Js[1] + y*(Js[2] + y*(Js[3] + y*(Js[4] + y*Js[5])))));
	}
	else
	{
		double z = 8 / x;
		double y = pow(z, 2);
		double xx = ax - 0.75 * M_PI;
		double t1 = sqrt(2 / ax / M_PI);
		double t2 = cos(xx)*(Jp[0] + y*(Jp[1] + y*(Jp[2] + y*(Jp[3] + y*Jp[4])))) - z*sin(xx)*(Jq[0] + y*(Jq[1] + y*(Jq[2] + y*(Jq[3] + y*Jq[4]))));
		int sgn;  //作用相当于matlab中的sign(x)
		if (x > 1e-10)
			sgn = 1;
		else if (x < -1e-10)
			sgn == -1;
		else
			sgn = 0;
		bj1 = t1 * t2 * sgn;
	}
	return bj1;
}


double ApertureIntegral::Bessel_mdf1Y1(double x)
{
	double by1;
	if (x < 8)
	{
		double bj11 = Bessel_J1(x);
		double y = pow(x,2);
		double t = x * (r[0] + y * (r[1] + y * (r[2] + y * (r[3] + y * (r[4] + y * r[5])))));
		by1 = t / (s[0] + y * (s[1] + y * (s[2] + y * (s[3] + y * (s[4] + y * (s[5] + s[6] * y))))));
		if (x > 1e-12)
		{
			by1 = by1 + 0.636619772367581 * (log(x) * (bj11 - 0.5 * x));  //log以e为底
		}
	}
	else
	{
		double z = 8 / x;
		double y = pow(z, 2);
		double xx = x - 2.35619449019234;
		double t1 = sqrt(0.636619772367581 / x) * sin(xx)*(p[0] + y*(p[1] + y*(p[2] + y*(p[3] + y*p[4]))));
		double t2 = z*cos(xx)*(q[0] + y*(q[1] + y*(q[2] + y*(q[3] + y*q[4]))));
		by1 = t1 + t2 + (2 / x - x*log(x)) / M_PI;    
	}
	return by1;
}

double ApertureIntegral::Bessel_mdf2Y1(double x)
{
	double by1;
	if (x < 8)
	{
		double bj11 = Bessel_J1(x);
		double y = pow(x, 2);
		double t = x * (r[0] + y * (r[1] + y * (r[2] + y * (r[3] + y * (r[4] + y * r[5])))));
		by1 = t / (s[0] + y * (s[1] + y * (s[2] + y * (s[3] + y * (s[4] + y * (s[5] + s[6] * y))))));
		if (x > 1e-12)
		{
			by1 = by1 + 0.636619772367581*(log(x)*(bj11));   
		}
	}
	else
	{
		double z = 8 / x;
		double y = pow(z, 2);
		double xx = x - 2.35619449019234;
		double t1 = sqrt(0.636619772367581 / x) * sin(xx)*(p[0] + y*(p[1] + y*(p[2] + y*(p[3] + y*p[4]))));
		double t2 = z*cos(xx)*(q[0] + y*(q[1] + y*(q[2] + y*(q[3] + y*q[4]))));
		by1 = t1 + t2 + 0.636619772367581 / x;     
	}
	return by1;
}

vector<double> ApertureIntegral::BI1_Algorithm_III(int n, vector<double> P, vector<double> col, vector<double> BMRE)
{
	double h = P[2] - P[1];

	vector<double> work;
	work.resize(n);
	vector<double> BMRE2;
	BMRE2.resize(n);

	work[0] = 2 / (M_PI * h);
	work[1] = -(1 - log(2)) / (M_PI * h);
	if (n > 2)
	{
		for (int i = 3; i<=n;i++)
			work[i-1] = -log(pow(i - 1,2) / (pow(i - 1, 2) - 1)) / (M_PI * h);
	}

	for (int i = 0; i < n; i++)
		BMRE2[i] = BMRE[i] + work[i];

	return BMRE2;
}