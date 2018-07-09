#include "SingalCavity.h"


SingalCavity::SingalCavity(unsigned int cavityType) :Cavity(cavityType)
{
	int test = 0;
}

void SingalCavity::InitElectromagneticParameter(double k0, complex<double> epr, double theta)
{
	this->k0 = k0;
	this->epr = epr;
	this->k2 = k0*k0*epr;
	this->theta = theta;

	//Mark function InitElectromagneticParameter has been executed
	this->initalCheckKey = this->initalCheckKey | 4; // 4 means 0000 0100
}

void SingalCavity::InitAperture(double apertureLeft, double apertureRight, double apertureY)
{
	//set the parameter of aperture
	this->apertureLeft = apertureLeft;
	this->apertureRight = apertureRight;
	this->apertureY = apertureY;

	//Mark function InitCavityShapeParameter has been executed
	this->initalCheckKey = this->initalCheckKey | 8; // 8 means 0000 1000
}



bool SingalCavity::Solve()
{
	char *checkLog = "";
	bool checkResult;
	checkResult = InitialCheck(checkLog);

	if (checkResult == false)
	{
		printf(checkLog);
		return false;
	}

	MyTimer myTimer(1);

	myTimer.Start("setTri");
	TriangleMesh U(this->meshWidth, this->meshHeight);
	TriangleMesh L(this->meshWidth, this->meshHeight);
	/*int nn;
	vector<int> nu;
	vector<vector<double>> nbound;*/
	setTri(U, L, nn, nu, nbound);
	myTimer.EndAndPrint();


	this->G_aperture = ApertureIntegral::computeG(nbound[0].size(), this->k0, this->apertureLeft, this->apertureRight);

	this->g_aperture = SingalCavity::compute_g(G_aperture,nbound);

	myTimer.Start("setGrid");
	vector<vector<gridCell>> gridCell = setGrid(U, L, nbound, nu);
	myTimer.EndAndPrint();

	myTimer.Start("setRightHand");
	VectorXcd rh = setRightHand(U, L, nu);
	myTimer.EndAndPrint();

	myTimer.Start("setA");
	//SparseMatrix<complex<double>> A = setA(nn, nu, nbound, gridCell);
	mxArray * mx_A = setA_mx(nn, nu, nbound, gridCell);
	myTimer.EndAndPrint();

	myTimer.Start("setB");
	//VectorXcd B = setB(gridCell, rh, nn, nu);
	mxArray *mx_B = setB_mx(gridCell, rh, nn, nu);
	myTimer.EndAndPrint();

	myTimer.Start("solveX");
	//VectorXcd x = solveX(A, B);
	VectorXcd x = solveX_mx(mx_A, mx_B);
	myTimer.EndAndPrint();
	//cout << x << endl;

	myTimer.Start("assign");
	assign(x, U, L, nu);
	myTimer.EndAndPrint();

	myTimer.Start("getAperture");
	this->solutionOfAperture = getAperture(nbound, U, L);
	myTimer.EndAndPrint();
	//cout << solutionOfAperture << endl;
}


void SingalCavity::PlotTriangleMesh(string title, string xlabel, string ylabel)
{
	char *checkLog = "";
	bool checkResult;
	checkResult = InitialCheck(checkLog);

	if (checkResult == false)
	{
		printf(checkLog);
		throw exception("初始化检查未通过！");
		return;
	}

	MyTimer myTimer(1);

	// =====================================================
	myTimer.Start("setTri");
	TriangleMesh U(this->meshWidth, this->meshHeight);
	TriangleMesh L(this->meshWidth, this->meshHeight);

	setTri(U, L, nn, nu, nbound);
	myTimer.EndAndPrint();


	// =====================================================
	myTimer.Start("绘制三角形网格");
	//调用matalb引擎画图
	engEvalString(ep, "figure");
	engEvalString(ep, "hold on;");

	// 定义matlab数组mxArray
	mxArray *mx_xValue, *mx_yValue;

	//将C中数组转化为matlab数组mxArray
	mx_xValue = mxCreateDoubleMatrix(1, 4, mxREAL);
	mx_yValue = mxCreateDoubleMatrix(1, 4, mxREAL);

	// 定义matlab数组的数据指针
	double *xValuePr,  *yValuePr;

	int m = this->meshWidth;
	int n = this->meshHeight;
	for (int j = 0; j < m; j++)
	{
		for (int l = 0; l < n; l++)
		{
			// 处理U
			if (U.Get_sign(j, l))
			{
				// normal triangle
				Triangle_Normal a = U.Get_normal(j, l);
				if (a.sgn == -10)
				{
					//plot([a.x, a.x(1)], [a.y, a.y(1)], 'b');
					xValuePr = mxGetPr(mx_xValue);
					*xValuePr++ = a.x(0); *xValuePr++ = a.x(1); *xValuePr++ = a.x(2); *xValuePr++ = a.x(0);
					yValuePr = mxGetPr(mx_yValue);
					*yValuePr++ = a.y(0); *yValuePr++ = a.y(1); *yValuePr++ = a.y(2); *yValuePr++ = a.y(0);
					engPutVariable(ep, "xValue", mx_xValue); engPutVariable(ep, "yValue", mx_yValue);
					engEvalString(ep, "plot(xValue, yValue, 'b');");
				}
			}
			else
			{
				// all triangle
				Triangle_All a = U.Get_all(j, l);
				if (abs(a.sgn) < 10)
				{
					if (a.sgn > 0)
					{
						//plot([a.x2, a.x3, a.x4, a.x2], [a.y2, a.y3, a.y4, a.y2], 'b');
						xValuePr = mxGetPr(mx_xValue);
						*xValuePr++ = a.x2; *xValuePr++ = a.x3; *xValuePr++ = a.x4; *xValuePr++ = a.x2;
						yValuePr = mxGetPr(mx_yValue);
						*yValuePr++ = a.y2; *yValuePr++ = a.y3; *yValuePr++ = a.y4; *yValuePr++ = a.y2;
						engPutVariable(ep, "xValue", mx_xValue); engPutVariable(ep, "yValue", mx_yValue);
						engEvalString(ep, "plot(xValue, yValue, 'b');");

						//plot([a.x3, a.x5, a.x4, a.x3], [a.y3, a.y5, a.y4, a.y3], 'b');
						xValuePr = mxGetPr(mx_xValue);
						*xValuePr++ = a.x3; *xValuePr++ = a.x5; *xValuePr++ = a.x4; *xValuePr++ = a.x3;
						yValuePr = mxGetPr(mx_yValue);
						*yValuePr++ = a.y3; *yValuePr++ = a.y5; *yValuePr++ = a.y4; *yValuePr++ = a.y3;
						engPutVariable(ep, "xValue", mx_xValue); engPutVariable(ep, "yValue", mx_yValue);
						engEvalString(ep, "plot(xValue, yValue, 'b');");
					}
					else
					{
						//plot([a.x1, a.x4, a.x5, a.x1], [a.y1, a.y4, a.y5, a.y1], 'b');
						xValuePr = mxGetPr(mx_xValue);
						*xValuePr++ = a.x1; *xValuePr++ = a.x4; *xValuePr++ = a.x5; *xValuePr++ = a.x1;
						yValuePr = mxGetPr(mx_yValue);
						*yValuePr++ = a.y1; *yValuePr++ = a.y4; *yValuePr++ = a.y5; *yValuePr++ = a.y1;
						engPutVariable(ep, "xValue", mx_xValue); engPutVariable(ep, "yValue", mx_yValue);
						engEvalString(ep, "plot(xValue, yValue, 'b');");
					}
				}
				if (abs(a.sgn) > 10)
				{
					if (a.sgn > 0)
					{
						//plot([a.x2, a.x3, a.p5(1), a.x2], [a.y2, a.y3, a.p5(2), a.y2], 'b');
						xValuePr = mxGetPr(mx_xValue);
						*xValuePr++ = a.x2; *xValuePr++ = a.x3; *xValuePr++ = a.p5(0); *xValuePr++ = a.x2;
						yValuePr = mxGetPr(mx_yValue);
						*yValuePr++ = a.y2; *yValuePr++ = a.y3; *yValuePr++ = a.p5(1); *yValuePr++ = a.y2;
						engPutVariable(ep, "xValue", mx_xValue); engPutVariable(ep, "yValue", mx_yValue);
						engEvalString(ep, "plot(xValue, yValue, 'b');");
					}
					else
					{
						//plot([a.x1, a.x2, a.p5(1), a.x1], [a.y1, a.y2, a.p5(2), a.y1], 'b');
						xValuePr = mxGetPr(mx_xValue);
						*xValuePr++ = a.x1; *xValuePr++ = a.x2; *xValuePr++ = a.p5(0); *xValuePr++ = a.x1;
						yValuePr = mxGetPr(mx_yValue);
						*yValuePr++ = a.y1; *yValuePr++ = a.y2; *yValuePr++ = a.p5(1); *yValuePr++ = a.y1;
						engPutVariable(ep, "xValue", mx_xValue); engPutVariable(ep, "yValue", mx_yValue);
						engEvalString(ep, "plot(xValue, yValue, 'b');");
					}
				}
			}

			// 处理L
			if (L.Get_sign(j, l))
			{
				// normal triangle
				Triangle_Normal a = L.Get_normal(j, l);
				if (a.sgn == -10)
				{
					//plot([a.x, a.x(1)], [a.y, a.y(1)], 'b');
					xValuePr = mxGetPr(mx_xValue);
					*xValuePr++ = a.x(0); *xValuePr++ = a.x(1); *xValuePr++ = a.x(2); *xValuePr++ = a.x(0);
					yValuePr = mxGetPr(mx_yValue);
					*yValuePr++ = a.y(0); *yValuePr++ = a.y(1); *yValuePr++ = a.y(2); *yValuePr++ = a.y(0);
					engPutVariable(ep, "xValue", mx_xValue); engPutVariable(ep, "yValue", mx_yValue);
					engEvalString(ep, "plot(xValue, yValue, 'b');");
				}
			}
			else
			{
				// all triangle
				Triangle_All a = L.Get_all(j, l);
				if (abs(a.sgn) < 10)
				{
					if (a.sgn > 0)
					{
						//plot([a.x2, a.x3, a.x4, a.x2], [a.y2, a.y3, a.y4, a.y2], 'b');
						xValuePr = mxGetPr(mx_xValue);
						*xValuePr++ = a.x2; *xValuePr++ = a.x3; *xValuePr++ = a.x4; *xValuePr++ = a.x2;
						yValuePr = mxGetPr(mx_yValue);
						*yValuePr++ = a.y2; *yValuePr++ = a.y3; *yValuePr++ = a.y4; *yValuePr++ = a.y2;
						engPutVariable(ep, "xValue", mx_xValue); engPutVariable(ep, "yValue", mx_yValue);
						engEvalString(ep, "plot(xValue, yValue, 'b');");

						//plot([a.x3, a.x5, a.x4, a.x3], [a.y3, a.y5, a.y4, a.y3], 'b');
						xValuePr = mxGetPr(mx_xValue);
						*xValuePr++ = a.x3; *xValuePr++ = a.x5; *xValuePr++ = a.x4; *xValuePr++ = a.x3;
						yValuePr = mxGetPr(mx_yValue);
						*yValuePr++ = a.y3; *yValuePr++ = a.y5; *yValuePr++ = a.y4; *yValuePr++ = a.y3;
						engPutVariable(ep, "xValue", mx_xValue); engPutVariable(ep, "yValue", mx_yValue);
						engEvalString(ep, "plot(xValue, yValue, 'b');");
					}
					else
					{
						//plot([a.x1, a.x4, a.x5, a.x1], [a.y1, a.y4, a.y5, a.y1], 'b');
						xValuePr = mxGetPr(mx_xValue);
						*xValuePr++ = a.x1; *xValuePr++ = a.x4; *xValuePr++ = a.x5; *xValuePr++ = a.x1;
						yValuePr = mxGetPr(mx_yValue);
						*yValuePr++ = a.y1; *yValuePr++ = a.y4; *yValuePr++ = a.y5; *yValuePr++ = a.y1;
						engPutVariable(ep, "xValue", mx_xValue); engPutVariable(ep, "yValue", mx_yValue);
						engEvalString(ep, "plot(xValue, yValue, 'b');");
					}
				}
				if (abs(a.sgn) > 10)
				{
					if (a.sgn > 0)
					{
						//plot([a.x2, a.x3, a.p5(1), a.x2], [a.y2, a.y3, a.p5(2), a.y2], 'b');
						xValuePr = mxGetPr(mx_xValue);
						*xValuePr++ = a.x2; *xValuePr++ = a.x3; *xValuePr++ = a.p5(0); *xValuePr++ = a.x2;
						yValuePr = mxGetPr(mx_yValue);
						*yValuePr++ = a.y2; *yValuePr++ = a.y3; *yValuePr++ = a.p5(1); *yValuePr++ = a.y2;
						engPutVariable(ep, "xValue", mx_xValue); engPutVariable(ep, "yValue", mx_yValue);
						engEvalString(ep, "plot(xValue, yValue, 'b');");
					}
					else
					{
						//plot([a.x1, a.x2, a.p5(1), a.x1], [a.y1, a.y2, a.p5(2), a.y1], 'b');
						xValuePr = mxGetPr(mx_xValue);
						*xValuePr++ = a.x1; *xValuePr++ = a.x2; *xValuePr++ = a.p5(0); *xValuePr++ = a.x1;
						yValuePr = mxGetPr(mx_yValue);
						*yValuePr++ = a.y1; *yValuePr++ = a.y2; *yValuePr++ = a.p5(1); *yValuePr++ = a.y1;
						engPutVariable(ep, "xValue", mx_xValue); engPutVariable(ep, "yValue", mx_yValue);
						engEvalString(ep, "plot(xValue, yValue, 'b');");
					}
				}
			}

		}
	}


	engEvalString(ep, "xlabel('x', 'FontSize', 12);");
	engEvalString(ep, "ylabel('y', 'FontSize', 12);");

	myTimer.EndAndPrint();
}




//---------------------------------protect---------------------------------


Matrix2d SingalCavity::beta(double x, double y)
{
	Matrix2d beta;
	beta(0, 0) = 1;
	beta(1, 0) = 0;
	beta(0, 1) = 0;
	beta(1, 1) = 1;
	return beta;
}

complex<double> SingalCavity::q(double x, double y)
{
	//here gives three different implementations
	//Need to check the correctness

	/*complex<double> q;
	q._Val[0] = -1 * k2.real();
	q._Val[1] = -1 * k2.imag();

	complex<double> q(this->k2 *(complex<double>)(-1));*/
	
	complex<double> q(-k2);

	return q;
}

complex<double> SingalCavity::f(double x, double y)
{
	complex<double> f(0, 0);
	return f;
}

complex<double> SingalCavity::u(double x, double y)
{
	//如果子类对应的具体算例有精确解，需要覆写本方法
	//否则本方法只会在获取腔体边界时用到，由于腔体是完全良导体，因此只需返回0
	complex<double> u(0, 0);
	return u;
}


void SingalCavity::setTri(TriangleMesh &U, TriangleMesh &L, int &nn, vector<int> &nu, vector<vector<double>> &nbound)
{
	int m = this->meshWidth;
	int n = this->meshHeight;
	
	double dx = this->stepX;
	double dy = this->stepY;

	double left = this->virtualBorderLeft;
	double bottom = this->virtualBorderBottom;
	
	// 分别计算腔体内非口径面与口径面上的变量个数
	int notApertureNum = 0;
	int ApertureNum = 0;
	for (int j = 1; j <= m; j++)
	{
		for (int k = 1; k <= n; k++)
		{
			double x4 = left + j*dx;
			double y4 = bottom + k*dy;
			double a4 = phi(x4, y4);
			if (a4 < -1e-10 && j < m)
			{
				if (k != n)
					notApertureNum++;
				else
					ApertureNum++;
			}
		}
	}
	nn = ApertureNum + notApertureNum;

	nu.resize(n*(m - 1));

	/*nbound.resize(ApertureNum);
	for (int i = 0; i < nbound.size(); i++)
	nbound[i].resize(2);*/
	nbound.resize(2);
	nbound[0].resize(ApertureNum);
	nbound[1].resize(ApertureNum);
	/*for (int i = 0; i < nbound.size(); i++)
		nbound[i].reserve(2);*/


	int ApertureIndex = 0;    // nn_aperture
	int notApertureIndex = 0;
	
	double x1, y1, x2, y2, x3, y3, x4, y4;
	double a1, a2, a3, a4;
	int id;
	for (int j = 0; j < m; j++)
	{
		for (int k = 0; k < n; k++)
		{
			// 3      4
			// ------- （j, k）
			// | \ U |
			// |  \  |
			// | L \ |
			// |    \|
			// -------
			// 1      2

			x1 = left + j*dx;
			y1 = bottom + k*dy;
			x2 = left + (j + 1)*dx;
			y2 = y1;
			x3 = x1;
			y3 = bottom + (k + 1)*dy;
			x4 = x2;
			y4 = y3;

			a1 = phi(x1, y1);
			a2 = phi(x2, y2);
			a3 = phi(x3, y3);
			a4 = phi(x4, y4);

			//_________________________________________
			// set nbound and nu
			//_________________________________________
			if (j < m - 1)
			{
				id = j + k*(m - 1);
				if (a4 <  -1e-10)
				{
					if (k == n - 1) {
						ApertureIndex = ApertureIndex + 1;
						nu[id] = notApertureNum + ApertureIndex;
						nbound[0][ApertureIndex - 1] = notApertureNum + ApertureIndex;
						nbound[1][ApertureIndex - 1] = x4;
					}
					else {
						notApertureIndex = notApertureIndex + 1;
						nu[id] = notApertureIndex;
					}
				}
			}

			//处理L
			int l_sgn;
			setTriangleType(x1, x2, x3, y1, y2, y3, l_sgn);

			if (l_sgn == 10)
			{
				Triangle_Normal cur_tri;
				cur_tri.sgn = l_sgn;
				tri(x1, x2, x3, y1, y2, y3, cur_tri);
				//solveTri(cur_tri);
				L.Set_normal(j, k, cur_tri);
			}
			else if (l_sgn == -10)
			{
				Triangle_Normal cur_tri;
				cur_tri.sgn = l_sgn;
				tri(x1, x2, x3, y1, y2, y3, cur_tri);
				//solveTri(cur_tri);
				L.Set_normal(j, k, cur_tri);
			}
			else if (abs(l_sgn) < 10)
			{
				Triangle_All cur_tri;
				cur_tri.sgn = l_sgn;
				tri(x1, x2, x3, y1, y2, y3, cur_tri);
				solveTri(cur_tri);
				L.Set_all(j, k, cur_tri);
			}
			else if (abs(l_sgn) > 10)
			{
				Triangle_All cur_tri;
				cur_tri.sgn = l_sgn;
				tri(x1, x2, x3, y1, y2, y3, cur_tri);
				solveTri(cur_tri);
				L.Set_all(j, k, cur_tri);
			}


			//处理U
			int u_sgn;

			double tempx4 = x4;
			double tempx3 = x3;
			double tempx2 = x2;
			double tempy4 = y4;
			double tempy3 = y3;
			double tempy2 = y2;
			x1 = tempx4;
			x2 = tempx3;
			x3 = tempx2;
			y1 = tempy4;
			y2 = tempy3;
			y3 = tempy2;

			setTriangleType(x1, x2, x3, y1, y2, y3, u_sgn);

			if (u_sgn == 10)
			{
				Triangle_Normal cur_tri;
				cur_tri.sgn = u_sgn;
				tri(x1, x2, x3, y1, y2, y3, cur_tri);
				//solveTri(cur_tri);
				U.Set_normal(j, k, cur_tri);
			}
			else if (u_sgn == -10)
			{
				Triangle_Normal cur_tri;
				cur_tri.sgn = u_sgn;
				tri(x1, x2, x3, y1, y2, y3, cur_tri);
				//solveTri(cur_tri);
				U.Set_normal(j, k, cur_tri);
			}
			else if (abs(u_sgn) < 10)
			{
				Triangle_All cur_tri;
				cur_tri.sgn = u_sgn;
				tri(x1, x2, x3, y1, y2, y3, cur_tri);
				solveTri(cur_tri);
				//U.Set(j, k, cur_tri);
				U.Set_all(j, k, cur_tri);
			}
			else if (abs(u_sgn) > 10)
			{
				Triangle_All cur_tri;
				cur_tri.sgn = u_sgn;
				tri(x1, x2, x3, y1, y2, y3, cur_tri);
				solveTri(cur_tri);
				U.Set_all(j, k, cur_tri);
			}


			#pragma region MyRegion
			////_________________________________________
			////set triangle mesh
			////_________________________________________
			//Triangle l, u;
			//l.x = Vector3d(x1, x2, x3);
			//l.y = Vector3d(y1, y2, y3);
			//u.x = Vector3d(x4, x3, x2);
			//u.y = Vector3d(y4, y3, y2);

			////-------------------------------------
			////set value for l
			////-------------------------------------
			//if (abs(a1) < 1e-10 && abs(a2) < 1e-10 && abs(a3) < 1e-10)
			//{
			//	//三角形三点都在界面上
			//	double ac = phi((x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3);
			//	if (ac > 1e-10)
			//		l.sgn = 10;
			//	if (ac < -1e-10)
			//		l.sgn = -10;
			//}
			//if ((abs(a1) < 1e-10 && a2 > 1e-10 && a3 > 1e-10) ||
			//	(a1 > 1e-10 && abs(a2) < 1e-10 && a3 > 1e-10) ||
			//	(a1 > 1e-10 && a2 > 1e-10 && abs(a3) < 1e-10) ||
			//	(abs(a1) < 1e-10 && abs(a2) < 1e-10 && a3 > 1e-10) ||
			//	(abs(a1) < 1e-10 && a2 > 1e-10 && abs(a3) < 1e-10) ||
			//	(a1 > 1e-10 && abs(a2) < 1e-10 && abs(a3) < 1e-10) ||
			//	(a1 > 1e-10 && a2 > 1e-10 && a3 > 1e-10))
			//{
			//	// 对应三角形位于Ω + 中的情况：
			//	//一个点在界面上，两个点在Ω +
			//	//或两个点在界面上，一个点在Ω +
			//	//或三个点都在Ω +
			//	l.sgn = 10;
			//}
			//if ((abs(a1) < 1e-10 && a2 < -1e-10 && a3 < -1e-10) ||
			//	(a1 < -1e-10 && abs(a2) < 1e-10 && a3 < -1e-10) ||
			//	(a1 < -1e-10 && a2 < -1e-10 && abs(a3) < 1e-10) ||
			//	(abs(a1) < 1e-10 && abs(a2) < 1e-10 && a3 < -1e-10) ||
			//	(abs(a1) < 1e-10 && a2 < -1e-10 && abs(a3) < 1e-10) ||
			//	(a1 < -1e-10 && abs(a2) < 1e-10 && abs(a3) < 1e-10) ||
			//	(a1 < -1e-10 && a2 < -1e-10 && a3 < -1e-10))
			//{
			//	// 对应三角形位于Ω - 中的情况：
			//	//一个点在界面上，两个点在Ω -
			//	//或两个点在界面上，一个点在Ω -
			//	//或三个点都在Ω -
			//	l.sgn = -10;
			//}
			//if (a1 > 1e-10 && a2 < -1e-10 && a3 < -1e-10)
			//{
			//	//三角形越过界面
			//	// (x1, y1)位于Ω + (x2, y2)、(x3, y3)位于Ω -
			//	l.sgn = 1;
			//	l.r1 = findzero(x1, y1, x2, y2);
			//	l.r2 = findzero(x1, y1, x3, y3);
			//}
			//if (a1 < -1e-10 && a2>1e-10 && a3 > 1e-10)
			//{
			//	// 三角形越过界面
			//	//a1位于Ω - a2、a3位于Ω +
			//	l.sgn = -1;
			//	l.r1 = findzero(x1, y1, x2, y2);
			//	l.r2 = findzero(x1, y1, x3, y3);
			//}
			//if (a2 > 1e-10 && a1 < -1e-10 && a3 < -1e-10)
			//{
			//	//三角形越过界面
			//	//a2位于Ω + a1、a3位于Ω -
			//	l.sgn = 2;
			//	l.r1 = findzero(x2, y2, x3, y3);
			//	l.r2 = findzero(x2, y2, x1, y1);
			//}
			//if (a2<-1e-10 && a1>1e-10 && a3 > 1e-10)
			//{
			//	l.sgn = -2;
			//	l.r1 = findzero(x2, y2, x3, y3);
			//	l.r2 = findzero(x2, y2, x1, y1);
			//}
			//if (a3 > 1e-10 && a2 < -1e-10 && a1 < -1e-10)
			//{
			//	l.sgn = 3;
			//	l.r1 = findzero(x3, y3, x1, y1);
			//	l.r2 = findzero(x3, y3, x2, y2);
			//}
			//if (a3 < -1e-10 && a2>1e-10 && a1 > 1e-10)
			//{
			//	l.sgn = -3;
			//	l.r1 = findzero(x3, y3, x1, y1);
			//	l.r2 = findzero(x3, y3, x2, y2);
			//}
			////一个点在界面上，一个在Ω + ，一个在Ω -
			//if (abs(a2) < 1e-10 && a1 > 1e-10 && a3 < -1e-10)
			//{
			//	// a2在界面上，a1在Ω + ，a3在Ω -
			//	l.sgn = 11;
			//	l.A = Vector2d(x1, y1);  //与B按照逆时针相对排列
			//	l.B = Vector2d(x2, y2);  //位于Γ上的点
			//	l.C = Vector2d(x3, y3);  //与B按照逆时针相对排列
			//	l.v = Vector3i(1, 2, 3); //标识ABC三点分别对应的点
			//	l.r1 = 1;
			//	l.r2 = findzero(x1, y1, x3, y3);
			//}
			//if (abs(a2) < 1e-10 && a1 < -1e-10 && a3>1e-10)
			//{
			//	l.sgn = -11;
			//	l.A = Vector2d(x1, y1);
			//	l.B = Vector2d(x2, y2);
			//	l.C = Vector2d(x3, y3);
			//	l.v = Vector3i(1, 2, 3);
			//	l.r1 = 1;
			//	l.r2 = findzero(x1, y1, x3, y3);
			//}
			//if (abs(a3) < 1e-10 && a2 > 1e-10 && a1 < -1e-10)
			//{
			//	l.sgn = 12;
			//	l.A = Vector2d(x2, y2);
			//	l.B = Vector2d(x3, y3);
			//	l.C = Vector2d(x1, y1);
			//	l.v = Vector3i(2, 3, 1);
			//	l.r1 = 1;
			//	l.r2 = findzero(x2, y2, x1, y1);
			//}
			//if (abs(a3) < 1e-10 && a2 < -1e-10 && a1>1e-10)
			//{
			//	l.sgn = -12;
			//	l.A = Vector2d(x2, y2);
			//	l.B = Vector2d(x3, y3);
			//	l.C = Vector2d(x1, y1);
			//	l.v = Vector3i(2, 3, 1);
			//	l.r1 = 1;
			//	l.r2 = findzero(x2, y2, x1, y1);
			//}
			//if (abs(a1) < 1e-10 && a2 < -1e-10 && a3>1e-10)
			//{
			//	l.sgn = 13;
			//	l.A = Vector2d(x3, y3);
			//	l.B = Vector2d(x1, y1);
			//	l.C = Vector2d(x2, y2);
			//	l.v = Vector3i(3, 1, 2);
			//	l.r1 = 1;
			//	l.r2 = findzero(x3, y3, x2, y2);
			//}
			//if (abs(a1) < 1e-10 && a2 > 1e-10 && a3 < -1e-10)
			//{
			//	l.sgn = -13;
			//	l.A = Vector2d(x3, y3);
			//	l.B = Vector2d(x1, y1);
			//	l.C = Vector2d(x2, y2);
			//	l.v = Vector3i(3, 1, 2);
			//	l.r1 = 1;
			//	l.r2 = findzero(x3, y3, x2, y2);
			//}
			//

			////-------------------------------------
			////set value for u
			////-------------------------------------
			//if (abs(a4) < 1e-10 && abs(a3) < 1e-10 && abs(a2) < 1e-10)
			//{
			//	double ac = phi((x4 + x3 + x2) / 3, (y4 + y3 + y2) / 3);
			//	if (ac > 1e-10)
			//		u.sgn = 10;
			//	if (ac < -1e-10)
			//		u.sgn = -10;
			//}
			//if ((abs(a4) < 1e-10 && a3 > 1e-10 && a2 > 1e-10) ||
			//	(a4 > 1e-10 && abs(a3) < 1e-10 && a2 > 1e-10) ||
			//	(a4 > 1e-10 && a3 > 1e-10 && abs(a2) < 1e-10) ||
			//	(abs(a4) < 1e-10 && abs(a3) < 1e-10 && a2 > 1e-10) ||
			//	(abs(a4) < 1e-10 && a3 > 1e-10 && abs(a2) < 1e-10) ||
			//	(a4 > 1e-10 && abs(a3) < 1e-10 && abs(a2) < 1e-10) ||
			//	(a4 > 1e-10 && a3 > 1e-10 && a2 > 1e-10))
			//{
			//	u.sgn = 10;
			//}
			//if ((abs(a4) < 1e-10 && a3 < -1e-10 && a2 < -1e-10) ||
			//	(a4 < -1e-10 && abs(a3) < 1e-10 && a2 < -1e-10) ||
			//	(a4 < -1e-10 && a3 < -1e-10 && abs(a2) < 1e-10) ||
			//	(abs(a4) < 1e-10 && abs(a3) < 1e-10 && a2 < -1e-10) ||
			//	(abs(a4) < 1e-10 && a3 < -1e-10 && abs(a2) < 1e-10) ||
			//	(a4 < -1e-10 && abs(a3) < 1e-10 && abs(a2) < 1e-10) ||
			//	(a4 < -1e-10 && a3 < -1e-10 && a2 < -1e-10))
			//{
			//	u.sgn = -10;
			//}
			//if (a4 > 1e-10 && a3 < -1e-10 && a2 < -1e-10)
			//{
			//	u.sgn = 1;
			//	u.r1 = findzero(x4, y4, x3, y3);
			//	u.r2 = findzero(x4, y4, x2, y2);
			//}
			//if (a4 < -1e-10 && a3>1e-10 && a2 > 1e-10)
			//{
			//	u.sgn = -1;
			//	u.r1 = findzero(x4, y4, x3, y3);
			//	u.r2 = findzero(x4, y4, x2, y2);
			//}
			//if (a3 > 1e-10 && a2 < -1e-10 && a4 < -1e-10)
			//{
			//	u.sgn = 2;
			//	u.r1 = findzero(x3, y3, x2, y2);
			//	u.r2 = findzero(x3, y3, x4, y4);
			//}
			//if (a3 < -1e-10 && a2>1e-10 && a4 > 1e-10)
			//{
			//	u.sgn = -2;
			//	u.r1 = findzero(x3, y3, x2, y2);
			//	u.r2 = findzero(x3, y3, x4, y4);
			//}
			//if (a2 > 1e-10 && a4 < -1e-10 && a3 < -1e-10)
			//{
			//	u.sgn = 3;
			//	u.r1 = findzero(x2, y2, x4, y4);
			//	u.r2 = findzero(x2, y2, x3, y3);
			//}
			//if (a2 < -1e-10 && a4>1e-10 && a3 > 1e-10)
			//{
			//	u.sgn = -3;
			//	u.r1 = findzero(x2, y2, x4, y4);
			//	u.r2 = findzero(x2, y2, x3, y3);
			//}
			////一个点在界面上，一个在Ω + ，一个在Ω -
			//if (abs(a3) < 1e-10 && a4 > 1e-10 && a2 < -1e-10)
			//{
			//	u.sgn = 11;
			//	u.A = Vector2d(x4, y4);
			//	u.B = Vector2d(x3, y3);
			//	u.C = Vector2d(x2, y2);
			//	u.v = Vector3i(1, 2, 3);
			//	u.r1 = 1;
			//	u.r2 = findzero(x4, y4, x2, y2);
			//}
			//if (abs(a3) < 1e-10 && a4 < -1e-10 && a2>1e-10)
			//{
			//	u.sgn = -11;
			//	u.A = Vector2d(x4, y4);
			//	u.B = Vector2d(x3, y3);
			//	u.C = Vector2d(x2, y2);
			//	u.v = Vector3i(1, 2, 3);
			//	u.r1 = 1;
			//	u.r2 = findzero(x4, y4, x2, y2);
			//}
			//if (abs(a2) < 1e-10 && a3 > 1e-10 && a4 < -1e-10)
			//{
			//	u.sgn = 12;
			//	u.A = Vector2d(x3, y3);
			//	u.B = Vector2d(x2, y2);
			//	u.C = Vector2d(x4, y4);
			//	u.v = Vector3i(2, 3, 1);
			//	u.r1 = 1;
			//	u.r2 = findzero(x3, y3, x4, y4);
			//}
			//if (abs(a2) < 1e-10 && a3 < -1e-10 && a4>1e-10)
			//{
			//	u.sgn = -12;
			//	u.A = Vector2d(x3, y3);
			//	u.B = Vector2d(x2, y2);
			//	u.C = Vector2d(x4, y4);
			//	u.v = Vector3i(2, 3, 1);
			//	u.r1 = 1;
			//	u.r2 = findzero(x3, y3, x4, y4);
			//}
			//if (abs(a4) < 1e-10 && a3 < -1e-10 && a2>1e-10)
			//{
			//	u.sgn = 13;
			//	u.A = Vector2d(x2, y2);
			//	u.B = Vector2d(x4, y4);
			//	u.C = Vector2d(x3, y3);
			//	u.v = Vector3i(3, 1, 2);
			//	u.r1 = 1;
			//	u.r2 = findzero(x2, y2, x3, y3);
			//}
			//if (abs(a4) < 1e-10 && a3 > 1e-10 && a2 < -1e-10)
			//{
			//	u.sgn = -13;
			//	u.A = Vector2d(x2, y2);
			//	u.B = Vector2d(x4, y4);
			//	u.C = Vector2d(x3, y3);
			//	u.v = Vector3i(3, 1, 2);
			//	u.r1 = 1;
			//	u.r2 = findzero(x2, y2, x3, y3);
			//}


			////-------------------------------------
			////sort out
			////-------------------------------------
			//if (abs(l.sgn) < 10)
			//{
			//	int p1 = abs(l.sgn);
			//	int p2 = p1%3 + 1;
			//	int p3 = p2%3 + 1;
			//	l.x1 = l.x(p1-1); l.y1 = l.y(p1-1);
			//	l.x2 = l.x(p2-1); l.y2 = l.y(p2-1);
			//	l.x3 = l.x(p3-1); l.y3 = l.y(p3-1);
			//	l.x4 = l.x1 + (l.x2 - l.x1)*l.r1; l.y4 = l.y1 + (l.y2 - l.y1)*l.r1;
			//	l.x5 = l.x1 + (l.x3 - l.x1)*l.r2; l.y5 = l.y1 + (l.y3 - l.y1)*l.r2;
			//}
			//if (abs(u.sgn) < 10)
			//{
			//	int p1 = abs(u.sgn);
			//	int p2 = p1%3 + 1;
			//	int p3 = p2%3 + 1;
			//	u.x1 = u.x(p1-1); u.y1 = u.y(p1-1);
			//	u.x2 = u.x(p2-1); u.y2 = u.y(p2-1);
			//	u.x3 = u.x(p3-1); u.y3 = u.y(p3-1);
			//	u.x4 = u.x1 + (u.x2 - u.x1)*u.r1; u.y4 = u.y1 + (u.y2 - u.y1)*u.r1;
			//	u.x5 = u.x1 + (u.x3 - u.x1)*u.r2; u.y5 = u.y1 + (u.y3 - u.y1)*u.r2;
			//}
			////-------------------------------------
			//if (abs(l.sgn) > 10)
			//{
			//	l.p5 = l.A + (l.C - l.A)*l.r2;
			//	l.x1 = l.A(0); l.y1 = l.A(1);
			//	l.x2 = l.B(0); l.y2 = l.B(1);
			//	l.x3 = l.C(0); l.y3 = l.C(1);
			//}
			//if (abs(u.sgn) > 10)
			//{
			//	u.p5 = u.A + (u.C - u.A)*u.r2;
			//	u.x1 = u.A(0); u.y1 = u.A(1);
			//	u.x2 = u.B(0); u.y2 = u.B(1);
			//	u.x3 = u.C(0); u.y3 = u.C(1);
			//}
			//

			//// -------------------------------------
			//U.Set(j, k, u);
			//L.Set(j, k, l);
			#pragma endregion			
		}
	}
}


VectorXcd SingalCavity::compute_g(MatrixXcd &G, vector<vector<double>> &nbound)
{
	double alpha = this->k0 * sin(this->theta);
	double beta = this->k0 * cos(this->theta);

	int size = G.rows();
	VectorXcd g(size);
	//g = -2 * 1i * beta * exp(1i * alpha * nbound(2, :)') ; 
	complex<double> i(0, 1);
	for (int j = 0; j < size; j++)
	{
		complex<double> temp(-2, 0);
		temp *= i;
		temp *= complex<double>(beta, 0);
		complex<double> temp2(alpha, 0);
		temp2 *= i;
		temp2 *= complex<double>(nbound[1][j], 0);
		g(j) = temp*pow(M_E, temp2);
	}
	return g;
}


vector<vector<gridCell>> SingalCavity::setGrid(TriangleMesh &U, TriangleMesh &L, vector<vector<double>> &nbound, vector<int> &nu)
{
	int m = this->meshWidth;
	int n = this->meshHeight;
	double dx = this->stepX;
	double dy = this->stepY;

	double left = this->virtualBorderLeft;
	double bottom = this->virtualBorderBottom;

	vector<vector<gridCell>> gridCellVector;
	gridCellVector.resize(m - 1);
	for (int i = 0; i < m - 1; i++)
		gridCellVector[i].resize(n);
	
	for (int j = 0; j < m - 1; j++)
	{
		for (int l = 0; l < n; l++)
		{
			// 3      4
			// ------- （j, k）
			// | \ U |
			// |  \  |
			// | L \ |
			// |    \|
			// -------
			// 1      2

			int id = j + l*(m - 1);
			gridCell grid;
			if (nu[id] > 0) 
			{
				// 只关心Ω - （腔体内）的点

				vector <complex<double>> t1;
				vector <complex<double>> t2;
				vector <complex<double>> t3;
				vector <complex<double>> t4;
				vector <complex<double>> t5;
				vector <complex<double>> t6;
				
				vector <complex<double>> t1m;
				vector <complex<double>> t2m;
				vector <complex<double>> t3m;
				vector <complex<double>> t4m;
				vector <complex<double>> t5m;
				vector <complex<double>> t6m;

				complex<double> zero(0, 0);
				int nboundSize = nbound[0].size();
				if (l == n - 1) // 口径面
				{
					initializeComplexVector(t1, zero, 4);
					initializeComplexVector(t1m, zero, nboundSize);
					initializeComplexVector(t2, zero, 4);
					initializeComplexVector(t2m, zero, nboundSize);
					initializeComplexVector(t3, zero, 4);
					initializeComplexVector(t3m, zero, nboundSize);
					if (U.Get_sign(j, l))
						weak(U.Get_normal(j, l), Vector3d(1, 0, 0), 1, nbound, t4, t4m);
					else
						weak(U.Get_all(j, l), Vector3d(1, 0, 0), 1, nbound, t4, t4m);

					if (L.Get_sign(j + 1, l))
						weak(L.Get_normal(j + 1, l), Vector3d(0, 0, 1), 0, nbound, t5, t5m);
					else
						weak(L.Get_all(j + 1, l), Vector3d(0, 0, 1), 0, nbound, t5, t5m);

					if (U.Get_sign(j + 1, l))
						weak(U.Get_normal(j + 1, l), Vector3d(0, 1, 0), 1, nbound, t6, t6m);
					else
						weak(U.Get_all(j + 1, l), Vector3d(0, 1, 0), 1, nbound, t6, t6m);
				}
				else if (l == n - 2) // 口径面下面一排
				{
					//[1, 0, 0][0, 1, 0][0, 0, 1]用于标识当前点是三角形的哪个点
					if (L.Get_sign(j + 1, l + 1))
						weak(L.Get_normal(j + 1, l + 1), Vector3d(1, 0, 0), 0, nbound, t1, t1m);
					else
						weak(L.Get_all(j + 1, l + 1), Vector3d(1, 0, 0), 0, nbound, t1, t1m);

					if (U.Get_sign(j, l + 1))
						weak(U.Get_normal(j, l + 1), Vector3d(0, 0, 1), 1, nbound, t2, t2m);
					else
						weak(U.Get_all(j, l + 1), Vector3d(0, 0, 1), 1, nbound, t2, t2m);

					if (L.Get_sign(j, l + 1))
						weak(L.Get_normal(j, l + 1), Vector3d(0, 1, 0), 0, nbound, t3, t3m);
					else
						weak(L.Get_all(j, l + 1), Vector3d(0, 1, 0), 0, nbound, t3, t3m);

					if (U.Get_sign(j, l))
						weak(U.Get_normal(j, l), Vector3d(1, 0, 0), 0, nbound, t4, t4m);
					else
						weak(U.Get_all(j, l), Vector3d(1, 0, 0), 0, nbound, t4, t4m);

					if (L.Get_sign(j + 1, l))
						weak(L.Get_normal(j + 1, l), Vector3d(0, 0, 1), 0, nbound, t5, t5m);
					else
						weak(L.Get_all(j + 1, l), Vector3d(0, 0, 1), 0, nbound, t5, t5m);

					if (U.Get_sign(j + 1, l))
						weak(U.Get_normal(j + 1, l), Vector3d(0, 1, 0), 0, nbound, t6, t6m);
					else
						weak(U.Get_all(j + 1, l), Vector3d(0, 1, 0), 0, nbound, t6, t6m);
				}
				else
				{
					if (L.Get_sign(j + 1, l + 1))
						weak(L.Get_normal(j + 1, l + 1), Vector3d(1, 0, 0), 0, nbound, t1, t1m);
					else
						weak(L.Get_all(j + 1, l + 1), Vector3d(1, 0, 0), 0, nbound, t1, t1m);

					if (U.Get_sign(j, l + 1))
						weak(U.Get_normal(j, l + 1), Vector3d(0, 0, 1), 0, nbound, t2, t2m);
					else
						weak(U.Get_all(j, l + 1), Vector3d(0, 0, 1), 0, nbound, t2, t2m);

					if (L.Get_sign(j, l + 1))
						weak(L.Get_normal(j, l + 1), Vector3d(0, 1, 0), 0, nbound, t3, t3m);
					else
						weak(L.Get_all(j, l + 1), Vector3d(0, 1, 0), 0, nbound, t3, t3m);

					if (U.Get_sign(j, l))
						weak(U.Get_normal(j, l), Vector3d(1, 0, 0), 0, nbound, t4, t4m);
					else
						weak(U.Get_all(j, l), Vector3d(1, 0, 0), 0, nbound, t4, t4m);

					if (L.Get_sign(j + 1, l))
						weak(L.Get_normal(j + 1, l), Vector3d(0, 0, 1), 0, nbound, t5, t5m);
					else
						weak(L.Get_all(j + 1, l), Vector3d(0, 0, 1), 0, nbound, t5, t5m);

					if (U.Get_sign(j + 1, l))
						weak(U.Get_normal(j + 1, l), Vector3d(0, 1, 0), 0, nbound, t6, t6m);
					else
						weak(U.Get_all(j + 1, l), Vector3d(0, 1, 0), 0, nbound, t6, t6m);
				}
				grid._c = t1[0] + t2[2] + t3[1] + t4[0] + t5[2] + t6[1];
				grid._w = t3[0] + t4[1];
				grid._s = t4[2] + t5[0];
				grid._n = t1[2] + t2[0];
				grid._e = t1[1] + t6[0];
				grid._nw = t2[1] + t3[2];
				grid._se = t5[1] + t6[2];
				grid._const = t1[3] + t2[3] + t3[3] + t4[3] + t5[3] + t6[3];

				//grid.tm = t1m + t2m + t3m + t4m + t5m + t6m;
				vector <complex<double>> temp_tm;
				initializeComplexVector(temp_tm, zero, nboundSize);
				complexVectorAdd(temp_tm, t1m);
				complexVectorAdd(temp_tm, t2m);
				complexVectorAdd(temp_tm, t3m);
				complexVectorAdd(temp_tm, t4m);
				complexVectorAdd(temp_tm, t5m);
				complexVectorAdd(temp_tm, t6m);
				grid.tm = temp_tm;

				gridCellVector[j][l] = grid;
			}
			else
			{
				//out{ j,l } = [];
				gridCellVector[j][l] = grid;
			}
		}
	}
	return gridCellVector;
}


void SingalCavity::weak(Triangle_Normal &T, Vector3d v, int topsign, vector<vector<double>> &nbound, vector <complex<double>> &out, vector <complex<double>> &out5m)
{
	int nboundSize = nbound[0].size();
	complex<double> zero(0, 0);

	out.resize(4);
	initializeComplexVector(out, zero, 4);
	out5m.resize(nboundSize);
	initializeComplexVector(out5m, zero, nboundSize);

	//phi = para.phi;
	//beta = para.beta;
	//u = para.g;
	double top = this->apertureY;
	if (T.sgn == -10)
	{
		Vector3cd w1 = weak1(T, v) + weak3(T, v);
		complex<double> out4(0, 0);
		//若三角形有顶点位于不规则边界上，计算不规则边界的积分
		if (abs(phi(T.x(0), T.y(0))) < 1e-10)
		{
			out4 = out4 + w1(0)*u(T.x(0), T.y(0));
			w1(0) = 0;
		}
		if (abs(phi(T.x(1), T.y(1))) < 1e-10)
		{
			out4 = out4 + w1(1)*u(T.x(1), T.y(1));
			w1(1) = 0;
		}
		if (abs(phi(T.x(2), T.y(2))) < 1e-10)
		{
			out4 = out4 + w1(2)*u(T.x(2), T.y(2));
			w1(2) = 0;
		}
		complex<double> out5(0, 0);
		if (topsign == 1)
		{
			if (abs(T.y(0) - top) < 1e-10&&abs(T.y(1) - top) < 1e-10)
				weak5(T.x(0), T.y(0), T.x(1), T.y(1), v(0), v(1), nbound, out5m, out5);
			else if (abs(T.y(1) - top) < 1e-10&&abs(T.y(2) - top) < 1e-10)
				weak5(T.x(1), T.y(1), T.x(2), T.y(2), v(1), v(2), nbound, out5m, out5);
			else if (abs(T.y(2) - top) < 1e-10&&abs(T.y(0) - top) < 1e-10)
				weak5(T.x(2), T.y(2), T.x(0), T.y(0), v(2), v(0), nbound, out5m, out5);
		}
		out[0] = w1(0);
		out[1] = w1(1);
		out[2] = w1(2);
		out[3] = out4 - out5;
		return;
	}
}


void SingalCavity::weak(Triangle_All &T, Vector3d v, int topsign, vector<vector<double>> &nbound, vector <complex<double>> &out, vector <complex<double>> &out5m)
{
	int nboundSize = nbound[0].size();
	complex<double> zero(0, 0);

	out.resize(4);
	initializeComplexVector(out, zero, 4);
	out5m.resize(nboundSize);
	initializeComplexVector(out5m, zero, nboundSize);

	//phi = para.phi;
	//beta = para.beta;
	//u = para.g;
	double top = this->apertureY;
	
	if (abs(T.sgn) < 10)
	{
		int s;
		if (T.sgn > 0)
			s = T.sgn;
		else
			s = -T.sgn;

		double r1 = T.r1;
		double r2 = T.r2;

		int j = s;
		int k = j % 3 + 1;
		int l = k % 3 + 1;

		int vp1 = v(j-1);
		int vp2 = v(k-1);
		int vp3 = v(l-1);
		double vp4 = vp1 + (vp2 - vp1)*r1;
		double vp5 = vp1 + (vp3 - vp1)*r2;

		Vector3d v1(vp1, vp4, vp5);
		Vector3d v2(vp2, vp3, vp4);
		Vector3d v3(vp3, vp5, vp4);

		double x1 = T.x1;
		double y1 = T.y1;
		double x2 = T.x2;
		double y2 = T.y2;
		double x3 = T.x3;
		double y3 = T.y3;
		double x4 = T.x4;
		double y4 = T.y4;
		double x5 = T.x5;
		double y5 = T.y5;
		double x6 = (x4 + x5) / 2;
		double y6 = (y4 + y5) / 2;
		double x7 = (x2 + x3) / 2;
		double y7 = (y2 + y3) / 2;

		Vector2d p4(x4, y4);
		Vector2d p5(x5, y5);
		Vector2d p6(x6, y6);
		Vector2d p7(x7, y7);

		Triangle_All T1, T2, T3;
		T1.x = Vector3d(x1, x4, x5);
		T1.y = Vector3d(y1, y4, y5);
		T2.x = Vector3d(x2, x3, x4);
		T2.y = Vector3d(y2, y3, y4);
		T3.x = Vector3d(x3, x5, x4);
		T3.y = Vector3d(y3, y5, y4);

		complex<double> u4 = u(x4, y4);
		complex<double> u5 = u(x5, y5);

		complex<double> c1, c2, c3, c4;
		if (T.sgn > 0)
		{
			Vector3cd out2 = weak1(T2, v2) + weak3(T2, v2);
			Vector3cd out3 = weak1(T3, v3) + weak3(T3, v3);
			Vector3cd out4 = weak4(p5, p4, p7, vp5, vp4, p6);
			
			complex<double> out5(0, 0);
			if (topsign == 1)
			{
				if (abs(y1 - top) < 1e-10&&abs(y2 - top) < 1e-10)
					weak5(x4, y4, x2, y2, vp4, vp2, nbound, out5m, out5);
				else if (abs(y2 - top) < 1e-10&&abs(y3 - top) < 1e-10)
					weak5(x2, y2, x3, y3, vp2, vp3, nbound, out5m, out5);
				else if (abs(y3 - top) < 1e-10&&abs(y1 - top) < 1e-10)
					weak5(x3, y3, x5, y5, vp3, vp5, nbound, out5m, out5);
			}
			c1 = 0;
			c2 = out2(0) + out4(2)*0.5;
			c3 = out2(1) + out3(0) + out4(2)*0.5;
			c4 = out2(2)*u4 + out3(1)*u5 + out3(2)*u4 + out4(0)*u5 + out4(1)*u4 - out5;
		}
		else
		{
			Vector3cd out1 = weak1(T1, v1) + weak3(T1, v1);
			Vector3cd out4 = weak4(p4, p5, Vector2d(x1, y1), vp4, vp5, p6);
			complex<double> out5(0, 0);
			if (topsign == 1)
			{
				if (abs(y1 - top) < 1e-10&&abs(y2 - top) < 1e-10)
					weak5(x1, y1, x4, y4, vp1, vp4, nbound, out5m, out5);
				else if (abs(y3 - top) < 1e-10&&abs(y1 - top) < 1e-10)
					weak5(x5, y5, x1, y1, vp5, vp1, nbound, out5m, out5);
			}
			c1 = out1(0) + out4(2);
			c2 = 0;
			c3 = 0;
			c4 = out1(1)*u4 + out1(2)*u5 + out4(0)*u4 + out4(1)*u5 - out5;
		}

		if (s == 1)
		{
			out[0] = c1;
			out[1] = c2;
			out[2] = c3;
			out[3] = c4;
		}
		if (s == 2)
		{
			out[0] = c3;
			out[1] = c1;
			out[2] = c2;
			out[3] = c4;
		}
		if (s == 3)
		{
			out[0] = c2;
			out[1] = c3;
			out[2] = c1;
			out[3] = c4;
		}
		return;
	}

	if (abs(T.sgn) > 10)
	{
		double vp1 = v(T.v(0)-1);
		double vp2 = v(T.v(1)-1);
		double vp3 = v(T.v(2)-1);
		double vp5 = vp1 + (vp3 - vp1)*T.r2;

		Vector3d v1(vp1, vp2, vp5);
		Vector3d v2(vp2, vp3, vp5);

		double x1 = T.x1;
		double y1 = T.y1;
		double x2 = T.x2;
		double y2 = T.y2;
		double x3 = T.x3;
		double y3 = T.y3;
		double x5 = T.p5(0);
		double y5 = T.p5(1);
		double x6 = (x2 + x5) / 2; 
		double y6 = (y2 + y5) / 2;

		Vector2d p2(x2, y2);
		Vector2d p5(x5, y5);
		Vector2d p6(x6, y6);

		Triangle_All T1, T2;
		T1.x = Vector3d(x1, x2, x5);
		T1.y = Vector3d(y1, y2, y5);
		T2.x = Vector3d(x2, x3, x5);
		T2.y = Vector3d(y2, y3, y5);

		complex<double> u2 = u(x2, y2);
		complex<double> u5 = u(x5, y5);

		complex<double> c1, c2, c3, c4;
		if (T.sgn > 0)
		{
			Vector3cd out2 = weak1(T2, v2) + weak3(T2, v2);
			Vector3cd out4 = weak4(p5, p2, Vector2d(x3, y3), vp5, vp2, p6);
			complex<double> out5(0, 0);
			if (topsign == 1)
			{
				if (abs(y2 - top) < 1e-10&&abs(y3 - top) < 1e-10)
					weak5(x2, y2, x3, y3, vp2, vp3, nbound, out5m, out5);
				else if( abs(y3 - top) < 1e-10&&abs(y1 - top) < 1e-10)
					weak5(x3, y3, x5, y5, vp3, vp5, nbound, out5m, out5);
			}
			c1 = 0;
			c2 = 0;
			c3 = out2(1) + out4(2);
			c4 = out2(0)*u2 + out2(2)*u5 + out4(0)*u5 + out4(1)*u2 - out5;
		}
		if (T.sgn < 0)
		{
			Vector3cd out1 = weak1(T1,v1) + weak3(T1, v1);
			Vector3cd out4 = weak4(p2, p5, Vector2d(x1, y1), vp2, vp5, p6);
			complex<double> out5(0, 0);
			if (topsign == 1)
			{
				if (abs(y1 - top) < 1e-10&&abs(y2 - top) < 1e-10)
					weak5(x1, y1, x2, y2, vp1, vp2, nbound, out5m, out5);
				else if (abs(y3 - top) < 1e-10&&abs(y1 - top) < 1e-10)
					weak5(x5, y5, x1, y1, vp5, vp1, nbound, out5m, out5);
			}
			c1 = out1(0) + out4(2);
			c2 = 0;
			c3 = 0;
			c4 = out1(1)*u2 + out1(2)*u5 + out4(0)*u2 + out4(1)*u5 - out5;
		}

		out[T.v(0)-1] = c1;
		out[T.v(1)-1] = c2;
		out[T.v(2)-1] = c3;
		out[3] = c4;

		return;
	}
}


Vector3cd SingalCavity::weak1(Triangle_Normal &T, Vector3d v)
{
	// 求拉普拉斯u的积分

	double s = calculateTriangleArea(T.x, T.y);//求三角形面积
	Matrix3d M;
	M << T.x(0), T.y(0), 1, T.x(1), T.y(1), 1, T.x(2), T.y(2), 1;

	//M = [T.x', T.y', [1; 1; 1]];
	Matrix3d invM = M.inverse(); //求M的逆

								 //M = M([1, 2], :);
	MatrixXd MM(2, 3);
	MM = invM.block(0, 0, 2, 3);

	Vector2d vdiff = MM*v;

	double x4 = (T.x(1) + T.x(2)) / 2;
	double x5 = (T.x(2) + T.x(0)) / 2;
	double x6 = (T.x(0) + T.x(1)) / 2;
	double y4 = (T.y(1) + T.y(2)) / 2;
	double y5 = (T.y(2) + T.y(0)) / 2;
	double y6 = (T.y(0) + T.y(1)) / 2;

	Matrix2d b4 = beta(x4, y4);
	Matrix2d b5 = beta(x5, y5);
	Matrix2d b6 = beta(x6, y6);

	Vector3cd out = vdiff.transpose() *(b4 + b5 + b6)*MM*s / 3;

	return out;
}


Vector3cd SingalCavity::weak1(Triangle_All &T, Vector3d v)
{
	// 求拉普拉斯u的积分
	
	double s = calculateTriangleArea(T.x, T.y);//求三角形面积
	Matrix3d M;
	M << T.x(0), T.y(0), 1, T.x(1), T.y(1), 1, T.x(2), T.y(2), 1 ;

	//M = [T.x', T.y', [1; 1; 1]];
	Matrix3d invM = M.inverse(); //求M的逆
		
	//M = M([1, 2], :);
	MatrixXd MM(2, 3);
	MM = invM.block(0, 0, 2, 3);

	Vector2d vdiff = MM*v;

	double x4 = (T.x(1) + T.x(2)) / 2;
	double x5 = (T.x(2) + T.x(0)) / 2;
	double x6 = (T.x(0) + T.x(1)) / 2;
	double y4 = (T.y(1) + T.y(2)) / 2;
	double y5 = (T.y(2) + T.y(0)) / 2;
	double y6 = (T.y(0) + T.y(1)) / 2;

	Matrix2d b4 = beta(x4, y4);
	Matrix2d b5 = beta(x5, y5);
	Matrix2d b6 = beta(x6, y6);

	Vector3cd out = vdiff.transpose() *(b4+b5+b6)*MM*s/3;

	return out;
}


Vector3cd SingalCavity::weak3(Triangle_Normal &T, Vector3d v)
{
	Vector3cd out;

	double s = calculateTriangleArea(T.x, T.y);//求三角形面积

	double x4 = (T.x(1) + T.x(2)) / 2;
	double x5 = (T.x(2) + T.x(0)) / 2;
	double x6 = (T.x(0) + T.x(1)) / 2;
	double y4 = (T.y(1) + T.y(2)) / 2;
	double y5 = (T.y(2) + T.y(0)) / 2;
	double y6 = (T.y(0) + T.y(1)) / 2;

	complex<double> q4 = q(x4, y4);
	complex<double> q5 = q(x5, y5);
	complex<double> q6 = q(x6, y6);

	double v4 = (v(1) + v(2)) / 2;
	double v5 = (v(2) + v(0)) / 2;
	double v6 = (v(0) + v(1)) / 2;

	complex<double> u1 = (q5*v5 + q6*v6)*s *(1.0 / 6);
	complex<double> u2 = (q6*v6 + q4*v4)*s *(1.0 / 6);
	complex<double> u3 = (q4*v4 + q5*v5)*s *(1.0 / 6);

	out = Vector3cd(u1, u2, u3);

	return out;
}


Vector3cd SingalCavity::weak3(Triangle_All &T, Vector3d v)
{
	Vector3cd out;

	double s = calculateTriangleArea(T.x, T.y);//求三角形面积

	double x4 = (T.x(1) + T.x(2)) / 2;
	double x5 = (T.x(2) + T.x(0)) / 2;
	double x6 = (T.x(0) + T.x(1)) / 2;
	double y4 = (T.y(1) + T.y(2)) / 2;
	double y5 = (T.y(2) + T.y(0)) / 2;
	double y6 = (T.y(0) + T.y(1)) / 2;

	complex<double> q4 = q(x4, y4);
	complex<double> q5 = q(x5, y5);
	complex<double> q6 = q(x6, y6);

	double v4 = (v(1) + v(2)) / 2;
	double v5 = (v(2) + v(0)) / 2;
	double v6 = (v(0) + v(1)) / 2;

	complex<double> u1 = (q5*v5 + q6*v6)*s *(1.0 / 6);
	complex<double> u2 = (q6*v6 + q4*v4)*s *(1.0 / 6);
	complex<double> u3 = (q4*v4 + q5*v5)*s *(1.0 / 6);

	out = Vector3cd(u1, u2, u3);

	return out;
}


Vector3cd SingalCavity::weak4(Vector2d &p1, Vector2d &p2, Vector2d &p3, double v1, double v2, Vector2d &p6)
{
	Vector3cd out;

	Matrix3d M;
	M << p1(0), p1(1), 1, p2(0), p2(1), 1, p3(0), p3(1), 1;

	Matrix3d invM = M.inverse(); //求M的逆

	MatrixXd MM(2, 3);
	MM = invM.block(0, 0, 2, 3);

	double x1 = p1(0);
	double y1 = p1(1);
	double x2 = p2(0);
	double y2 = p2(1);
	double s = (p2 - p1).norm();//s = norm(p2 - p1).;
	double n1 = y2 - y1;
	double n2 = x1 - x2;
	Vector2d nb(n1 / s, n2 / s); 
	Matrix2d b6 = beta(p6(0), p6(1));

	out = -nb.transpose()*b6*MM*s*(v1 + v2) / 2;

	return out;
}

void SingalCavity::weak5(double x1, double y1, double x2, double y2, double v1, double v2, vector<vector<double>> &nbound, vector<complex<double>> &out5m, complex<double> &out5)
{
	double s = sqrt(pow(x2 - x1,2)+ pow(y2 - y1, 2));
	int nboundSize = nbound[0].size();

	if (abs(nbound[1][0] - max(x1, x2)) < 1e-10)
	{
		MatrixXcd G_row(1, nboundSize);
		G_row = G_aperture.block(0, 0, 1, nboundSize);
		VectorXcd temp5m = -G_row.transpose()*(v1 + v2)*s / 2;
		out5m = eigenVector2stdVector(temp5m);
		out5 = g_aperture(0)*(v1 + v2)*s * 0.5;
	}
	else if (abs(nbound[1][nboundSize-1] - min(x1, x2)) < 1e-10)
	{
		MatrixXcd G_row(1, nboundSize);
		G_row = G_aperture.block(nboundSize-1, 0, 1, nboundSize);
		VectorXcd temp5m = -G_row.transpose()*(v1 + v2)*s / 2;
		out5m = eigenVector2stdVector(temp5m);
		out5 = g_aperture(nboundSize-1)*(v1 + v2)*s * 0.5;
	}
	else
	{
		int sq = findXinNbound(nbound, min(x1, x2));
		MatrixXcd G_row1= G_aperture.block(sq, 0, 1, nboundSize);
		MatrixXcd G_row2 = G_aperture.block(sq+1, 0, 1, nboundSize);
		VectorXcd temp5m = -(G_row1.transpose() +G_row2.transpose())*(v1 + v2)*s / 4;
		out5m = eigenVector2stdVector(temp5m);
		out5 = (g_aperture(sq) + g_aperture(sq + 1)) * (v1 + v2)*s *0.25;
	}
}

VectorXcd SingalCavity::setRightHand(TriangleMesh &U, TriangleMesh &L, vector<int> &nu)
{
	int size = nu.size();
	complex<double> zero(0, 0);
	VectorXcd right = VectorXcd::Constant(size,zero);


	return right;
}


SparseMatrix<complex<double>> SingalCavity::setA(double nn, vector<int> &nu, vector<vector<double>> &nbound, vector<vector<gridCell>> &grid)
{
	int m = this->meshWidth;
	int n = this->meshHeight;

	double dx = this->stepX;
	double dy = this->stepY;

	double left = this->virtualBorderLeft;
	double bottom = this->virtualBorderBottom;

	int apertureNumber = nbound[0].size();  

	double* A, *pA;
	A = (double*)malloc(sizeof(double)*(nn * 7 + apertureNumber*apertureNumber) * 4);
	pA = A;
	if (A == NULL) {
		//printf("allocate memory for data failed.");
		//return(A);
		throw "allocate memory for data failed.";
	}

	int list = 0;
	int id;
	double x, y;
	int sq;
	double cur_nu;
	for (int j = 0; j < m - 1; j++)
	{
		for (int l = 0; l < n; l++)
		{
			// j=j+1;l=l+1;
			//x = left + (j + 1)*dx;
			//y = bottom + (l + 1)*dy;
			//以当前点（j,l）为中心，分别将其在周围6点处的积分值组装到对应位置
			id = j + l*(m - 1);
			cur_nu = nu[id];
			if (cur_nu > 0)
			{
				gridCell curGrid = grid[j][l];

				if (l == n - 1)
					sq = findXIndexinNbound(nbound, cur_nu);

				// west
				if (nu[id - 1] > 0)
				{
					if (l == n - 1)
					{
						//sq = findApertureIndex(nbound, apertureNumber, cur_nu);
						curGrid._w += curGrid.tm[sq - 1];
					}
					*(pA++) = cur_nu;
					*(pA++) = nu[id - 1];
					*(pA++) = curGrid._w.real();
					*(pA++) = curGrid._w.imag();
				}
				// north - west
				if (l != n - 1 && nu[id + m - 2] > 0)
				{
					*(pA++) = cur_nu;
					*(pA++) = nu[id + m - 2];
					*(pA++) = curGrid._nw.real();
					*(pA++) = curGrid._nw.imag();
				}
				// north
				if (l != n - 1 && nu[id + m - 1] > 0)
				{
					*(pA++) = cur_nu;
					*(pA++) = nu[id + m - 1];
					*(pA++) = curGrid._n.real();
					*(pA++) = curGrid._n.imag();
				}
				// east
				if (nu[id + 1] > 0)
				{
					if (l == n - 1)
					{
						//sq = findApertureIndex(nbound, apertureNumber, cur_nu);
						curGrid._e += curGrid.tm[sq + 1];
					}
					*(pA++) = cur_nu;
					*(pA++) = nu[id + 1];
					*(pA++) = curGrid._e.real();
					*(pA++) = curGrid._e.imag();
				}
				// south
				if (nu[id - m + 1] > 0)
				{
					*(pA++) = cur_nu;
					*(pA++) = nu[id - m + 1];
					*(pA++) = curGrid._s.real();
					*(pA++) = curGrid._s.imag();
				}
				// south - east
				if (nu[id - m + 2] > 0)
				{
					*(pA++) = cur_nu;
					*(pA++) = nu[id - m + 2];
					*(pA++) = curGrid._se.real();
					*(pA++) = curGrid._se.imag();
				}
				// 将 l = n-1 时，在口径面剩余点处产生的积分值进行组装
				if (l == n - 1)
				{
					//sq = findApertureIndex(nbound, apertureNumber, cur_nu);
					curGrid._c += curGrid.tm[sq];
					for (int apeIndex = 0; apeIndex < apertureNumber; apeIndex++)
					{
						if (apeIndex<sq - 1 || apeIndex>sq + 1)//这个条件确保不会再次考虑已经组装的部分
						{
							*(pA++) = cur_nu;
							*(pA++) = nbound[0][apeIndex];
							*(pA++) = curGrid.tm[apeIndex].real();
							*(pA++) = curGrid.tm[apeIndex].imag();
						}
					}
				}
				*(pA++) = cur_nu;
				*(pA++) = cur_nu;
				*(pA++) = curGrid._c.real();
				*(pA++) = curGrid._c.imag();
			}

		}
	}

	//截去多余内存区域
	int ArowNum = (pA - A)/4;
	double *A_result = (double*)malloc(ArowNum * 4 * sizeof(double));
	//memcpy(A_result, A, (pA - A)- sizeof(double));
	memcpy(A_result, A, ArowNum * 4 * sizeof(double));

	//将三元组形式数据转换为系数矩阵形式
	SparseMatrix<complex<double>> matrixA = buildSparseMatrixA(A_result, nn, ArowNum);

	return matrixA;
}


mxArray* SingalCavity::setA_mx(double nn, vector<int> &nu, vector<vector<double>> &nbound, vector<vector<gridCell>> &grid)
{
	int m = this->meshWidth;
	int n = this->meshHeight;

	double dx = this->stepX;
	double dy = this->stepY;

	double left = this->virtualBorderLeft;
	double bottom = this->virtualBorderBottom;

	int apertureNumber = nbound[0].size();

	double* A, *pA;
	A = (double*)malloc(sizeof(double)*(nn * 7 + apertureNumber*apertureNumber) * 4);
	pA = A;
	if (A == NULL) {
		//printf("allocate memory for data failed.");
		//return(A);
		throw "allocate memory for data failed.";
	}

	int list = 0;
	int id;
	double x, y;
	int sq;
	double cur_nu;
	for (int j = 0; j < m - 1; j++)
	{
		for (int l = 0; l < n; l++)
		{
			// j=j+1;l=l+1;
			//x = left + (j + 1)*dx;
			//y = bottom + (l + 1)*dy;
			//以当前点（j,l）为中心，分别将其在周围6点处的积分值组装到对应位置
			id = j + l*(m - 1);
			cur_nu = nu[id];
			if (cur_nu > 0)
			{
				gridCell curGrid = grid[j][l];

				if (l == n - 1)
					sq = findXIndexinNbound(nbound, cur_nu);

				// west
				if (nu[id - 1] > 0)
				{
					if (l == n - 1)
					{
						//sq = findApertureIndex(nbound, apertureNumber, cur_nu);
						curGrid._w += curGrid.tm[sq - 1];
					}
					*(pA++) = cur_nu;
					*(pA++) = nu[id - 1];
					*(pA++) = curGrid._w.real();
					*(pA++) = curGrid._w.imag();
				}
				// north - west
				if (l != n - 1 && nu[id + m - 2] > 0)
				{
					*(pA++) = cur_nu;
					*(pA++) = nu[id + m - 2];
					*(pA++) = curGrid._nw.real();
					*(pA++) = curGrid._nw.imag();
				}
				// north
				if (l != n - 1 && nu[id + m - 1] > 0)
				{
					*(pA++) = cur_nu;
					*(pA++) = nu[id + m - 1];
					*(pA++) = curGrid._n.real();
					*(pA++) = curGrid._n.imag();
				}
				// east
				if (nu[id + 1] > 0)
				{
					if (l == n - 1)
					{
						//sq = findApertureIndex(nbound, apertureNumber, cur_nu);
						curGrid._e += curGrid.tm[sq + 1];
					}
					*(pA++) = cur_nu;
					*(pA++) = nu[id + 1];
					*(pA++) = curGrid._e.real();
					*(pA++) = curGrid._e.imag();
				}
				// south
				if (nu[id - m + 1] > 0)
				{
					*(pA++) = cur_nu;
					*(pA++) = nu[id - m + 1];
					*(pA++) = curGrid._s.real();
					*(pA++) = curGrid._s.imag();
				}
				// south - east
				if (nu[id - m + 2] > 0)
				{
					*(pA++) = cur_nu;
					*(pA++) = nu[id - m + 2];
					*(pA++) = curGrid._se.real();
					*(pA++) = curGrid._se.imag();
				}
				// 将 l = n-1 时，在口径面剩余点处产生的积分值进行组装
				if (l == n - 1)
				{
					//sq = findApertureIndex(nbound, apertureNumber, cur_nu);
					curGrid._c += curGrid.tm[sq];
					for (int apeIndex = 0; apeIndex < apertureNumber; apeIndex++)
					{
						if (apeIndex<sq - 1 || apeIndex>sq + 1)//这个条件确保不会再次考虑已经组装的部分
						{
							*(pA++) = cur_nu;
							*(pA++) = nbound[0][apeIndex];
							*(pA++) = curGrid.tm[apeIndex].real();
							*(pA++) = curGrid.tm[apeIndex].imag();
						}
					}
				}
				*(pA++) = cur_nu;
				*(pA++) = cur_nu;
				*(pA++) = curGrid._c.real();
				*(pA++) = curGrid._c.imag();
			}

		}
	}

	//截去多余内存区域
	int ArowNum = (pA - A) / 4;
	double *A_result = (double*)malloc(ArowNum * 4 * sizeof(double));
	//memcpy(A_result, A, (pA - A)- sizeof(double));
	memcpy(A_result, A, ArowNum * 4 * sizeof(double));

	//构造matlab中对应的setA函数的返回值
	mxArray *mx_A;

	//将C中数组转化为matlab数组mxArray
	//注意matlab默认是列优先的
	//我们要生产ArowNum*3的矩阵，可以先构造3*ArowNum的矩阵再转置
	mx_A = mxCreateDoubleMatrix(3, ArowNum, mxCOMPLEX);
	double *AvaluePr = mxGetPr(mx_A);
	double *AvaluePi = mxGetPi(mx_A);
	for (int i = 0; i < ArowNum; i++)
	{
		//行索引
		*AvaluePr++ = *A_result++; ////这个索引是从1计数的(因此传给matlab不需要-1)
		*AvaluePi++ = 0;
		//列索引
		*AvaluePr++ = *A_result++;
		*AvaluePi++ = 0;
		//取值
		*AvaluePr++ = *A_result++;
		*AvaluePi++ = *A_result++;
	}
	return mx_A;
}


VectorXcd SingalCavity::setB(vector<vector<gridCell>> &grid, VectorXcd &rh, int nn, vector<int> &nu)
{
	int m = this->meshWidth;
	int n = this->meshHeight;
	double dx = this->stepX;
	double dy = this->stepY;
	double left = this->virtualBorderLeft;
	double bottom = this->virtualBorderBottom;

	VectorXcd B(nn);
	//B = zeros(nn, 1);
	
	gridCell curGrid;
	int id;
	for (int j = 0 ; j<m - 1;j++)
	{
		for (int l=0; l< n;l++)
		{
			curGrid = grid[j][l];
			id = j + l*(m - 1);
			if (nu[id] > 0)
			{
				B(nu[id]-1) = rh(id) - curGrid._const;
			}
		}
	}
	return B;
}


mxArray* SingalCavity::setB_mx(vector<vector<gridCell>> &grid, VectorXcd &rh, int nn, vector<int> &nu)
{
	int m = this->meshWidth;
	int n = this->meshHeight;
	double dx = this->stepX;
	double dy = this->stepY;
	double left = this->virtualBorderLeft;
	double bottom = this->virtualBorderBottom;

	VectorXcd B(nn);
	//B = zeros(nn, 1);

	gridCell curGrid;
	int id;
	for (int j = 0; j<m - 1; j++)
	{
		for (int l = 0; l< n; l++)
		{
			curGrid = grid[j][l];
			id = j + l*(m - 1);
			if (nu[id] > 0)
			{
				B(nu[id] - 1) = rh(id) - curGrid._const;
			}
		}
	}
	
	//构造matlab中对应的setB函数的返回值
	mxArray *mx_B;
	//将C中数组转化为matlab数组mxArray
	mx_B = mxCreateDoubleMatrix(B.size(), 1, mxCOMPLEX);
	double *BvaluePr = mxGetPr(mx_B);
	double *BvaluePi = mxGetPi(mx_B);
	for (int i = 0; i < B.size(); i++)
	{
		*BvaluePr++ = B(i).real();
		*BvaluePi++ = B(i).imag();
	}
	return mx_B;
}


void SingalCavity::assign(VectorXcd solution, TriangleMesh &U, TriangleMesh &L, vector<int> &nu)
{
	int m = this->meshWidth;
	int n = this->meshHeight;
	//g = para.g;
	//phi = para.phi;
	int id;
	for (int j = 0; j < m; j++)
	{
		for (int l = 0; l < n; l++)
		{
			int id = j + l*(m - 1);

			if (U.Get_sign(j, l))
			{
				Triangle_Normal *A = U.GetPt_normal(j, l);
				if (A->sgn == -10)
				{
					complex<double> a1, a2, a3;
					if (abs(phi(A->x(0), A->y(0))) < 1e-10)
						a1 = this->u(A->x(0), A->y(0));
					else
						a1 = solution(nu[id] - 1);

					if (abs(phi(A->x(1), A->y(1))) < 1e-10)
						a2 = this->u(A->x(1), A->y(1));
					else
						a2 = solution(nu[id - 1] - 1);

					if (abs(phi(A->x(2), A->y(2))) < 1e-10)
						a3 = this->u(A->x(2), A->y(2));
					else
						a3 = solution(nu[id - m + 1] - 1);

					A->z_int = Vector3cd(a1, a2, a3);
				}
			}
			else
			{
				Triangle_All *A = U.GetPt_all(j, l);
				if (abs(A->sgn) < 10)
				{
					int r = abs(A->sgn);
					int s = r % 3 + 1;
					int t = s % 3 + 1;
					if (A->sgn < 0)
					{
						if (r == 1)
							A->z1 = solution(nu[id] - 1);

						if (r == 2)
							A->z1 = solution(nu[id - 1] - 1);

						if (r == 3)
							A->z1 = solution(nu[id - m + 1] - 1);
					}
					else
					{
						if (A->sgn > 0)
						{
							if (s == 1)
								A->z2 = solution(nu[id] - 1);

							if (s == 2)
								A->z2 = solution(nu[id - 1] - 1);

							if (s == 3)
								A->z2 = solution(nu[id - m + 1] - 1);

							if (t == 1)
								A->z3 = solution(nu[id] - 1);

							if (t == 2)
								A->z3 = solution(nu[id - 1] - 1);

							if (t == 3)
								A->z3 = solution(nu[id - m + 1] - 1);
						}
					}
					A->z4 = this->u(A->x4, A->y4);
					A->z5 = this->u(A->x5, A->y5);
				}

				if (abs(A->sgn) > 10)
				{
					if (A->sgn < 0)
					{
						if (A->v(0) == 1)
							A->z1 = solution(nu[id] - 1);

						if (A->v(0) == 2)
							A->z1 = solution(nu[id - 1] - 1);

						if (A->v(0) == 3)
							A->z1 = solution(nu[id - m + 1] - 1);
					}
					if (A->sgn > 0)
					{
						if (A->v(2) == 1)
							A->z3 = solution(nu[id] - 1);

						if (A->v(2) == 2)
							A->z3 = solution(nu[id - 1] - 1);

						if (A->v(2) == 3)
							A->z3 = solution(nu[id - m + 1] - 1);
					}
					A->z2 = this->u(A->x2, A->y2);
					A->z5 = this->u(A->p5(0), A->p5(1));
				}
			}


			if (L.Get_sign(j, l))
			{
				Triangle_Normal *B = L.GetPt_normal(j, l);

				if (B->sgn == -10)
				{
					complex<double> a1, a2, a3;
					if (abs(phi(B->x(0), B->y(0))) < 1e-10)
						a1 = this->u(B->x(0), B->y(0));
					else
						a1 = solution(nu[id - m] - 1);

					if (abs(phi(B->x(1), B->y(1))) < 1e-10)
						a2 = this->u(B->x(1), B->y(1));
					else
						a2 = solution(nu[id - m + 1] - 1);

					if (abs(phi(B->x(2), B->y(2))) < 1e-10)
						a3 = this->u(B->x(2), B->y(2));
					else
						a3 = solution(nu[id - 1] - 1);

					B->z_int = Vector3cd(a1, a2, a3);
				}
			}
			else
			{
				Triangle_All *B = L.GetPt_all(j, l);

				if (abs(B->sgn) < 10)
				{
					int r = abs(B->sgn);
					int s = r % 3 + 1;
					int t = s % 3 + 1;
					if (B->sgn < 0)
					{
						if (r == 1)
							B->z1 = solution(nu[id - m] - 1);

						if (r == 2)
							B->z1 = solution(nu[id - m + 1] - 1);

						if (r == 3)
							B->z1 = solution(nu[id - 1] - 1);
					}
					else
					{
						if (B->sgn > 0)
						{
							if (s == 1)
								B->z2 = solution(nu[id - m] - 1);

							if (s == 2)
								B->z2 = solution(nu[id - m + 1] - 1);

							if (s == 3)
								B->z2 = solution(nu[id - 1] - 1);

							if (t == 1)
								B->z3 = solution(nu[id - m] - 1);

							if (t == 2)
								B->z3 = solution(nu[id - m + 1] - 1);

							if (t == 3)
								B->z3 = solution(nu[id - 1] - 1);
						}
					}
					B->z4 = this->u(B->x4, B->y4);
					B->z5 = this->u(B->x5, B->y5);
				}

				if (abs(B->sgn) > 10)
				{
					if (B->sgn < 0)
					{
						if (B->v(0) == 1)
							B->z1 = solution(nu[id - m] - 1);

						if (B->v(0) == 2)
							B->z1 = solution(nu[id - m + 1] - 1);

						if (B->v(0) == 3)
							B->z1 = solution(nu[id - 1] - 1);
					}
					if (B->sgn > 0)
					{
						if (B->v(2) == 1)
							B->z3 = solution(nu[id - m] - 1);

						if (B->v(2) == 2)
							B->z3 = solution(nu[id - m + 1] - 1);

						if (B->v(2) == 3)
							B->z3 = solution(nu[id - 1] - 1);
					}
					B->z2 = this->u(B->x2, B->y2);
					B->z5 = this->u(B->p5(0), B->p5(1));
				}
			}
			

		}// end for (int l = 0; l < n; l++)
	}
}


VectorXcd SingalCavity::getAperture(vector<vector<double>> &nbound, TriangleMesh &U, TriangleMesh &L)
{
	int m = this->meshWidth;
	int n = this->meshHeight;
	int aperture_Num = nbound[0].size();
	VectorXcd ApertureValue(aperture_Num);

	int index = 0;
	for (int j = 0; j < m; j++)
	{
		// l == n -1
		if (U.Get_sign(j, n - 1))
		{
			Triangle_Normal tri = U.Get_normal(j, n - 1);

			complex<double> z_int0;
			if (tri.x(0) == nbound[1][index])
			{
				z_int0 = tri.z_int(0);
				ApertureValue(index) = z_int0;
				index++;
				if (index >= aperture_Num)
					break;
			}
		}
		else
		{
			Triangle_All tri = U.Get_all(j, n - 1);

			complex<double> z_int0;
			if (tri.x(0) == nbound[1][index])
			{
				z_int0 = tri.z_int(0);
				ApertureValue(index) = z_int0;
				index++;
				if (index >= aperture_Num)
					break;
			}
		}
	}

	return ApertureValue;
}


void SingalCavity::drawAperture(VectorXd &plotX, VectorXd &plotY, vector<vector<double>> &nbound, int sign)
{
	int plotArraySize = this->solutionOfAperture.size() + 2;
	plotX = VectorXd(plotArraySize);
	plotY = VectorXd(plotArraySize);

	plotX(0) = this->apertureLeft;
	plotX(plotArraySize - 1) = this->apertureRight;
	double *nboundPt = &nbound[1][0];
	for (int i = 1; i < plotArraySize - 1; i++)
	{
		plotX(i) = *nboundPt++;
	}

	plotY(0) = 0;
	plotY(plotArraySize - 1) = 0;
	if (sign == 0)
	{
		complex<double> *aperPt = &solutionOfAperture(0);
		for (int i = 1; i < plotArraySize-1; i++)
		{
			plotY(i) = abs(*aperPt++);
		}
	}
	else if (sign == 1)
	{
		complex<double> *aperPt = &solutionOfAperture(0);
		for (int i = 1; i < plotArraySize - 1; i++)
		{
			plotY(i) = real(*aperPt++);
		}
	}
	else if (sign == -1)
	{
		complex<double> *aperPt = &solutionOfAperture(0);
		for (int i = 1; i < plotArraySize - 1; i++)
		{
			plotY(i) = imag(*aperPt++);
		}
	}
	else if (sign == 10)
	{
		complex<double> *aperPt = &solutionOfAperture(0);
		for (int i = 1; i < plotArraySize - 1; i++)
		{
			double re = real(*aperPt);
			double im = imag(*aperPt);
			plotY(i) = (-180 / M_PI)*atan2(im, re);
		}
	}
	else
	{
		;
	}
}


//int PDE::findApertureIndex(double *nbound, int nboundLen, double value)
//{
//int index;
//for (int i = 0; i < nboundLen; i++)
//{
//if (*(nbound++) == value)
//{
//index = i;
//break;
//}
//}
//return index;
//}


//判断并设置三角形类型
void SingalCavity::setTriangleType(double x1, double x2, double x3, double y1, double y2, double y3, int &sgn)
{
	double a1 = phi(x1, y1);
	double a2 = phi(x2, y2);
	double a3 = phi(x3, y3);

#pragma region sgn（标识三角形网格与腔体之间关系）

	if (abs(a1) < 1e-10 && abs(a2) < 1e-10 && abs(a3) < 1e-10)
	{
		//三角形三点都在界面上
		double ac = phi((x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3);
		if (ac > 1e-10)
			sgn = 10;
		if (ac < -1e-10)
			sgn = -10;
	}

	if ((abs(a1) < 1e-10 && a2 > 1e-10 && a3 > 1e-10) ||
		(a1 > 1e-10 && abs(a2) < 1e-10 && a3 > 1e-10) ||
		(a1 > 1e-10 && a2 > 1e-10 && abs(a3) < 1e-10) ||
		(abs(a1) < 1e-10 && abs(a2) < 1e-10 && a3 > 1e-10) ||
		(abs(a1) < 1e-10 && a2 > 1e-10 && abs(a3) < 1e-10) ||
		(a1 > 1e-10 && abs(a2) < 1e-10 && abs(a3) < 1e-10) ||
		(a1 > 1e-10 && a2 > 1e-10 && a3 > 1e-10))
	{
		// 对应三角形位于Ω + 中的情况：
		//一个点在界面上，两个点在Ω +
		//或两个点在界面上，一个点在Ω +
		//或三个点都在Ω +
		sgn = 10;
	}

	if ((abs(a1) < 1e-10 && a2 < -1e-10 && a3 < -1e-10) ||
		(a1 < -1e-10 && abs(a2) < 1e-10 && a3 < -1e-10) ||
		(a1 < -1e-10 && a2 < -1e-10 && abs(a3) < 1e-10) ||
		(abs(a1) < 1e-10 && abs(a2) < 1e-10 && a3 < -1e-10) ||
		(abs(a1) < 1e-10 && a2 < -1e-10 && abs(a3) < 1e-10) ||
		(a1 < -1e-10 && abs(a2) < 1e-10 && abs(a3) < 1e-10) ||
		(a1 < -1e-10 && a2 < -1e-10 && a3 < -1e-10))
	{
		// 对应三角形位于Ω - 中的情况：
		//一个点在界面上，两个点在Ω -
		//或两个点在界面上，一个点在Ω -
		//或三个点都在Ω -
		sgn = -10;
	}

	if (a1 > 1e-10 && a2 < -1e-10 && a3 < -1e-10)
	{
		//三角形越过界面
		// (x1, y1)位于Ω + (x2, y2)、(x3, y3)位于Ω -
		sgn = 1;
	}
	if (a1 < -1e-10 && a2>1e-10 && a3 > 1e-10)
	{
		// 三角形越过界面
		//a1位于Ω - a2、a3位于Ω +
		sgn = -1;
	}
	if (a2 > 1e-10 && a1 < -1e-10 && a3 < -1e-10)
	{
		//三角形越过界面
		//a2位于Ω + a1、a3位于Ω -
		sgn = 2;
	}
	if (a2<-1e-10 && a1>1e-10 && a3 > 1e-10)
	{
		sgn = -2;
	}
	if (a3 > 1e-10 && a2 < -1e-10 && a1 < -1e-10)
	{
		sgn = 3;
	}
	if (a3 < -1e-10 && a2>1e-10 && a1 > 1e-10)
	{
		sgn = -3;
	}

	//一个点在界面上，一个在Ω + ，一个在Ω -
	if (abs(a2) < 1e-10 && a1 > 1e-10 && a3 < -1e-10)
	{
		// a2在界面上，a1在Ω + ，a3在Ω -
		sgn = 11;
	}
	if (abs(a2) < 1e-10 && a1 < -1e-10 && a3>1e-10)
	{
		sgn = -11;
	}
	if (abs(a3) < 1e-10 && a2 > 1e-10 && a1 < -1e-10)
	{
		sgn = 12;
	}
	if (abs(a3) < 1e-10 && a2 < -1e-10 && a1>1e-10)
	{
		sgn = -12;
	}
	if (abs(a1) < 1e-10 && a2 < -1e-10 && a3>1e-10)
	{
		sgn = 13;
	}
	if (abs(a1) < 1e-10 && a2 > 1e-10 && a3 < -1e-10)
	{
		sgn = -13;
	}
#pragma endregion
}


void SingalCavity::tri(double x1, double x2, double x3, double y1, double y2, double y3, Triangle_Normal& L)
{
	L.x = Vector3d(x1, x2, x3);
	L.y = Vector3d(y1, y2, y3);

	//-------------------------------------
	// l=10 || l=-10 无需其他操作
	//-------------------------------------
}

void SingalCavity::tri(double x1, double x2, double x3, double y1, double y2, double y3, Triangle_All& L)
{
	L.x = Vector3d(x1, x2, x3);
	L.y = Vector3d(y1, y2, y3);

	//-------------------------------------
	// l=10 || l=-10 无需其他操作
	//-------------------------------------


	//-------------------------------------
	// abs(l)<10 
	//-------------------------------------
	if (abs(L.sgn) < 10)
	{
		if (L.sgn == 1)
		{
			//三角形越过界面
			// (x1, y1)位于Ω + (x2, y2)、(x3, y3)位于Ω -
			L.r1 = findzero(x1, y1, x2, y2);
			L.r2 = findzero(x1, y1, x3, y3);
		}
		if (L.sgn == -1)
		{
			// 三角形越过界面
			//a1位于Ω - a2、a3位于Ω +
			L.r1 = findzero(x1, y1, x2, y2);
			L.r2 = findzero(x1, y1, x3, y3);
		}
		if (L.sgn == 2)
		{
			//三角形越过界面
			//a2位于Ω + a1、a3位于Ω -
			L.r1 = findzero(x2, y2, x3, y3);
			L.r2 = findzero(x2, y2, x1, y1);
		}
		if (L.sgn == -2)
		{
			L.r1 = findzero(x2, y2, x3, y3);
			L.r2 = findzero(x2, y2, x1, y1);
		}
		if (L.sgn == 3)
		{
			L.r1 = findzero(x3, y3, x1, y1);
			L.r2 = findzero(x3, y3, x2, y2);
		}
		if (L.sgn == -3)
		{
			L.r1 = findzero(x3, y3, x1, y1);
			L.r2 = findzero(x3, y3, x2, y2);
		}
	}


	//-------------------------------------
	// abs(l)>10 
	//一个点在界面上，一个在Ω + ，一个在Ω -
	//-------------------------------------
	if (abs(L.sgn) > 10)
	{
		if (L.sgn == 11)
		{
			// a2在界面上，a1在Ω + ，a3在Ω -
			L.A = Vector2d(x1, y1);  //与B按照逆时针相对排列
			L.B = Vector2d(x2, y2);  //位于Γ上的点
			L.C = Vector2d(x3, y3);  //与B按照逆时针相对排列
			L.v = Vector3i(1, 2, 3); //标识ABC三点分别对应的点
			L.r1 = 1;
			L.r2 = findzero(x1, y1, x3, y3);
		}
		if (L.sgn == -11)
		{
			L.A = Vector2d(x1, y1);
			L.B = Vector2d(x2, y2);
			L.C = Vector2d(x3, y3);
			L.v = Vector3i(1, 2, 3);
			L.r1 = 1;
			L.r2 = findzero(x1, y1, x3, y3);
		}
		if (L.sgn == 12)
		{
			L.A = Vector2d(x2, y2);
			L.B = Vector2d(x3, y3);
			L.C = Vector2d(x1, y1);
			L.v = Vector3i(2, 3, 1);
			L.r1 = 1;
			L.r2 = findzero(x2, y2, x1, y1);
		}
		if (L.sgn == -12)
		{
			L.A = Vector2d(x2, y2);
			L.B = Vector2d(x3, y3);
			L.C = Vector2d(x1, y1);
			L.v = Vector3i(2, 3, 1);
			L.r1 = 1;
			L.r2 = findzero(x2, y2, x1, y1);
		}
		if (L.sgn == 13)
		{
			L.A = Vector2d(x3, y3);
			L.B = Vector2d(x1, y1);
			L.C = Vector2d(x2, y2);
			L.v = Vector3i(3, 1, 2);
			L.r1 = 1;
			L.r2 = findzero(x3, y3, x2, y2);
		}
		if (L.sgn == -13)
		{
			L.A = Vector2d(x3, y3);
			L.B = Vector2d(x1, y1);
			L.C = Vector2d(x2, y2);
			L.v = Vector3i(3, 1, 2);
			L.r1 = 1;
			L.r2 = findzero(x3, y3, x2, y2);
		}
	}
}

//按照统一规则标准化三角形坐标信息
void SingalCavity::solveTri(Triangle_All &L)
{
	// 按照统一规则标准化三角形坐标信息

	// sgn == 10时本函数不需单独处理
	// sgn == -10时的情况不需要考虑

	if (abs(L.sgn) < 10)
	{
		int p1 = abs(L.sgn);
		int p2 = p1 % 3 + 1;
		int p3 = p2 % 3 + 1;
		L.x1 = L.x(p1 - 1); L.y1 = L.y(p1 - 1);
		L.x2 = L.x(p2 - 1); L.y2 = L.y(p2 - 1);
		L.x3 = L.x(p3 - 1); L.y3 = L.y(p3 - 1);
		L.x4 = L.x1 + (L.x2 - L.x1)*L.r1; L.y4 = L.y1 + (L.y2 - L.y1)*L.r1;
		L.x5 = L.x1 + (L.x3 - L.x1)*L.r2; L.y5 = L.y1 + (L.y3 - L.y1)*L.r2;
		return;
	}
	if (abs(L.sgn) > 10)
	{
		L.p5 = L.A + (L.C - L.A)*L.r2;
		L.x5 = L.p5(0); L.y5 = L.p5(1);
		L.x1 = L.A(0); L.y1 = L.A(1);
		L.x2 = L.B(0); L.y2 = L.B(1);
		L.x3 = L.C(0); L.y3 = L.C(1);
		return;
	}
}