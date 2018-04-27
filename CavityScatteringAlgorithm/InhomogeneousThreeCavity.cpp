#include "InhomogeneousThreeCavity.h"


InhomogeneousThreeCavity::InhomogeneousThreeCavity(unsigned int cavityType) :Cavity(cavityType)
{
	int test = 0;
}

void InhomogeneousThreeCavity::InitElectromagneticParameter(double k0, complex<double>epr1, complex<double> epr1_int, complex<double>epr2, complex<double> epr2_int, complex<double>epr3, complex<double> epr3_int, double theta)
{
	this->k0 = k0;

	this->epr1 = epr1;
	this->epr1_int = epr1_int;
	this->epr2 = epr2;
	this->epr2_int = epr2_int;
	this->epr3 = epr3;
	this->epr3_int = epr3_int;
	this->k1_2 = k0*k0*epr1;
	this->k1_2int = k0*k0*epr1_int;
	this->k2_2 = k0*k0*epr2;
	this->k2_2int = k0*k0*epr2_int;
	this->k3_2 = k0*k0*epr3;
	this->k3_2int = k0*k0*epr3_int;
	this->theta = theta;

	//Mark function InitElectromagneticParameter has been executed
	this->initalCheckKey = this->initalCheckKey | 4; // 4 means 0000 0100
}


void InhomogeneousThreeCavity::InitAperture(double apertureY, double aperture1Left, double aperture1Right, double aperture2Left, double aperture2Right, double aperture3Left, double aperture3Right)
{
	this->cavity_number = 3;

	//set the parameter of aperture
	this->apertureY = apertureY;

	this->aperture1Left = aperture1Left;
	this->aperture1Right = aperture1Right;
	this->aperture2Left = aperture2Left;
	this->aperture2Right = aperture2Right;
	this->aperture3Left = aperture3Left;
	this->aperture3Right = aperture3Right;
	this->separator1_2 = (aperture1Right + aperture2Left) / 2;
	this->separator2_3 = (aperture2Right + aperture3Left) / 2;

	//Mark function InitCavityShapeParameter has been executed
	this->initalCheckKey = this->initalCheckKey | 8; // 8 means 0000 1000
}


bool InhomogeneousThreeCavity::Solve()
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
	//int nn;
	//vector<int> nu;
	//vector<vector<double>> nbound;
	setTri(U, L, nn, nu, nbound);
	myTimer.EndAndPrint();


	this->G1_aperture = ApertureIntegral::computeG(nbound[0].size() / 3, this->k0, this->aperture1Left, this->aperture1Right);
	this->G2_aperture = ApertureIntegral::computeG(nbound[0].size() / 3, this->k0, this->aperture2Left, this->aperture2Right);
	this->G3_aperture = ApertureIntegral::computeG(nbound[0].size() / 3, this->k0, this->aperture3Left, this->aperture3Right);

	this->g_aperture = compute_g(G1_aperture, nbound);   //注意这里会变

	myTimer.Start("setGrid");
	vector<vector<gridCell>> gridCell = setGrid(U, L, nbound, nu);
	myTimer.EndAndPrint();

	myTimer.Start("setRightHand");
	VectorXcd rh = setRightHand(U, L, nu);
	//cout << rh;
	myTimer.EndAndPrint();

	myTimer.Start("setA");
	SparseMatrix<complex<double>> A = setA(nn, nu, nbound, gridCell);
	myTimer.EndAndPrint();

	myTimer.Start("setB");
	VectorXcd B = setB(gridCell, rh, nn, nu);
	myTimer.EndAndPrint();

	myTimer.Start("solveX");
	VectorXcd x = solveX(A, B);
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






/*  --------------------------------------
--------------------------------------
protect
--------------------------------------
--------------------------------------
*/


Matrix2d InhomogeneousThreeCavity::beta(double x, double y)
{
	Matrix2d beta;
	beta(0, 0) = 1;
	beta(1, 0) = 0;
	beta(0, 1) = 0;
	beta(1, 1) = 1;
	return beta;
}

Matrix2d InhomogeneousThreeCavity::beta_(double x, double y)
{
	Matrix2d beta;
	beta(0, 0) = 1;
	beta(1, 0) = 0;
	beta(0, 1) = 0;
	beta(1, 1) = 1;
	return beta;
}

complex<double> InhomogeneousThreeCavity::q(double x, double y)
{
	complex<double> q;

	if( x<separator1_2)
		q = -k1_2;
	else if( x >= separator1_2 && x <= separator2_3)
		q = -k2_2;
	else
		q = -k3_2;

	return q;
}

complex<double> InhomogeneousThreeCavity::q_(double x, double y)
{
	complex<double> q_;

	if (x<separator1_2)
		q_ = -k1_2int;
	else if (x >= separator1_2 && x <= separator2_3)
		q_ = -k2_2int;
	else
		q_ = -k3_2int;

	return q_;
}

complex<double> InhomogeneousThreeCavity::f(double x, double y)
{
	complex<double> f(0, 0);
	return f;
}

complex<double> InhomogeneousThreeCavity::f_(double x, double y)
{
	complex<double> f_(0, 0);
	return f_;
}

complex<double> InhomogeneousThreeCavity::u(double x, double y)
{
	//腔体内精确解的值，用于获取边界条件或考察算法精确度
	complex<double> u(0, 0);
	return u;
}

complex<double> InhomogeneousThreeCavity::u_(double x, double y)
{
	//腔体内非均匀介质精确解的值，用于获取边界条件或考察算法精确度
	complex<double> u_(0, 0);
	return u_;
}

double InhomogeneousThreeCavity::a(double x, double y)
{
	double out = 0;
	return out;
}

Vector2d InhomogeneousThreeCavity::b(double x, double y)
{
	Vector2d out(0, 0);
	return out;
}

void InhomogeneousThreeCavity::setTri(TriangleMesh &U, TriangleMesh &L, int &nn, vector<int> &nu, vector<vector<double>> &nbound)
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
			if (a4 < 0 && j < m)
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

	nbound.resize(2);
	nbound[0].resize(ApertureNum);
	nbound[1].resize(ApertureNum);

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
			// |  \   |
			// | L \ |
			// |     \|
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
				if (a4 < 0)
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

			//_________________________________________
			//利用腔体边界的水平集函数，设置当前三角形与腔体边界关系的相关信息（腔体内？与腔体边界相交？相交位置？）
			//_________________________________________
			//判断并设置三角形类型，与界面交点位置等信息
			Triangle l = tri(a1, a2, a3, x1, x2, x3, y1, y2, y3);
			Triangle u = tri(a4, a3, a2, x4, x3, x2, y4, y3, y2);
			//按照统一规则标准化三角形坐标信息
			solveTri(l);
			solveTri(u);

			//_________________________________________
			//利用非均匀介质的水平集函数，设置当前三角形与非均匀介质边界关系的相关信息
			//_________________________________________
			if (l.sgn != 10)  // l.sgn!= 10代表三角形在关心的腔体区域内或边界处，l.sgn == 10则完全不需要考虑
			{
				double aa1, aa2, aa3, aa4;
				if (l.sgn == -10)// 三角形完全在腔体内部
				{
					aa1 = phi_int(x1, y1);
					aa2 = phi_int(x2, y2);
					aa3 = phi_int(x3, y3);
					tri_int(aa1, aa2, aa3, x1, x2, x3, y1, y2, y3, l);
					l.sgn_intex = 0;
					//l.sgn_intex对应三角形被界面切割后，切下的四边形位于腔体内部时，形成的两个三角形中，标识第二个三角形的类型信息
					//因此除了，abs(l.sgn) < 10 && l.sgn > 0 的情况，
					//l.sgn_intex均为0，即不存在第二个需要标识的三角形
				}
				if (abs(l.sgn) < 10) // 界面割三角形两条边
				{
					//三角形被界面切割后会形成一个三角形和一个四边形
					if (l.sgn < 0) // 三角形被界面切割后，三角形位于腔体内部
					{
						aa1 = phi_int(l.x1, l.y1);
						aa2 = phi_int(l.x4, l.y4);
						aa3 = phi_int(l.x5, l.y5);
						tri_int(aa1, aa2, aa3, l.x1, l.x4, l.x5, l.y1, l.y4, l.y5, l);
						l.sgn_intex = 0;
					}
					if (l.sgn > 0) // 三角形被界面切割后，四边形位于腔体内部
					{
						//四边形又分为两个三角形处理
						aa1 = phi_int(l.x2, l.y2);
						aa2 = phi_int(l.x3, l.y3);
						aa3 = phi_int(l.x4, l.y4);
						aa4 = phi_int(l.x5, l.y5);
						tri_int(aa1, aa2, aa3, l.x2, l.x3, l.x4, l.y2, l.y3, l.y4, l);
						tri_int1(aa2, aa4, aa3, l.x3, l.x5, l.x4, l.y3, l.y5, l.y4, l);
					}
				}
				if (abs(l.sgn) > 10) // 界面割三角形一条边，一个顶点
				{
					if (l.sgn < 0)
					{
						aa1 = phi_int(l.x1, l.y1);
						aa2 = phi_int(l.x2, l.y2);
						aa3 = phi_int(l.x5, l.y5);
						tri_int(aa1, aa2, aa3, l.x1, l.x2, l.x5, l.y1, l.y2, l.y5, l);
						l.sgn_intex = 0;
					}
					if (l.sgn > 0)
					{
						aa1 = phi_int(l.x3, l.y3);
						aa2 = phi_int(l.x5, l.y5);
						aa3 = phi_int(l.x2, l.y2);
						tri_int(aa1, aa2, aa3, l.x3, l.x5, l.x2, l.y3, l.y5, l.y2, l);
						l.sgn_intex = 0;
					}
				}
				solveTri_int(l);
				if (l.sgn_intex != 0)
				{
					solveTri_int1(l);
				}
			} //end of:  l.sgn!= 10

			if (u.sgn != 10)
			{
				double aa1, aa2, aa3, aa4;
				if (u.sgn == -10)
				{
					aa1 = phi_int(x4, y4);
					aa2 = phi_int(x3, y3);
					aa3 = phi_int(x2, y2);
					tri_int(aa1, aa2, aa3, x4, x3, x2, y4, y3, y2, u);
					u.sgn_intex = 0;
				}
				if (abs(u.sgn) < 10)
				{
					if (u.sgn < 0)
					{
						aa1 = phi_int(u.x1, u.y1);
						aa2 = phi_int(u.x4, u.y4);
						aa3 = phi_int(u.x5, u.y5);
						tri_int(aa1, aa2, aa3, u.x1, u.x4, u.x5, u.y1, u.y4, u.y5, u);
						u.sgn_intex = 0;
					}
					if (u.sgn > 0)
					{
						aa1 = phi_int(u.x2, u.y2);
						aa2 = phi_int(u.x3, u.y3);
						aa3 = phi_int(u.x4, u.y4);
						aa4 = phi_int(u.x5, u.y5);
						tri_int(aa1, aa2, aa3, u.x2, u.x3, u.x4, u.y2, u.y3, u.y4, u);
						tri_int1(aa2, aa4, aa3, u.x3, u.x5, u.x4, u.y3, u.y5, u.y4, u);
					}
				}
				if (abs(u.sgn) > 10)
				{
					if (u.sgn < 0)
					{
						aa1 = phi_int(u.x1, u.y1);
						aa2 = phi_int(u.x2, u.y2);
						aa3 = phi_int(u.x5, u.y5);
						tri_int(aa1, aa2, aa3, u.x1, u.x2, u.x5, u.y1, u.y2, u.y5, u);
						u.sgn_intex = 0;
					}
					if (u.sgn > 0)
					{
						aa1 = phi_int(u.x3, u.y3);
						aa2 = phi_int(u.x5, u.y5);
						aa3 = phi_int(u.x2, u.y2);
						tri_int(aa1, aa2, aa3, u.x3, u.x5, u.x2, u.y3, u.y5, u.y2, u);
						u.sgn_intex = 0;
					}
				}
				solveTri_int(u);
				if (u.sgn_intex != 0)
				{
					solveTri_int1(u);
				}
			}//end of:  u.sgn!= 10


			 // -------------------------------------
			U.Set(j, k, u);
			L.Set(j, k, l);

		}//end of: for (int k = 0; k < n; k++)
	}//end of: for (int j = 0; j < m; j++)
}


VectorXcd InhomogeneousThreeCavity::compute_g(MatrixXcd &G, vector<vector<double>> &nbound)
{
	//InhomogeneousSingalCavity中实现的是无精确解时的求解方法
	//当子类对应于有精确解的方法时，需覆写本方法

	double alpha = this->k0 * sin(this->theta);
	double beta = this->k0 * cos(this->theta);

	int size = nbound[1].size();
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


vector<vector<gridCell>> InhomogeneousThreeCavity::setGrid(TriangleMesh &U, TriangleMesh &L, vector<vector<double>> &nbound, vector<int> &nu)
{
	int m = this->meshWidth;
	int n = this->meshHeight;
	/*double dx = this->stepX;
	double dy = this->stepY;

	double left = this->virtualBorderLeft;
	double bottom = this->virtualBorderBottom;*/

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
					weak(U.Get(j, l), Vector3d(1, 0, 0), 1, nbound, t4, t4m);
					weak(L.Get(j + 1, l), Vector3d(0, 0, 1), 0, nbound, t5, t5m);
					weak(U.Get(j + 1, l), Vector3d(0, 1, 0), 1, nbound, t6, t6m);
				}
				else if (l == n - 2) // 口径面下面一排
				{
					//[1, 0, 0][0, 1, 0][0, 0, 1]用于标识当前点是三角形的哪个点
					weak(L.Get(j + 1, l + 1), Vector3d(1, 0, 0), 0, nbound, t1, t1m);
					weak(U.Get(j, l + 1), Vector3d(0, 0, 1), 1, nbound, t2, t2m);
					weak(L.Get(j, l + 1), Vector3d(0, 1, 0), 0, nbound, t3, t3m);
					weak(U.Get(j, l), Vector3d(1, 0, 0), 0, nbound, t4, t4m);
					weak(L.Get(j + 1, l), Vector3d(0, 0, 1), 0, nbound, t5, t5m);
					weak(U.Get(j + 1, l), Vector3d(0, 1, 0), 0, nbound, t6, t6m);
				}
				else
				{
					weak(L.Get(j + 1, l + 1), Vector3d(1, 0, 0), 0, nbound, t1, t1m);
					weak(U.Get(j, l + 1), Vector3d(0, 0, 1), 0, nbound, t2, t2m);
					weak(L.Get(j, l + 1), Vector3d(0, 1, 0), 0, nbound, t3, t3m);
					weak(U.Get(j, l), Vector3d(1, 0, 0), 0, nbound, t4, t4m);
					weak(L.Get(j + 1, l), Vector3d(0, 0, 1), 0, nbound, t5, t5m);
					weak(U.Get(j + 1, l), Vector3d(0, 1, 0), 0, nbound, t6, t6m);
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

void InhomogeneousThreeCavity::weak(Triangle &T, Vector3d v, int topsign, vector<vector<double>> &nbound, vector <complex<double>> &out, vector <complex<double>> &out5m)
{
	int nboundSize = nbound[0].size();
	complex<double> zero(0, 0);

	out.resize(4);
	initializeComplexVector(out, zero, 4);
	out5m.resize(nboundSize);
	initializeComplexVector(out5m, zero, nboundSize);

	double top = this->apertureY;

	if (T.sgn == -10)
	{
		Vector4cd w1 = weak_int(T, v);

		complex<double> out4(0, 0);
		//若三角形有顶点位于不规则边界上，计算不规则边界的积分
		if (abs(phi(T.x(0), T.y(0))) < 1e-10)
		{
			out4 = out4 + w1(0)*value(T.x(0), T.y(0));
			w1(0) = 0;
		}
		if (abs(phi(T.x(1), T.y(1))) < 1e-10)
		{
			out4 = out4 + w1(1)*value(T.x(1), T.y(1));
			w1(1) = 0;
		}
		if (abs(phi(T.x(2), T.y(2))) < 1e-10)
		{
			out4 = out4 + w1(2)*value(T.x(2), T.y(2));
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
		out[3] = w1(3) + out4 - out5;
		return;
	}//end if :  (T.sgn == -10)

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

		int vp1 = v(j - 1);
		int vp2 = v(k - 1);
		int vp3 = v(l - 1);
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

		Triangle T1, T2, T3;
		T1.x = Vector3d(x1, x4, x5);
		T1.y = Vector3d(y1, y4, y5);
		T2.x = Vector3d(x2, x3, x4);
		T2.y = Vector3d(y2, y3, y4);
		T3.x = Vector3d(x3, x5, x4);
		T3.y = Vector3d(y3, y5, y4);

		complex<double> u4 = value(x4, y4);
		complex<double> u5 = value(x5, y5);

		complex<double> c1, c2, c3, c4;
		if (T.sgn > 0)
		{
			reloadtri(T2, T);
			reloadtriex(T3, T);
			Vector4cd out2 = weak_int(T2, v2);
			Vector4cd out3 = weak_int(T3, v3) + solveweak4(T3, vp5, vp4);

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
			c2 = out2(0);
			c3 = out2(1) + out3(0);
			c4 = out2(2)*u4 + out3(1)*u5 + out3(2)*u4 + out2(3) + out3(3) - out5;
		}
		else
		{
			reloadtri(T1, T);
			Vector4cd out1 = weak_int(T1, v1) + solveweak4(T1, vp4, vp5);
			complex<double> out5(0, 0);
			if (topsign == 1)
			{
				if (abs(y1 - top) < 1e-10&&abs(y2 - top) < 1e-10)
					weak5(x1, y1, x4, y4, vp1, vp4, nbound, out5m, out5);
				else if (abs(y3 - top) < 1e-10&&abs(y1 - top) < 1e-10)
					weak5(x5, y5, x1, y1, vp5, vp1, nbound, out5m, out5);
			}
			c1 = out1(0);
			c2 = 0;
			c3 = 0;
			c4 = out1(1)*u4 + out1(2)*u5 + +out1(3) - out5;
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
	}//end of : (abs(T.sgn) < 10)

	if (abs(T.sgn) > 10)
	{
		double vp1 = v(T.v(0) - 1);
		double vp2 = v(T.v(1) - 1);
		double vp3 = v(T.v(2) - 1);
		double vp5 = vp1 + (vp3 - vp1)*T.r2;

		Vector3d v1(vp1, vp2, vp5);
		Vector3d v2(vp3, vp5, vp2);

		double x1 = T.x1;
		double y1 = T.y1;
		double x2 = T.x2;
		double y2 = T.y2;
		double x3 = T.x3;
		double y3 = T.y3;
		double x5 = T.x5;
		double y5 = T.y5;

		Triangle T1, T2;
		T1.x = Vector3d(x1, x2, x5);
		T1.y = Vector3d(y1, y2, y5);
		T2.x = Vector3d(x3, x5, x2);
		T2.y = Vector3d(y3, y5, y2);

		complex<double> u2 = value(x2, y2);
		complex<double> u5 = value(x5, y5);

		complex<double> c1, c2, c3, c4;
		if (T.sgn > 0)
		{
			reloadtri(T2, T);
			Vector4cd out2 = weak_int(T2, v2) + solveweak4(T2, vp5, vp2);
			complex<double> out5(0, 0);
			if (topsign == 1)
			{
				if (abs(y2 - top) < 1e-10&&abs(y3 - top) < 1e-10)
					weak5(x2, y2, x3, y3, vp2, vp3, nbound, out5m, out5);
				else if (abs(y3 - top) < 1e-10&&abs(y1 - top) < 1e-10)
					weak5(x3, y3, x5, y5, vp3, vp5, nbound, out5m, out5);
			}
			c1 = 0;
			c2 = 0;
			c3 = out2(0);
			c4 = out2(1)*u5 + out2(2)*u2 + out2(3) - out5;
		}
		if (T.sgn < 0)
		{
			reloadtri(T1, T);
			Vector4cd out1 = weak_int(T1, v1) + solveweak4(T1, vp2, vp5);
			complex<double> out5(0, 0);
			if (topsign == 1)
			{
				if (abs(y1 - top) < 1e-10&&abs(y2 - top) < 1e-10)
					weak5(x1, y1, x2, y2, vp1, vp2, nbound, out5m, out5);
				else if (abs(y3 - top) < 1e-10&&abs(y1 - top) < 1e-10)
					weak5(x5, y5, x1, y1, vp5, vp1, nbound, out5m, out5);
			}
			c1 = out1(0);
			c2 = 0;
			c3 = 0;
			c4 = out1(1)*u2 + out1(2)*u5 + out1(3) - out5;
		}

		out[T.v(0) - 1] = c1;
		out[T.v(1) - 1] = c2;
		out[T.v(2) - 1] = c3;
		out[3] = c4;

		return;
	}//end of :  (abs(T.sgn) > 10)
}

Vector4cd InhomogeneousThreeCavity::weak_int(Triangle &T, Vector3d v)
{
	if (T.sgn_int == 10)
	{
		Vector4cd out = weak1_int(T, v, true) + weak3_int(T, v, true);
		return out;
	}

	if (T.sgn_int == -10)
	{
		Vector4cd w1 = weak1_int(T, v, false) + weak3_int(T, v, false);
		complex<double> out4 = 0;
		if (abs(phi_int(T.x_int(0), T.y_int(0))) < 1e-10 && abs(phi_int(T.x_int(1), T.y_int(1))) < 1e-10)
		{
			out4 = weak4_int(Vector2d(T.x_int(0), T.y_int(0)), Vector2d(T.x_int(1), T.y_int(1)), b((T.x_int(1) + T.x_int(0)) / 2, (T.y_int(1) + T.y_int(0)) / 2), v(0), v(1));
		}
		if (abs(phi_int(T.x_int(1), T.y_int(1))) < 1e-10 && abs(phi_int(T.x_int(2), T.y_int(2))) < 1e-10)
		{
			out4 = weak4_int(Vector2d(T.x_int(1), T.y_int(1)), Vector2d(T.x_int(2), T.y_int(2)), b((T.x_int(1) + T.x_int(2)) / 2, (T.y_int(1) + T.y_int(2)) / 2), v(1), v(2));
		}
		if (abs(phi_int(T.x_int(2), T.y_int(2))) < 1e-10 && abs(phi_int(T.x_int(0), T.y_int(0))) < 1e-10)
		{
			out4 = weak4_int(Vector2d(T.x_int(2), T.y_int(2)), Vector2d(T.x_int(0), T.y_int(0)), b((T.x_int(2) + T.x_int(0)) / 2, (T.y_int(2) + T.y_int(0)) / 2), v(2), v(0));
		}
		Vector4cd out(w1(0), w1(1), w1(2), w1(3) + out4);
		return out;
	}//end if :  (T.sgn_int == -10)

	if (abs(T.sgn_int) < 10)
	{
		//Matrix2d (*beta1)(double x, double y) ;            //定义与函数beta参数类型匹配的函数指针
		//Matrix2d (*beta2)(double x, double y) ;            //定义与函数beta_参数类型匹配的函数指针
		//complex<double> (*q1)(double x, double y);    //定义与函数q参数类型匹配的函数指针
		//complex<double> (*q2)(double x, double y);    //定义与函数q_参数类型匹配的函数指针

		bool beta1, beta2;
		bool q1, q2;
		int s;
		if (T.sgn_int > 0)
		{
			beta1 = true;
			beta2 = false;
			q1 = true;
			q2 = false;
			s = T.sgn_int;
		}
		else
		{
			beta1 = false;
			beta2 = true;
			q1 = false;
			q2 = true;
			s = -T.sgn_int;
		}

		double r1 = T.r1_int;
		double r2 = T.r2_int;

		int j = s;
		int k = j % 3 + 1;
		int l = k % 3 + 1;

		int vp1 = v(j - 1);
		int vp2 = v(k - 1);
		int vp3 = v(l - 1);
		double vp4 = vp1 + (vp2 - vp1)*r1;
		double vp5 = vp1 + (vp3 - vp1)*r2;

		Vector3d v1(vp1, vp4, vp5);
		Vector3d v2(vp2, vp3, vp4);
		Vector3d v3(vp3, vp5, vp4);

		double x1 = T.x1_int;
		double y1 = T.y1_int;
		double x2 = T.x2_int;
		double y2 = T.y2_int;
		double x3 = T.x3_int;
		double y3 = T.y3_int;
		double x4 = T.x4_int;
		double y4 = T.y4_int;
		double x5 = T.x5_int;
		double y5 = T.y5_int;
		double x6 = (x4 + x5) / 2;
		double y6 = (y4 + y5) / 2;


		Triangle T1, T2, T3;
		T1.x = Vector3d(x1, x4, x5);
		T1.y = Vector3d(y1, y4, y5);
		T2.x = Vector3d(x2, x3, x4);
		T2.y = Vector3d(y2, y3, y4);
		T3.x = Vector3d(x3, x5, x4);
		T3.y = Vector3d(y3, y5, y4);
		T1.sgn_int = T.sgn_int;
		T2.sgn_int = T.sgn_int;
		T3.sgn_int = T.sgn_int;

		Vector4d u4 = T.u4_int;
		Vector4d u5 = T.u5_int;
		Vector4d u4_ = T.u4_int_;
		Vector4d u5_ = T.u5_int_;

		Vector4cd out1 = weak1_int(T1, v1, beta1) + weak3_int(T1, v1, q1);
		Vector4cd out2 = weak1_int(T2, v2, beta2) + weak3_int(T2, v2, q2);
		Vector4cd out3 = weak1_int(T3, v3, beta2) + weak3_int(T3, v3, q2);
		double out4;
		if (T.sgn_int > 0)
		{
			out4 = weak4_int(Vector2d(x5, y5), Vector2d(x4, y4), b(x6, y6), vp5, vp4);
		}
		else
		{
			out4 = weak4_int(Vector2d(x4, y4), Vector2d(x5, y5), b(x6, y6), vp4, vp5);
		}
		complex<double> c1, c2, c3, c4;
		c1 = out1(0) + out1(1)*u4(0) + out1(2)*u5(0) + out2(2)*u4_(0) + out3(1)*u5_(0) + out3(2)*u4_(0);
		c2 = out1(1)*u4(1) + out1(2)*u5(1) + out2(0) + out2(2)*u4_(1) + out3(1)*u5_(1) + out3(2)*u4_(1);
		c3 = out1(1)*u4(2) + out1(2)*u5(2) + out2(1) + out2(2)*u4_(2) + out3(0) + out3(1)*u5_(2) + out3(2)*u4_(2);
		c4 = out1(1)*u4(3) + out1(2)*u5(3) + out2(2)*u4_(3) + out3(1)*u5_(3) + out3(2)*u4_(3) + out4;

		if (s == 1)
		{
			Vector4cd out(c1, c2, c3, c4);
			return out;
		}
		if (s == 2)
		{
			Vector4cd out(c3, c1, c2, c4);
			return out;
		}
		if (s == 3)
		{
			Vector4cd out(c2, c3, c1, c4);
			return out;
		}
	}//end if :  (abs(T.sgn_int) < 10)

	if (abs(T.sgn_int) > 10)
	{

		int vp1 = v(T.v_int(0) - 1);
		int vp2 = v(T.v_int(1) - 1);
		int vp3 = v(T.v_int(2) - 1);
		double vp5 = vp1 + (vp3 - vp1)*T.r2_int;

		Vector3d v1(vp1, vp2, vp5);
		Vector3d v2(vp3, vp5, vp2);

		double x1 = T.x1_int;
		double y1 = T.y1_int;
		double x2 = T.x2_int;
		double y2 = T.y2_int;
		double x3 = T.x3_int;
		double y3 = T.y3_int;
		double x5 = T.x5_int;
		double y5 = T.y5_int;
		double x6 = (x2 + x5) / 2;
		double y6 = (y2 + y5) / 2;

		Triangle T1, T2;
		T1.x = Vector3d(x1, x2, x5);
		T1.y = Vector3d(y1, y2, y5);
		T2.x = Vector3d(x3, x5, x2);
		T2.y = Vector3d(y3, y5, y2);

		T1.sgn_int = T.sgn_int;
		T2.sgn_int = -T.sgn_int;
		T1.n = 1;
		T2.n = 2;

		Vector4d u5 = T.u5_int;
		Vector4d u5_ = T.u5_int_;

		complex<double> c1, c2, c3, c4;
		if (T.sgn_int > 0)
		{
			Vector4cd out1 = weak1_int(T1, v1, true) + weak3_int(T1, v1, true);
			Vector4cd out2 = weak1_int(T2, v2, false) + weak3_int(T2, v2, false);
			double out4 = weak4_int(Vector2d(x5, y5), Vector2d(x2, y2), b(x6, y6), vp5, vp2);

			c1 = out1(0) + out1(2)*u5(0) + out2(1)*u5_(0);
			c2 = out1(1) + out1(2)*u5(1) + out2(2) + out2(1)*u5_(1);
			c3 = out1(2)*u5(2) + out2(0) + out2(1)*u5_(2);
			c4 = out1(2)*u5(3) + out2(1)*u5_(3) + out1(3) + out2(3) + out4;
		}
		if (T.sgn_int < 0)
		{
			Vector4cd out1 = weak1_int(T1, v1, false) + weak3_int(T1, v1, false);
			Vector4cd out2 = weak1_int(T2, v2, true) + weak3_int(T2, v2, true);
			double out4 = weak4_int(Vector2d(x2, y2), Vector2d(x5, y5), b(x6, y6), vp2, vp5);

			c1 = out1(0) + out1(2)*u5_(0) + out2(1)*u5(0);
			c2 = out1(1) + out1(2)*u5_(1) + out2(2) + out2(1)*u5(1);
			c3 = out1(2)*u5_(2) + out2(0) + out2(1)*u5(2);
			c4 = out1(2)*u5_(3) + out2(1)*u5(3) + out1(3) + out2(3) + out4;
		}
		Vector4cd out;
		out(T.v_int(0) - 1) = c1;
		out(T.v_int(1) - 1) = c2;
		out(T.v_int(2) - 1) = c3;
		out(3) = c4;
		return out;
	}//end if :  (abs(T.sgn_int) > 10)
}

Vector3cd InhomogeneousThreeCavity::weak1(Triangle &T, Vector3d v)
{
	throw("wait to finish");
}

Vector4cd InhomogeneousThreeCavity::weak1_int(Triangle &T, Vector3d v, bool beta)
{
	// beta:    true->beta     false->beta_

	//phi = para.phi_int;
	double s = calculateTriangleArea(T.x, T.y);//求三角形面积
	Matrix3d M;
	M << T.x(0), T.y(0), 1, T.x(1), T.y(1), 1, T.x(2), T.y(2), 1;

	Matrix3d invM = M.inverse(); //求M的逆

	MatrixXd MM(2, 3);
	MM = invM.block(0, 0, 2, 3);

	Vector2d vdiff = MM*v;

	double x4 = (T.x(1) + T.x(2)) / 2;
	double x5 = (T.x(2) + T.x(0)) / 2;
	double x6 = (T.x(0) + T.x(1)) / 2;
	double y4 = (T.y(1) + T.y(2)) / 2;
	double y5 = (T.y(2) + T.y(0)) / 2;
	double y6 = (T.y(0) + T.y(1)) / 2;

	Matrix2d b4, b5, b6;
	if (beta)
	{
		b4 = this->beta(x4, y4);
		b5 = this->beta(x5, y5);
		b6 = this->beta(x6, y6);
	}
	else
	{
		b4 = this->beta_(x4, y4);
		b5 = this->beta_(x5, y5);
		b6 = this->beta_(x6, y6);
	}

	Vector3cd inte = vdiff.transpose() *(b4 + b5 + b6)*MM*s / 3;

	complex<double> con = 0;
	if (T.sgn_int == 10)
	{
		Vector4cd out(inte(0), inte(1), inte(2), 0);
		return out;
	}

	if (T.sgn_int == -10)
	{
		if (abs(phi_int(T.x(0), T.y(0))) < 1e-10)
		{
			con = con - inte(0)*a(T.x(0), T.y(0));
		}
		if (abs(phi_int(T.x(1), T.y(1))) < 1e-10)
		{
			con = con - inte(1)*a(T.x(1), T.y(1));
		}
		if (abs(phi_int(T.x(2), T.y(2))) < 1e-10)
		{
			con = con - inte(2)*a(T.x(2), T.y(2));
		}
		Vector4cd out(inte(0), inte(1), inte(2), con);
		return out;
	}

	if (abs(T.sgn_int) < 10)
	{
		Vector4cd out(inte(0), inte(1), inte(2), 0);
		return out;
	}

	if (abs(T.sgn_int) > 10)
	{
		if (T.sgn_int > 0)
		{
			Vector4cd out(inte(0), inte(1), inte(2), 0);
			return out;
		}
		if (T.sgn_int < 0)
		{
			if (T.n == 1) {
				con = con - inte(1)*a(T.x(1), T.y(1));
			}
			if (T.n == 2) {
				con = con - inte(2)*a(T.x(2), T.y(2));
			}
			Vector4cd out(inte(0), inte(1), inte(2), con);
			return out;
		}
	}
}

Vector3cd InhomogeneousThreeCavity::weak3(Triangle &T, Vector3d v)
{
	throw("wait to finish");
}

Vector4cd InhomogeneousThreeCavity::weak3_int(Triangle &T, Vector3d v, bool q)
{
	// q:    true->q     false->q_

	//phi = para.phi_int;

	double s = calculateTriangleArea(T.x, T.y);//求三角形面积

	double x4 = (T.x(1) + T.x(2)) / 2;
	double x5 = (T.x(2) + T.x(0)) / 2;
	double x6 = (T.x(0) + T.x(1)) / 2;
	double y4 = (T.y(1) + T.y(2)) / 2;
	double y5 = (T.y(2) + T.y(0)) / 2;
	double y6 = (T.y(0) + T.y(1)) / 2;

	complex<double> q4, q5, q6;
	if (q)
	{
		q4 = this->q(x4, y4);
		q5 = this->q(x5, y5);
		q6 = this->q(x6, y6);
	}
	else
	{
		q4 = this->q_(x4, y4);
		q5 = this->q_(x5, y5);
		q6 = this->q_(x6, y6);
	}

	double v4 = (v(1) + v(2)) / 2;
	double v5 = (v(2) + v(0)) / 2;
	double v6 = (v(0) + v(1)) / 2;

	complex<double> u1 = (q5*v5 + q6*v6)*s *(1.0 / 6);
	complex<double> u2 = (q6*v6 + q4*v4)*s *(1.0 / 6);
	complex<double> u3 = (q4*v4 + q5*v5)*s *(1.0 / 6);

	Vector3cd inte(u1, u2, u3);

	complex<double> con = 0;
	if (T.sgn_int == 10)
	{
		Vector4cd out(inte(0), inte(1), inte(2), 0);
		return out;
	}

	if (T.sgn_int == -10)
	{
		if (abs(phi_int(T.x(0), T.y(0))) < 1e-10)
		{
			con = con - inte(0)*a(T.x(0), T.y(0));
		}
		if (abs(phi_int(T.x(1), T.y(1))) < 1e-10)
		{
			con = con - inte(1)*a(T.x(1), T.y(1));
		}
		if (abs(phi_int(T.x(2), T.y(2))) < 1e-10)
		{
			con = con - inte(2)*a(T.x(2), T.y(2));
		}
		Vector4cd out(inte(0), inte(1), inte(2), con);
		return out;
	}

	if (abs(T.sgn_int) < 10)
	{
		Vector4cd out(inte(0), inte(1), inte(2), 0);
		return out;
	}

	if (abs(T.sgn_int) > 10)
	{
		if (T.sgn_int > 0)
		{
			Vector4cd out(inte(0), inte(1), inte(2), 0);
			return out;
		}
		if (T.sgn_int < 0)
		{
			if (T.n == 1)
			{
				con = con - inte(1)*a(T.x(1), T.y(1));
			}
			if (T.n == 2)
			{
				con = con - inte(2)*a(T.x(2), T.y(2));
			}
			Vector4cd out(inte(0), inte(1), inte(2), con);
			return out;
		}
	}
}

Vector3cd InhomogeneousThreeCavity::weak4(Vector2d &p1, Vector2d &p2, Vector2d &p3, double v1, double v2, Vector2d &p6)
{
	throw("wait to finish");
}

Vector3cd InhomogeneousThreeCavity::weak4(Triangle &T, double v2, double v3, bool beta)
{
	// beta:    true->beta     false->beta_

	Vector2d p1(T.x(0), T.y(0));
	Vector2d p2(T.x(1), T.y(1));
	Vector2d p3(T.x(2), T.y(2));
	Vector2d p6 = (p2 + p3) / 2;

	Matrix3d M;
	M << p1(0), p1(1), 1, p2(0), p2(1), 1, p3(0), p3(1), 1;

	Matrix3d invM = M.inverse(); //求M的逆

	MatrixXd MM(2, 3);
	MM = invM.block(0, 0, 2, 3);

	double x1 = p2(0);
	double y1 = p2(1);
	double x2 = p3(0);
	double y2 = p3(1);
	double s = (p3 - p2).norm();//s = norm(p3 - p2);

	double n1 = y2 - y1;
	double n2 = x1 - x2;
	Vector2d nb(n1 / s, n2 / s);
	Matrix2d b6;
	if (beta) {
		b6 = this->beta(p6(0), p6(1));
	}
	else {
		b6 = this->beta_(p6(0), p6(1));
	}

	Vector3cd out = -nb.transpose()*b6*MM*s*(v2 + v3) / 2;
	return out;
}

double InhomogeneousThreeCavity::weak4_int(Vector2d p1, Vector2d p2, Vector2d b, double v1, double v2)
{
	double x1 = p1(0);
	double y1 = p1(1);
	double x2 = p2(0);
	double y2 = p2(1);
	double s = (p2 - p1).norm();

	double n1 = y2 - y1;
	double n2 = x1 - x2;
	double nb = (n1*b(0) + n2*b(1)) / s;
	double out = nb*s*(v1 + v2) / 2;
	return out;
}

void InhomogeneousThreeCavity::weak5(double x1, double y1, double x2, double y2, double v1, double v2, vector<vector<double>> &nbound, vector<complex<double>> &out5m, complex<double> &out5)
{
	double s = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
	int nboundSize = nbound[0].size();
	int matrixGsize = nboundSize / 3;

	int separator_index_12, separator_index_23;
	separator_index_12 = nboundSize / this->cavity_number;         //第一个腔体口径面信息的结束索引，用于区分第一个和第二个腔体
	separator_index_23 = 2.0 * nboundSize / this->cavity_number;   //第二个腔体口径面信息的结束索引，用于区分第二个和第三个腔体


	if (x1 < separator1_2)
	{
		if (abs(nbound[1][0] - max(x1, x2)) < 1e-10)
		{
			MatrixXcd G_row(1, matrixGsize);
			G_row = G1_aperture.block(0, 0, 1, matrixGsize);
			VectorXcd temp5m1 = -G_row.transpose()*(v1 + v2)*s / 2;
			G_row = G2_aperture.block(0, 0, 1, matrixGsize);
			VectorXcd temp5m2 = -G_row.transpose()*(v1 + v2)*s / 2;
			G_row = G3_aperture.block(0, 0, 1, matrixGsize);
			VectorXcd temp5m3 = -G_row.transpose()*(v1 + v2)*s / 2;

			VectorXcd temp5m(nboundSize);
			temp5m << temp5m1, temp5m2, temp5m3;

			out5m = eigenVector2stdVector(temp5m);
			out5 = g_aperture(0)*(v1 + v2)*s * 0.5;
		}
		else if (abs(nbound[1][separator_index_12 - 1] - min(x1, x2)) < 1e-10)
		{
			MatrixXcd G_row(1, matrixGsize);
			G_row = G1_aperture.block(separator_index_12 - 1, 0, 1, matrixGsize);
			VectorXcd temp5m1 = -G_row.transpose()*(v1 + v2)*s / 2;
			G_row = G2_aperture.block(separator_index_12 - 1, 0, 1, matrixGsize);
			VectorXcd temp5m2 = -G_row.transpose()*(v1 + v2)*s / 2;
			G_row = G3_aperture.block(separator_index_12 - 1, 0, 1, matrixGsize);
			VectorXcd temp5m3 = -G_row.transpose()*(v1 + v2)*s / 2;

			VectorXcd temp5m(nboundSize);
			temp5m << temp5m1, temp5m2, temp5m3;

			out5m = eigenVector2stdVector(temp5m);
			out5 = g_aperture(separator_index_12 - 1)*(v1 + v2)*s * 0.5;
		}
		else
		{
			int sq = findXinNbound_multiple(nbound, min(x1, x2),0, separator_index_12 - 1);
			MatrixXcd G_row1 = G1_aperture.block(sq, 0, 1, matrixGsize);
			MatrixXcd G_row2 = G1_aperture.block(sq + 1, 0, 1, matrixGsize);
			VectorXcd temp5m1 = -(G_row1.transpose() + G_row2.transpose())*(v1 + v2)*s / 4;
			G_row1 = G2_aperture.block(sq, 0, 1, matrixGsize);
			G_row2 = G2_aperture.block(sq + 1, 0, 1, matrixGsize);
			VectorXcd temp5m2 = -(G_row1.transpose() + G_row2.transpose())*(v1 + v2)*s / 4;
			G_row1 = G3_aperture.block(sq, 0, 1, matrixGsize);
			G_row2 = G3_aperture.block(sq + 1, 0, 1, matrixGsize);
			VectorXcd temp5m3 = -(G_row1.transpose() + G_row2.transpose())*(v1 + v2)*s / 4;

			VectorXcd temp5m(nboundSize);
			temp5m << temp5m1, temp5m2, temp5m3;

			out5m = eigenVector2stdVector(temp5m);
			out5 = (g_aperture(sq) + g_aperture(sq + 1)) * (v1 + v2)*s *0.25;
		}
	}
	else if (x1 >= separator1_2 && x1 <= separator2_3)
	{
		if (abs(nbound[1][separator_index_12] - max(x1, x2)) < 1e-10)
		{
			MatrixXcd G_row(1, matrixGsize);
			G_row = G1_aperture.block(0, 0, 1, matrixGsize);
			VectorXcd temp5m1 = -G_row.transpose()*(v1 + v2)*s / 2;
			G_row = G2_aperture.block(0, 0, 1, matrixGsize);
			VectorXcd temp5m2 = -G_row.transpose()*(v1 + v2)*s / 2;
			G_row = G3_aperture.block(0, 0, 1, matrixGsize);
			VectorXcd temp5m3 = -G_row.transpose()*(v1 + v2)*s / 2;

			VectorXcd temp5m(nboundSize);
			temp5m << temp5m1, temp5m2, temp5m3;

			out5m = eigenVector2stdVector(temp5m);
			out5 = g_aperture(separator_index_12)*(v1 + v2)*s * 0.5;
		}
		else if (abs(nbound[1][separator_index_23 - 1] - min(x1, x2)) < 1e-10)
		{
			MatrixXcd G_row(1, matrixGsize);
			G_row = G1_aperture.block(matrixGsize - 1, 0, 1, matrixGsize);
			VectorXcd temp5m1 = -G_row.transpose()*(v1 + v2)*s / 2;
			G_row = G2_aperture.block(matrixGsize - 1, 0, 1, matrixGsize);
			VectorXcd temp5m2 = -G_row.transpose()*(v1 + v2)*s / 2;
			G_row = G3_aperture.block(matrixGsize - 1, 0, 1, matrixGsize);
			VectorXcd temp5m3 = -G_row.transpose()*(v1 + v2)*s / 2;

			VectorXcd temp5m(nboundSize);
			temp5m << temp5m1, temp5m2, temp5m3;

			out5m = eigenVector2stdVector(temp5m);
			out5 = g_aperture(separator_index_23 - 1)*(v1 + v2)*s * 0.5;
		}
		else
		{
			int sq = findXinNbound_multiple(nbound, min(x1, x2), separator_index_12, separator_index_23-1);
			MatrixXcd G_row1 = G1_aperture.block(sq, 0, 1, matrixGsize);
			MatrixXcd G_row2 = G1_aperture.block(sq + 1, 0, 1, matrixGsize);
			VectorXcd temp5m1 = -(G_row1.transpose() + G_row2.transpose())*(v1 + v2)*s / 4;
			G_row1 = G2_aperture.block(sq, 0, 1, matrixGsize);
			G_row2 = G2_aperture.block(sq + 1, 0, 1, matrixGsize);
			VectorXcd temp5m2 = -(G_row1.transpose() + G_row2.transpose())*(v1 + v2)*s / 4;
			G_row1 = G3_aperture.block(sq, 0, 1, matrixGsize);
			G_row2 = G3_aperture.block(sq + 1, 0, 1, matrixGsize);
			VectorXcd temp5m3 = -(G_row1.transpose() + G_row2.transpose())*(v1 + v2)*s / 4;

			VectorXcd temp5m(nboundSize);
			temp5m << temp5m1, temp5m2, temp5m3;

			out5m = eigenVector2stdVector(temp5m);
			out5 = (g_aperture(separator_index_12+sq) + g_aperture(separator_index_12+sq + 1)) * (v1 + v2)*s *0.25;
		}
	}
	else
	{
		if (abs(nbound[1][separator_index_23] - max(x1, x2)) < 1e-10)
		{
			MatrixXcd G_row(1, matrixGsize);
			G_row = G1_aperture.block(0, 0, 1, matrixGsize);
			VectorXcd temp5m1 = -G_row.transpose()*(v1 + v2)*s / 2;
			G_row = G2_aperture.block(0, 0, 1, matrixGsize);
			VectorXcd temp5m2 = -G_row.transpose()*(v1 + v2)*s / 2;
			G_row = G3_aperture.block(0, 0, 1, matrixGsize);
			VectorXcd temp5m3 = -G_row.transpose()*(v1 + v2)*s / 2;

			VectorXcd temp5m(nboundSize);
			temp5m << temp5m1, temp5m2, temp5m3;

			out5m = eigenVector2stdVector(temp5m);
			out5 = g_aperture(separator_index_23)*(v1 + v2)*s * 0.5;
		}
		else if (abs(nbound[1][nboundSize - 1] - min(x1, x2)) < 1e-10)
		{
			MatrixXcd G_row(1, matrixGsize);
			G_row = G1_aperture.block(matrixGsize - 1, 0, 1, matrixGsize);
			VectorXcd temp5m1 = -G_row.transpose()*(v1 + v2)*s / 2;
			G_row = G2_aperture.block(matrixGsize - 1, 0, 1, matrixGsize);
			VectorXcd temp5m2 = -G_row.transpose()*(v1 + v2)*s / 2;
			G_row = G3_aperture.block(matrixGsize - 1, 0, 1, matrixGsize);
			VectorXcd temp5m3 = -G_row.transpose()*(v1 + v2)*s / 2;

			VectorXcd temp5m(nboundSize);
			temp5m << temp5m1, temp5m2, temp5m3;

			out5m = eigenVector2stdVector(temp5m);
			out5 = g_aperture(nboundSize - 1)*(v1 + v2)*s * 0.5;
		}
		else
		{
			int sq = findXinNbound_multiple(nbound, min(x1, x2), separator_index_23, nboundSize - 1);
			MatrixXcd G_row1 = G1_aperture.block(sq, 0, 1, matrixGsize);
			MatrixXcd G_row2 = G1_aperture.block(sq + 1, 0, 1, matrixGsize);
			VectorXcd temp5m1 = -(G_row1.transpose() + G_row2.transpose())*(v1 + v2)*s / 4;
			G_row1 = G2_aperture.block(sq, 0, 1, matrixGsize);
			G_row2 = G2_aperture.block(sq + 1, 0, 1, matrixGsize);
			VectorXcd temp5m2 = -(G_row1.transpose() + G_row2.transpose())*(v1 + v2)*s / 4;
			G_row1 = G3_aperture.block(sq, 0, 1, matrixGsize);
			G_row2 = G3_aperture.block(sq + 1, 0, 1, matrixGsize);
			VectorXcd temp5m3 = -(G_row1.transpose() + G_row2.transpose())*(v1 + v2)*s / 4;

			VectorXcd temp5m(nboundSize);
			temp5m << temp5m1, temp5m2, temp5m3;

			out5m = eigenVector2stdVector(temp5m);
			out5 = (g_aperture(separator_index_23+sq) + g_aperture(separator_index_23+sq + 1)) * (v1 + v2)*s *0.25;
		}
	}

}


VectorXcd InhomogeneousThreeCavity::setRightHand(TriangleMesh &U, TriangleMesh &L, vector<int> &nu)
{
	int m = this->meshWidth;
	int n = this->meshHeight;

	int size = nu.size();
	complex<double> zero(0, 0);
	VectorXcd right = VectorXcd::Constant(size, zero);

	//func = @fInt_f; !!!!!!!!!!!!

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
			if (nu[id] > 0)
			{
				complex<double> t1, t2, t3, t4, t5, t6;
				if (l == n - 1)
				{
					t1 = 0;
					t2 = 0;
					t3 = 0;
					t4 = fInt_f(U.Get(j, l), Vector3d(1, 0, 0));
					t5 = fInt_f(L.Get(j + 1, l), Vector3d(0, 0, 1));
					t6 = fInt_f(U.Get(j + 1, l), Vector3d(0, 1, 0));
				}
				else
				{
					t1 = fInt_f(L.Get(j + 1, l + 1), Vector3d(1, 0, 0));
					t2 = fInt_f(U.Get(j, l + 1), Vector3d(0, 0, 1));
					t3 = fInt_f(L.Get(j, l + 1), Vector3d(0, 1, 0));
					t4 = fInt_f(U.Get(j, l), Vector3d(1, 0, 0));
					t5 = fInt_f(L.Get(j + 1, l), Vector3d(0, 0, 1));
					t6 = fInt_f(U.Get(j + 1, l), Vector3d(0, 1, 0));
				}
				right[id] = t1 + t2 + t3 + t4 + t5 + t6;
			}
			else
			{
				right[id] = 0;
			}//endif:  (nu[id] > 0)
		}
	}
	return right;
}

complex<double> InhomogeneousThreeCavity::fInt_f(Triangle &T, Vector3d v)
{
	complex<double> out;
	if (T.sgn == -10)
	{
		out = fInt_int(T, v);
		return out;
	}

	if (abs(T.sgn) < 10)
	{
		Triangle T1, T2, T3;
		T1.x = Vector3d(T.x1, T.x4, T.x5);
		T2.x = Vector3d(T.x2, T.x3, T.x4);
		T3.x = Vector3d(T.x3, T.x5, T.x4);
		T1.y = Vector3d(T.y1, T.y4, T.y5);
		T2.y = Vector3d(T.y2, T.y3, T.y4);
		T3.y = Vector3d(T.y3, T.y5, T.y4);

		int s = abs(T.sgn);

		double r1 = T.r1;
		double r2 = T.r2;

		int j = s;
		int k = j % 3 + 1;
		int l = k % 3 + 1;

		double vp1 = v(j - 1);
		double vp2 = v(k - 1);
		double vp3 = v(l - 1);
		double vp4 = vp1 + (vp2 - vp1)*r1;
		double vp5 = vp1 + (vp3 - vp1)*r2;

		Vector3d v1(vp1, vp4, vp5);
		Vector3d v2(vp2, vp3, vp4);
		Vector3d v3(vp3, vp5, vp4);

		if (T.sgn > 0)
		{
			reloadtri(T2, T);
			reloadtri(T3, T);
			out = fInt_int(T2, v2) + fInt_int(T3, v3); //tr(T2, f, v2) + tr(T3, f, v3);
		}
		else
		{
			reloadtri(T1, T);
			out = fInt_int(T1, v1); //tr(T1, f, v1);
		}
		return out;
	}

	if (abs(T.sgn) > 10)
	{
		Triangle T1, T2;
		T1.x = Vector3d(T.x1, T.x2, T.x5);
		T1.y = Vector3d(T.y1, T.y2, T.y5);
		T2.x = Vector3d(T.x3, T.x5, T.x2);
		T2.y = Vector3d(T.y3, T.y5, T.y2);

		double vp1 = v(T.v(0) - 1);
		double vp2 = v(T.v(1) - 1);
		double vp3 = v(T.v(2) - 1);
		double vp5 = vp1 + (vp3 - vp1)*T.r2;

		Vector3d v1(vp1, vp2, vp5);
		Vector3d v2(vp3, vp5, vp2);

		if (T.sgn > 0)
		{
			reloadtri(T2, T);
			out = fInt_int(T2, v2); //tr(T2, f, v2);
		}
		else
		{
			reloadtri(T1, T);
			out = fInt_int(T1, v1); //tr(T1, f, v1);
		}
		return out;
	}
}

complex<double> InhomogeneousThreeCavity::fInt_int(Triangle &T, Vector3d v)
{
	complex<double> out;

	if (T.sgn_int == 10)
	{
		out = tr(T, v, true);
		return out;
	}
	if (T.sgn_int == -10)
	{
		out = tr(T, v, false);
		return out;
	}

	if (abs(T.sgn_int) < 10)
	{
		Triangle T1, T2, T3;
		T1.x = Vector3d(T.x1_int, T.x4_int, T.x5_int);
		T2.x = Vector3d(T.x2_int, T.x3_int, T.x4_int);
		T3.x = Vector3d(T.x3_int, T.x5_int, T.x4_int);
		T1.y = Vector3d(T.y1_int, T.y4_int, T.y5_int);
		T2.y = Vector3d(T.y2_int, T.y3_int, T.y4_int);
		T3.y = Vector3d(T.y3_int, T.y5_int, T.y4_int);

		double r1 = T.r1_int;
		double r2 = T.r2_int;

		int s = abs(T.sgn_int);
		int j = s;
		int k = j % 3 + 1;
		int l = k % 3 + 1;

		double vp1 = v(j - 1);
		double vp2 = v(k - 1);
		double vp3 = v(l - 1);
		double vp4 = vp1 + (vp2 - vp1)*r1;
		double vp5 = vp1 + (vp3 - vp1)*r2;

		Vector3d v1(vp1, vp4, vp5);
		Vector3d v2(vp2, vp3, vp4);
		Vector3d v3(vp3, vp5, vp4);

		if (T.sgn_int > 0)
		{
			//f1 = f;
			//f2 = f_;
			//f3 = f_;
			out = tr(T1, v1, true) + tr(T2, v2, false) + tr(T3, v3, false);
		}
		else
		{
			//f1 = f_;
			//f2 = f;
			//f3 = f;
			out = tr(T1, v1, false) + tr(T2, v2, true) + tr(T3, v3, true);
		}
		return out;
	}

	if (abs(T.sgn_int) > 10)
	{
		Triangle T1, T2;
		T1.x = Vector3d(T.x1_int, T.x2_int, T.x5_int);
		T1.y = Vector3d(T.y1_int, T.y2_int, T.y5_int);
		T2.x = Vector3d(T.x3_int, T.x5_int, T.x2_int);
		T2.y = Vector3d(T.y3_int, T.y5_int, T.y2_int);

		double vp1 = v(T.v_int(0) - 1);
		double vp2 = v(T.v_int(1) - 1);
		double vp3 = v(T.v_int(2) - 1);
		double vp5 = vp1 + (vp3 - vp1)*T.r2_int;

		Vector3d v1(vp1, vp2, vp5);
		Vector3d v2(vp3, vp5, vp2);

		if (T.sgn_int > 0)
		{
			//f1 = f;
			//f2 = f_;
			out = tr(T1, v1, true) + tr(T2, v2, false);
		}
		else
		{
			//f1 = f_;
			//f2 = f;
			out = tr(T1, v1, false) + tr(T2, v2, true);
		}
		return out;
	}

}

complex<double> InhomogeneousThreeCavity::tr(Triangle &T, Vector3d v, bool f)
{
	//f:        true->f     false->f_

	double s = calculateTriangleArea(T.x, T.y);//求三角形面积
	double x4 = (T.x(1) + T.x(2)) / 2;
	double x5 = (T.x(2) + T.x(0)) / 2;
	double x6 = (T.x(0) + T.x(1)) / 2;
	double y4 = (T.y(1) + T.y(2)) / 2;
	double y5 = (T.y(2) + T.y(0)) / 2;
	double y6 = (T.y(0) + T.y(1)) / 2;

	complex<double> f4, f5, f6;
	if (f) {
		f4 = this->f(x4, y4);
		f5 = this->f(x5, y5);
		f6 = this->f(x6, y6);
	}
	else {
		f4 = this->f_(x4, y4);
		f5 = this->f_(x5, y5);
		f6 = this->f_(x6, y6);
	}

	double v4 = (v(1) + v(2)) / 2;
	double v5 = (v(2) + v(0)) / 2;
	double v6 = (v(0) + v(1)) / 2;

	complex<double> out = (f4*v4 + f5*v5 + f6*v6)*s *(1.0 / 3);
	return out;
}

SparseMatrix<complex<double>> InhomogeneousThreeCavity::setA(double nn, vector<int> &nu, vector<vector<double>> &nbound, vector<vector<gridCell>> &grid)
{
	int m = this->meshWidth;
	int n = this->meshHeight;

	double dx = this->stepX;
	double dy = this->stepY;

	double left = this->virtualBorderLeft;
	double bottom = this->virtualBorderBottom;

	int apertureNumber = nbound[0].size();

	int separator_index_12, separator_index_23;
	separator_index_12 = apertureNumber / this->cavity_number;         //第一个腔体口径面信息的结束索引，用于区分第一个和第二个腔体
	separator_index_23 = 2 * apertureNumber / this->cavity_number;   //第二个腔体口径面信息的结束索引，用于区分第二个和第三个腔体


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
					if (sq < separator_index_12)
					{
						for (int apeIndex = 0 ; apeIndex <separator_index_12; apeIndex++)
						{
							if (apeIndex<sq - 1 || apeIndex>sq + 1)//这个条件确保不会再次考虑已经组装的部分
							{
								*(pA++) = cur_nu;
								*(pA++) = nbound[0][apeIndex];
								*(pA++) = curGrid.tm[apeIndex].real();
								*(pA++) = curGrid.tm[apeIndex].imag();
							}
						}
						for (int apeIndex = separator_index_12; apeIndex < apertureNumber; apeIndex++)
						{
							*(pA++) = cur_nu;
							*(pA++) = nbound[0][apeIndex];
							*(pA++) = curGrid.tm[apeIndex].real();
							*(pA++) = curGrid.tm[apeIndex].imag();
						}
					}
					else if(sq >= separator_index_12 && sq < separator_index_23)
					{
						for (int apeIndex = separator_index_12 ; apeIndex< separator_index_23;apeIndex++)
						{
							if (apeIndex<sq - 1 || apeIndex>sq + 1) // 这个条件确保不会再次考虑已经组装的部分
							{
								*(pA++) = cur_nu;
								*(pA++) = nbound[0][apeIndex];
								*(pA++) = curGrid.tm[apeIndex].real();
								*(pA++) = curGrid.tm[apeIndex].imag();
							}
						}
						for (int apeIndex = 0 ; apeIndex< separator_index_12;apeIndex++)
						{
							*(pA++) = cur_nu;
							*(pA++) = nbound[0][apeIndex];
							*(pA++) = curGrid.tm[apeIndex].real();
							*(pA++) = curGrid.tm[apeIndex].imag();
						}
						for (int apeIndex = separator_index_23; apeIndex< apertureNumber; apeIndex++)
						{
							*(pA++) = cur_nu;
							*(pA++) = nbound[0][apeIndex];
							*(pA++) = curGrid.tm[apeIndex].real();
							*(pA++) = curGrid.tm[apeIndex].imag();
						}
					}
					else
					{
						for (int apeIndex = separator_index_23; apeIndex < apertureNumber; apeIndex++)
						{
							if (apeIndex<sq - 1 || apeIndex>sq + 1) // 这个条件确保不会再次考虑已经组装的部分
							{
								*(pA++) = cur_nu;
								*(pA++) = nbound[0][apeIndex];
								*(pA++) = curGrid.tm[apeIndex].real();
								*(pA++) = curGrid.tm[apeIndex].imag();
							}
						}
						for (int apeIndex = 0 ; apeIndex < separator_index_23; apeIndex++)
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

	//将三元组形式数据转换为系数矩阵形式
	SparseMatrix<complex<double>> matrixA = buildSparseMatrixA(A_result, nn, ArowNum);

	return matrixA;
}

VectorXcd InhomogeneousThreeCavity::setB(vector<vector<gridCell>> &grid, VectorXcd &rh, int nn, vector<int> &nu)
{
	int m = this->meshWidth;
	int n = this->meshHeight;

	VectorXcd B(nn);

	for (int j = 0; j < m - 1; j++)
	{
		for (int l = 0; l < n; l++)
		{
			int id = j + l*(m - 1);
			int cur_nu = nu[id];
			if (cur_nu > 0)
			{
				//gridCell curGrid = grid[j][l];
				B(cur_nu - 1) = rh(id) - grid[j][l]._const;
			}
		}
	}
	return B;
}


///待重写！！！！！！！！！！！！！！！！
void InhomogeneousThreeCavity::assign(VectorXcd solution, TriangleMesh &U, TriangleMesh &L, vector<int> &nu)
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
			Triangle *a = U.GetPt(j, l);
			Triangle *b = L.GetPt(j, l);

			if (a->sgn == -10)
			{
				complex<double> a1, a2, a3;
				if (abs(phi(a->x(0), a->y(0))) < 1e-10)
					a1 = this->u(a->x(0), a->y(0));
				else
					a1 = solution(nu[id] - 1);

				if (abs(phi(a->x(1), a->y(1))) < 1e-10)
					a2 = this->u(a->x(1), a->y(1));
				else
					a2 = solution(nu[id - 1] - 1);

				if (abs(phi(a->x(2), a->y(2))) < 1e-10)
					a3 = this->u(a->x(2), a->y(2));
				else
					a3 = solution(nu[id - m + 1] - 1);

				a->z = Vector3cd(a1, a2, a3);
			}

			if (b->sgn == -10)
			{
				complex<double> a1, a2, a3;
				if (abs(phi(b->x(0), b->y(0))) < 1e-10)
					a1 = this->u(b->x(0), b->y(0));
				else
					a1 = solution(nu[id - m] - 1);

				if (abs(phi(b->x(1), b->y(1))) < 1e-10)
					a2 = this->u(b->x(1), b->y(1));
				else
					a2 = solution(nu[id - m + 1] - 1);

				if (abs(phi(b->x(2), b->y(2))) < 1e-10)
					a3 = this->u(b->x(2), b->y(2));
				else
					a3 = solution(nu[id - 1] - 1);

				b->z = Vector3cd(a1, a2, a3);
			}

			if (abs(a->sgn) < 10)
			{
				int r = abs(a->sgn);
				int s = r % 3 + 1;
				int t = s % 3 + 1;
				if (a->sgn < 0)
				{
					if (r == 1)
						a->z1 = solution(nu[id] - 1);

					if (r == 2)
						a->z1 = solution(nu[id - 1] - 1);

					if (r == 3)
						a->z1 = solution(nu[id - m + 1] - 1);
				}
				else
				{
					if (a->sgn > 0)
					{
						if (s == 1)
							a->z2 = solution(nu[id] - 1);

						if (s == 2)
							a->z2 = solution(nu[id - 1] - 1);

						if (s == 3)
							a->z2 = solution(nu[id - m + 1] - 1);

						if (t == 1)
							a->z3 = solution(nu[id] - 1);

						if (t == 2)
							a->z3 = solution(nu[id - 1] - 1);

						if (t == 3)
							a->z3 = solution(nu[id - m + 1] - 1);
					}
				}
				a->z4 = this->u(a->x4, a->y4);
				a->z5 = this->u(a->x5, a->y5);
			}

			if (abs(b->sgn) < 10)
			{
				int r = abs(b->sgn);
				int s = r % 3 + 1;
				int t = s % 3 + 1;
				if (b->sgn < 0)
				{
					if (r == 1)
						b->z1 = solution(nu[id - m] - 1);

					if (r == 2)
						b->z1 = solution(nu[id - m + 1] - 1);

					if (r == 3)
						b->z1 = solution(nu[id - 1] - 1);
				}
				else
				{
					if (b->sgn > 0)
					{
						if (s == 1)
							b->z2 = solution(nu[id - m] - 1);

						if (s == 2)
							b->z2 = solution(nu[id - m + 1] - 1);

						if (s == 3)
							b->z2 = solution(nu[id - 1] - 1);

						if (t == 1)
							b->z3 = solution(nu[id - m] - 1);

						if (t == 2)
							b->z3 = solution(nu[id - m + 1] - 1);

						if (t == 3)
							b->z3 = solution(nu[id - 1] - 1);
					}
				}
				b->z4 = this->u(b->x4, b->y4);
				b->z5 = this->u(b->x5, b->y5);
			}

			if (abs(a->sgn) > 10)
			{
				if (a->sgn < 0)
				{
					if (a->v(0) == 1)
						a->z1 = solution(nu[id] - 1);

					if (a->v(0) == 2)
						a->z1 = solution(nu[id - 1] - 1);

					if (a->v(0) == 3)
						a->z1 = solution(nu[id - m + 1] - 1);
				}
				if (a->sgn > 0)
				{
					if (a->v(2) == 1)
						a->z3 = solution(nu[id] - 1);

					if (a->v(2) == 2)
						a->z3 = solution(nu[id - 1] - 1);

					if (a->v(2) == 3)
						a->z3 = solution(nu[id - m + 1] - 1);
				}
				a->z2 = this->u(a->x2, a->y2);
				a->z5 = this->u(a->p5(0), a->p5(1));
			}

			if (abs(b->sgn) > 10)
			{
				if (b->sgn < 0)
				{
					if (b->v(0) == 1)
						b->z1 = solution(nu[id - m] - 1);

					if (b->v(0) == 2)
						b->z1 = solution(nu[id - m + 1] - 1);

					if (b->v(0) == 3)
						b->z1 = solution(nu[id - 1] - 1);
				}
				if (b->sgn > 0)
				{
					if (b->v(2) == 1)
						b->z3 = solution(nu[id - m] - 1);

					if (b->v(2) == 2)
						b->z3 = solution(nu[id - m + 1] - 1);

					if (b->v(2) == 3)
						b->z3 = solution(nu[id - 1] - 1);
				}
				b->z2 = this->u(b->x2, b->y2);
				b->z5 = this->u(b->p5(0), b->p5(1));
			}

		}// end for (int l = 0; l < n; l++)
	}
}

///待重写！！！！！！！！！！！！！！！！
VectorXcd InhomogeneousThreeCavity::getAperture(vector<vector<double>> &nbound, TriangleMesh &U, TriangleMesh &L)
{
	int m = this->meshWidth;
	int n = this->meshHeight;
	int aperture_Num = nbound[0].size();
	VectorXcd ApertureValue(aperture_Num);

	int index = 0;
	for (int j = 0; j < m; j++)
	{
		// l == n -1
		Triangle tri = U.Get(j, n - 1);

		if (tri.x(0) == nbound[1][index])
		{
			ApertureValue(index) = tri.z(0);
			index++;
			if (index >= aperture_Num)
				break;
		}
	}

	return ApertureValue;
}

///待重写！！！！！！！！！！！！！！！！
void InhomogeneousThreeCavity::drawAperture(VectorXd &plotX, VectorXd &plotY, vector<vector<double>> &nbound, int sign)
{
	int plotArraySize = this->solutionOfAperture.size() + 2;
	plotX = VectorXd(plotArraySize);
	plotY = VectorXd(plotArraySize);

	//plotX(0) = this->apertureLeft;
	//plotX(plotArraySize - 1) = this->apertureRight;
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
		for (int i = 1; i < plotArraySize - 1; i++)
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








/*  --------------------------------------
--------------------------------------
private
--------------------------------------
--------------------------------------
*/

//***setTri中需要使用的子函数***

//判断并设置三角形类型，与界面交点位置等信息
Triangle InhomogeneousThreeCavity::tri(double a1, double a2, double a3, double x1, double x2, double x3, double y1, double y2, double y3)
{
	Triangle l;
	l.x = Vector3d(x1, x2, x3);
	l.y = Vector3d(y1, y2, y3);

	//-------------------------------------
	//set value for l
	//-------------------------------------
	if (abs(a1) < 1e-10 && abs(a2) < 1e-10 && abs(a3) < 1e-10)
	{
		//三角形三点都在界面上
		double ac = phi((x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3);
		if (ac > 1e-10)
			l.sgn = 10;
		if (ac < -1e-10)
			l.sgn = -10;
		return l;
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
		l.sgn = 10;
		return l;
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
		l.sgn = -10;
		return l;
	}

	if (a1 > 1e-10 && a2 < -1e-10 && a3 < -1e-10)
	{
		//三角形越过界面
		// (x1, y1)位于Ω + (x2, y2)、(x3, y3)位于Ω -
		l.sgn = 1;
		l.r1 = findzero(x1, y1, x2, y2);
		l.r2 = findzero(x1, y1, x3, y3);
		return l;
	}
	if (a1 < -1e-10 && a2>1e-10 && a3 > 1e-10)
	{
		// 三角形越过界面
		//a1位于Ω - a2、a3位于Ω +
		l.sgn = -1;
		l.r1 = findzero(x1, y1, x2, y2);
		l.r2 = findzero(x1, y1, x3, y3);
		return l;
	}
	if (a2 > 1e-10 && a1 < -1e-10 && a3 < -1e-10)
	{
		//三角形越过界面
		//a2位于Ω + a1、a3位于Ω -
		l.sgn = 2;
		l.r1 = findzero(x2, y2, x3, y3);
		l.r2 = findzero(x2, y2, x1, y1);
		return l;
	}
	if (a2<-1e-10 && a1>1e-10 && a3 > 1e-10)
	{
		l.sgn = -2;
		l.r1 = findzero(x2, y2, x3, y3);
		l.r2 = findzero(x2, y2, x1, y1);
		return l;
	}
	if (a3 > 1e-10 && a2 < -1e-10 && a1 < -1e-10)
	{
		l.sgn = 3;
		l.r1 = findzero(x3, y3, x1, y1);
		l.r2 = findzero(x3, y3, x2, y2);
		return l;
	}
	if (a3 < -1e-10 && a2>1e-10 && a1 > 1e-10)
	{
		l.sgn = -3;
		l.r1 = findzero(x3, y3, x1, y1);
		l.r2 = findzero(x3, y3, x2, y2);
		return l;
	}

	//一个点在界面上，一个在Ω + ，一个在Ω -
	if (abs(a2) < 1e-10 && a1 > 1e-10 && a3 < -1e-10)
	{
		// a2在界面上，a1在Ω + ，a3在Ω -
		l.sgn = 11;
		l.A = Vector2d(x1, y1);  //与B按照逆时针相对排列
		l.B = Vector2d(x2, y2);  //位于Γ上的点
		l.C = Vector2d(x3, y3);  //与B按照逆时针相对排列
		l.v = Vector3i(1, 2, 3); //标识ABC三点分别对应的点
		l.r1 = 1;
		l.r2 = findzero(x1, y1, x3, y3);
		return l;
	}
	if (abs(a2) < 1e-10 && a1 < -1e-10 && a3>1e-10)
	{
		l.sgn = -11;
		l.A = Vector2d(x1, y1);
		l.B = Vector2d(x2, y2);
		l.C = Vector2d(x3, y3);
		l.v = Vector3i(1, 2, 3);
		l.r1 = 1;
		l.r2 = findzero(x1, y1, x3, y3);
		return l;
	}
	if (abs(a3) < 1e-10 && a2 > 1e-10 && a1 < -1e-10)
	{
		l.sgn = 12;
		l.A = Vector2d(x2, y2);
		l.B = Vector2d(x3, y3);
		l.C = Vector2d(x1, y1);
		l.v = Vector3i(2, 3, 1);
		l.r1 = 1;
		l.r2 = findzero(x2, y2, x1, y1);
		return l;
	}
	if (abs(a3) < 1e-10 && a2 < -1e-10 && a1>1e-10)
	{
		l.sgn = -12;
		l.A = Vector2d(x2, y2);
		l.B = Vector2d(x3, y3);
		l.C = Vector2d(x1, y1);
		l.v = Vector3i(2, 3, 1);
		l.r1 = 1;
		l.r2 = findzero(x2, y2, x1, y1);
		return l;
	}
	if (abs(a1) < 1e-10 && a2 < -1e-10 && a3>1e-10)
	{
		l.sgn = 13;
		l.A = Vector2d(x3, y3);
		l.B = Vector2d(x1, y1);
		l.C = Vector2d(x2, y2);
		l.v = Vector3i(3, 1, 2);
		l.r1 = 1;
		l.r2 = findzero(x3, y3, x2, y2);
		return l;
	}
	if (abs(a1) < 1e-10 && a2 > 1e-10 && a3 < -1e-10)
	{
		l.sgn = -13;
		l.A = Vector2d(x3, y3);
		l.B = Vector2d(x1, y1);
		l.C = Vector2d(x2, y2);
		l.v = Vector3i(3, 1, 2);
		l.r1 = 1;
		l.r2 = findzero(x3, y3, x2, y2);
		return l;
	}
}

void InhomogeneousThreeCavity::tri_int(double a1, double a2, double a3, double x1, double  x2, double  x3, double y1, double y2, double y3, Triangle &l)
{
	l.x_int = Vector3d(x1, x2, x3);
	l.y_int = Vector3d(y1, y2, y3);

	//-------------------------------------
	//set value for l
	//-------------------------------------
	if (abs(a1) < 1e-10 && abs(a2) < 1e-10 && abs(a3) < 1e-10)
	{
		//三角形三点都在界面上
		double ac = phi_int((x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3);
		if (ac > 1e-10)
			l.sgn_int = 10;
		if (ac < -1e-10)
			l.sgn_int = -10;
		return;
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
		l.sgn_int = 10;
		return;
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
		l.sgn_int = -10;
		return;
	}

	if (a1 > 1e-10 && a2 < -1e-10 && a3 < -1e-10)
	{
		//三角形越过界面
		// (x1, y1)位于Ω + (x2, y2)、(x3, y3)位于Ω -
		l.sgn_int = 1;
		l.r1_int = findzero_int(x1, y1, x2, y2);
		l.r2_int = findzero_int(x1, y1, x3, y3);
		return;
	}
	if (a1 < -1e-10 && a2>1e-10 && a3 > 1e-10)
	{
		// 三角形越过界面
		//a1位于Ω - a2、a3位于Ω +
		l.sgn_int = -1;
		l.r1_int = findzero_int(x1, y1, x2, y2);
		l.r2_int = findzero_int(x1, y1, x3, y3);
		return;
	}
	if (a2 > 1e-10 && a1 < -1e-10 && a3 < -1e-10)
	{
		//三角形越过界面
		//a2位于Ω + a1、a3位于Ω -
		l.sgn_int = 2;
		l.r1_int = findzero_int(x2, y2, x3, y3);
		l.r2_int = findzero_int(x2, y2, x1, y1);
		return;
	}
	if (a2<-1e-10 && a1>1e-10 && a3 > 1e-10)
	{
		l.sgn_int = -2;
		l.r1_int = findzero_int(x2, y2, x3, y3);
		l.r2_int = findzero_int(x2, y2, x1, y1);
		return;
	}
	if (a3 > 1e-10 && a2 < -1e-10 && a1 < -1e-10)
	{
		l.sgn_int = 3;
		l.r1_int = findzero_int(x3, y3, x1, y1);
		l.r2_int = findzero_int(x3, y3, x2, y2);
		return;
	}
	if (a3 < -1e-10 && a2>1e-10 && a1 > 1e-10)
	{
		l.sgn_int = -3;
		l.r1_int = findzero_int(x3, y3, x1, y1);
		l.r2_int = findzero_int(x3, y3, x2, y2);
		return;
	}

	//一个点在界面上，一个在Ω + ，一个在Ω -
	if (abs(a2) < 1e-10 && a1 > 1e-10 && a3 < -1e-10)
	{
		// a2在界面上，a1在Ω + ，a3在Ω -
		l.sgn_int = 11;
		l.A_int = Vector2d(x1, y1);  //与B按照逆时针相对排列
		l.B_int = Vector2d(x2, y2);  //位于Γ上的点
		l.C_int = Vector2d(x3, y3);  //与B按照逆时针相对排列
		l.v_int = Vector3i(1, 2, 3); //标识ABC三点分别对应的点
		l.r1_int = 1;
		l.r2_int = findzero_int(x1, y1, x3, y3);
		return;
	}
	if (abs(a2) < 1e-10 && a1 < -1e-10 && a3>1e-10)
	{
		l.sgn_int = -11;
		l.A_int = Vector2d(x1, y1);
		l.B_int = Vector2d(x2, y2);
		l.C_int = Vector2d(x3, y3);
		l.v_int = Vector3i(1, 2, 3);
		l.r1_int = 1;
		l.r2_int = findzero_int(x1, y1, x3, y3);
		return;
	}
	if (abs(a3) < 1e-10 && a2 > 1e-10 && a1 < -1e-10)
	{
		l.sgn_int = 12;
		l.A_int = Vector2d(x2, y2);
		l.B_int = Vector2d(x3, y3);
		l.C_int = Vector2d(x1, y1);
		l.v_int = Vector3i(2, 3, 1);
		l.r1_int = 1;
		l.r2_int = findzero_int(x2, y2, x1, y1);
		return;
	}
	if (abs(a3) < 1e-10 && a2 < -1e-10 && a1>1e-10)
	{
		l.sgn_int = -12;
		l.A_int = Vector2d(x2, y2);
		l.B_int = Vector2d(x3, y3);
		l.C_int = Vector2d(x1, y1);
		l.v_int = Vector3i(2, 3, 1);
		l.r1_int = 1;
		l.r2_int = findzero_int(x2, y2, x1, y1);
		return;
	}
	if (abs(a1) < 1e-10 && a2 < -1e-10 && a3>1e-10)
	{
		l.sgn_int = 13;
		l.A_int = Vector2d(x3, y3);
		l.B_int = Vector2d(x1, y1);
		l.C_int = Vector2d(x2, y2);
		l.v_int = Vector3i(3, 1, 2);
		l.r1_int = 1;
		l.r2_int = findzero_int(x3, y3, x2, y2);
		return;
	}
	if (abs(a1) < 1e-10 && a2 > 1e-10 && a3 < -1e-10)
	{
		l.sgn_int = -13;
		l.A_int = Vector2d(x3, y3);
		l.B_int = Vector2d(x1, y1);
		l.C_int = Vector2d(x2, y2);
		l.v_int = Vector3i(3, 1, 2);
		l.r1_int = 1;
		l.r2_int = findzero_int(x3, y3, x2, y2);
		return;
	}
}

void InhomogeneousThreeCavity::tri_int1(double a1, double a2, double a3, double x1, double  x2, double  x3, double y1, double y2, double y3, Triangle &l)
{
	l.x_intex = Vector3d(x1, x2, x3);
	l.y_intex = Vector3d(y1, y2, y3);

	//-------------------------------------
	//set value for l
	//-------------------------------------
	if (abs(a1) < 1e-10 && abs(a2) < 1e-10 && abs(a3) < 1e-10)
	{
		//三角形三点都在界面上
		double ac = phi_int((x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3);
		if (ac > 1e-10)
			l.sgn_intex = 10;
		if (ac < -1e-10)
			l.sgn_intex = -10;
		return;
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
		l.sgn_intex = 10;
		return;
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
		l.sgn_intex = -10;
		return;
	}

	if (a1 > 1e-10 && a2 < -1e-10 && a3 < -1e-10)
	{
		//三角形越过界面
		// (x1, y1)位于Ω + (x2, y2)、(x3, y3)位于Ω -
		l.sgn_intex = 1;
		l.r1_intex = findzero_int(x1, y1, x2, y2);
		l.r2_intex = findzero_int(x1, y1, x3, y3);
		return;
	}
	if (a1 < -1e-10 && a2>1e-10 && a3 > 1e-10)
	{
		// 三角形越过界面
		//a1位于Ω - a2、a3位于Ω +
		l.sgn_intex = -1;
		l.r1_intex = findzero_int(x1, y1, x2, y2);
		l.r2_intex = findzero_int(x1, y1, x3, y3);
		return;
	}
	if (a2 > 1e-10 && a1 < -1e-10 && a3 < -1e-10)
	{
		//三角形越过界面
		//a2位于Ω + a1、a3位于Ω -
		l.sgn_intex = 2;
		l.r1_intex = findzero_int(x2, y2, x3, y3);
		l.r2_intex = findzero_int(x2, y2, x1, y1);
		return;
	}
	if (a2<-1e-10 && a1>1e-10 && a3 > 1e-10)
	{
		l.sgn_intex = -2;
		l.r1_intex = findzero_int(x2, y2, x3, y3);
		l.r2_intex = findzero_int(x2, y2, x1, y1);
		return;
	}
	if (a3 > 1e-10 && a2 < -1e-10 && a1 < -1e-10)
	{
		l.sgn_intex = 3;
		l.r1_intex = findzero_int(x3, y3, x1, y1);
		l.r2_intex = findzero_int(x3, y3, x2, y2);
		return;
	}
	if (a3 < -1e-10 && a2>1e-10 && a1 > 1e-10)
	{
		l.sgn_intex = -3;
		l.r1_intex = findzero_int(x3, y3, x1, y1);
		l.r2_intex = findzero_int(x3, y3, x2, y2);
		return;
	}

	//一个点在界面上，一个在Ω + ，一个在Ω -
	if (abs(a2) < 1e-10 && a1 > 1e-10 && a3 < -1e-10)
	{
		// a2在界面上，a1在Ω + ，a3在Ω -
		l.sgn_intex = 11;
		l.A_intex = Vector2d(x1, y1);  //与B按照逆时针相对排列
		l.B_intex = Vector2d(x2, y2);  //位于Γ上的点
		l.C_intex = Vector2d(x3, y3);  //与B按照逆时针相对排列
		l.v_intex = Vector3i(1, 2, 3); //标识ABC三点分别对应的点
		l.r1_intex = 1;
		l.r2_intex = findzero_int(x1, y1, x3, y3);
		return;
	}
	if (abs(a2) < 1e-10 && a1 < -1e-10 && a3>1e-10)
	{
		l.sgn_intex = -11;
		l.A_intex = Vector2d(x1, y1);
		l.B_intex = Vector2d(x2, y2);
		l.C_intex = Vector2d(x3, y3);
		l.v_intex = Vector3i(1, 2, 3);
		l.r1_intex = 1;
		l.r2_intex = findzero_int(x1, y1, x3, y3);
		return;
	}
	if (abs(a3) < 1e-10 && a2 > 1e-10 && a1 < -1e-10)
	{
		l.sgn_intex = 12;
		l.A_intex = Vector2d(x2, y2);
		l.B_intex = Vector2d(x3, y3);
		l.C_intex = Vector2d(x1, y1);
		l.v_intex = Vector3i(2, 3, 1);
		l.r1_intex = 1;
		l.r2_intex = findzero_int(x2, y2, x1, y1);
		return;
	}
	if (abs(a3) < 1e-10 && a2 < -1e-10 && a1>1e-10)
	{
		l.sgn_intex = -12;
		l.A_intex = Vector2d(x2, y2);
		l.B_intex = Vector2d(x3, y3);
		l.C_intex = Vector2d(x1, y1);
		l.v_intex = Vector3i(2, 3, 1);
		l.r1_intex = 1;
		l.r2_intex = findzero_int(x2, y2, x1, y1);
		return;
	}
	if (abs(a1) < 1e-10 && a2 < -1e-10 && a3>1e-10)
	{
		l.sgn_intex = 13;
		l.A_intex = Vector2d(x3, y3);
		l.B_intex = Vector2d(x1, y1);
		l.C_intex = Vector2d(x2, y2);
		l.v_intex = Vector3i(3, 1, 2);
		l.r1_intex = 1;
		l.r2_intex = findzero_int(x3, y3, x2, y2);
		return;
	}
	if (abs(a1) < 1e-10 && a2 > 1e-10 && a3 < -1e-10)
	{
		l.sgn_intex = -13;
		l.A_intex = Vector2d(x3, y3);
		l.B_intex = Vector2d(x1, y1);
		l.C_intex = Vector2d(x2, y2);
		l.v_intex = Vector3i(3, 1, 2);
		l.r1_intex = 1;
		l.r2_intex = findzero_int(x3, y3, x2, y2);
		return;
	}
}

//按照统一规则标准化三角形坐标信息
void InhomogeneousThreeCavity::solveTri(Triangle &l)
{
	// 按照统一规则标准化三角形坐标信息

	// sgn == 10时本函数不需单独处理
	// sgn == -10时的情况不需要考虑

	if (abs(l.sgn) < 10)
	{
		int p1 = abs(l.sgn);
		int p2 = p1 % 3 + 1;
		int p3 = p2 % 3 + 1;
		l.x1 = l.x(p1 - 1); l.y1 = l.y(p1 - 1);
		l.x2 = l.x(p2 - 1); l.y2 = l.y(p2 - 1);
		l.x3 = l.x(p3 - 1); l.y3 = l.y(p3 - 1);
		l.x4 = l.x1 + (l.x2 - l.x1)*l.r1; l.y4 = l.y1 + (l.y2 - l.y1)*l.r1;
		l.x5 = l.x1 + (l.x3 - l.x1)*l.r2; l.y5 = l.y1 + (l.y3 - l.y1)*l.r2;
		return;
	}
	if (abs(l.sgn) > 10)
	{
		l.p5 = l.A + (l.C - l.A)*l.r2;
		l.x5 = l.p5(0); l.y5 = l.p5(1);
		l.x1 = l.A(0); l.y1 = l.A(1);
		l.x2 = l.B(0); l.y2 = l.B(1);
		l.x3 = l.C(0); l.y3 = l.C(1);
		return;
	}
}

void InhomogeneousThreeCavity::solveTri_int(Triangle &l)
{
	if (abs(l.sgn_int) < 10)
	{
		int p1 = abs(l.sgn_int);
		int p2 = p1 % 3 + 1;
		int p3 = p2 % 3 + 1;
		l.x1_int = l.x_int(p1 - 1); l.y1_int = l.y_int(p1 - 1);
		l.x2_int = l.x_int(p2 - 1); l.y2_int = l.y_int(p2 - 1);
		l.x3_int = l.x_int(p3 - 1); l.y3_int = l.y_int(p3 - 1);
		l.x4_int = l.x1_int + (l.x2_int - l.x1_int)*l.r1_int; l.y4_int = l.y1_int + (l.y2_int - l.y1_int)*l.r1_int;
		l.x5_int = l.x1_int + (l.x3_int - l.x1_int)*l.r2_int; l.y5_int = l.y1_int + (l.y3_int - l.y1_int)*l.r2_int;
		Vector4d u4, u4_, u5, u5_;
		ustar_int(l, u4, u5, u4_, u5_);
		l.u4_int = u4;
		l.u4_int_ = u4_;
		l.u5_int = u5;
		l.u5_int_ = u5_;
		return;
	}


	if (abs(l.sgn_int) > 10)
	{
		Vector2d p_int5 = l.A_int + (l.C_int - l.A_int)*l.r2_int;
		l.x5_int = p_int5(0);
		l.y5_int = p_int5(1);
		l.x1_int = l.A_int(0);
		l.y1_int = l.A_int(1);
		l.x2_int = l.B_int(0);
		l.y2_int = l.B_int(1);
		l.x3_int = l.C_int(0);
		l.y3_int = l.C_int(1);
		Vector4d u5, u5_;
		ustar2_int(l, u5, u5_);
		l.u5_int = u5;
		l.u5_int_ = u5_;
		return;
	}
}


void InhomogeneousThreeCavity::solveTri_int1(Triangle &l)
{
	if (abs(l.sgn_intex) < 10)
	{
		int p1 = abs(l.sgn_intex);
		int p2 = p1 % 3 + 1;
		int p3 = p2 % 3 + 1;
		l.x1_intex = l.x_intex(p1 - 1); l.y1_intex = l.y_intex(p1 - 1);
		l.x2_intex = l.x_intex(p2 - 1); l.y2_intex = l.y_intex(p2 - 1);
		l.x3_intex = l.x_intex(p3 - 1); l.y3_intex = l.y_intex(p3 - 1);
		l.x4_intex = l.x1_intex + (l.x2_intex - l.x1_intex)*l.r1_intex; l.y4_intex = l.y1_intex + (l.y2_intex - l.y1_intex)*l.r1_intex;
		l.x5_intex = l.x1_intex + (l.x3_intex - l.x1_intex)*l.r2_intex; l.y5_intex = l.y1_intex + (l.y3_intex - l.y1_intex)*l.r2_intex;
		Vector4d u4, u4_, u5, u5_;
		ustar_int1(l, u4, u5, u4_, u5_);
		l.u4_intex = u4;
		l.u4_intex_ = u4_;
		l.u5_intex = u5;
		l.u5_intex_ = u5_;
		return;
	}

	if (abs(l.sgn_intex) > 10)
	{
		Vector2d p_intex5 = l.A_intex + (l.C_intex - l.A_intex)*l.r2_intex;
		l.x5_intex = p_intex5(0);
		l.y5_intex = p_intex5(1);
		l.x1_intex = l.A_intex(0);
		l.y1_intex = l.A_intex(1);
		l.x2_intex = l.B_intex(0);
		l.y2_intex = l.B_intex(1);
		l.x3_intex = l.C_intex(0);
		l.y3_intex = l.C_intex(1);
		Vector4d u5, u5_;
		ustar2_int1(l, u5, u5_);
		l.u5_intex = u5;
		l.u5_intex_ = u5_;
		return;
	}
}

double InhomogeneousThreeCavity::findzero_int(double x1, double y1, double x2, double y2)
{
	// (x1,y1)与(x2,y2)位于interface的两侧，故假设interfaceΓ交线段(x1,y1)、(x2,y2)于点cut处
	// 由于Γ上的点水平集函数值为0
	// findzero函数使用二分法寻找点cut，直至满足精度要求 | phi(cut) | <eps
	// out返回交点的对应位置，从(x1,y1)->(x2,y2)分别对应0->1

	double eps = 1e-10;

	double left_x = x1;
	double left_y = y1;
	double right_x = x2;
	double right_y = y2;
	double left_phi = phi_int(x1, y1);    //左端点水平集
	double right_phi = phi_int(x2, y2);   //右端点水平集

										  // 初始化
	double cut_old_x = (left_x + right_x)*0.5;       // 初始化起始分割点cut为(x1,y1)、(x2,y2)中点
	double cut_old_y = (left_y + right_y)*0.5;
	double cut_phi = phi_int(cut_old_x, cut_old_y);      // 起始分割点cut的水平集函数值

														 // 迭代直至满足精度要求
	double cut_new_x, cut_new_y;
	while (abs(cut_phi) >= eps)
	{
		//不满足精度要求，继续二分寻找

		// 首先判断真正的交点位于哪一边？left——cut_old or cut_old——right
		if (left_phi * cut_phi < 0)
		{
			// 交点位于 left——cut_old
			cut_new_x = (cut_old_x + left_x)*0.5;       // 新的中点
			cut_new_y = (cut_old_y + left_y)*0.5;
			right_x = cut_old_x;                        // 修改右端点
			right_y = cut_old_y;
			right_phi = phi_int(right_x, right_y);          // 修改右端点水平集
		}
		else if (right_phi * cut_phi < 0)
		{
			// 交点位于cut_old——right_x
			cut_new_x = (cut_old_x + right_x)*0.5;     // 新的中点
			cut_new_y = (cut_old_y + right_y)*0.5;
			left_x = cut_old_x;                        // 修改右端点
			left_y = cut_old_y;
			left_phi = phi_int(left_x, left_y);            // 修改右端点水平集
		}
		else
		{
			throw "something wrong in findzero.m";
		}
		cut_old_x = cut_new_x;
		cut_old_y = cut_new_y;
		cut_phi = phi_int(cut_old_x, cut_old_y);
	}
	double d_old_1 = sqrt(pow(cut_old_x - x1, 2) + pow(cut_old_y - y1, 2));
	double d_1_2 = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));

	return d_old_1 / d_1_2;
}


void InhomogeneousThreeCavity::ustar_int(Triangle &T, Vector4d &u4, Vector4d &u5, Vector4d &u4_, Vector4d &u5_)
{
	double x1 = T.x1_int;
	double y1 = T.y1_int;
	double x2 = T.x2_int;
	double y2 = T.y2_int;
	double x3 = T.x3_int;
	double y3 = T.y3_int;
	double x4 = T.x4_int;
	double y4 = T.y4_int;
	double x5 = T.x5_int;
	double y5 = T.y5_int;
	double x6 = (x4 + x5) / 2;
	double y6 = (y4 + y5) / 2;

	double a4 = a(x4, y4)*sign(T.sgn_int);
	double a5 = a(x5, y5)*sign(T.sgn_int);
	Vector2d bvalue = b(x6, y6);
	double b1 = bvalue(0);
	double b2 = bvalue(1);
	double s = sqrt((x4 - x5)*(x4 - x5) + (y4 - y5) *(y4 - y5));
	double n1 = (y4 - y5) / s*sign(T.sgn_int);
	double n2 = -(x4 - x5) / s*sign(T.sgn_int);

	Matrix3d  M2;
	M2 << x2, y2, 1,
		x3, y3, 1,
		x4, y4, 1;
	M2 = M2.inverse().eval();

	Vector3d u234 = Vector3d(x5, y5, 1).transpose() * M2;      // u2, u3, u4     

	Matrix3d  M1;
	M1 << x1, y1, 1,
		x4, y4, 1,
		x5, y5, 1;
	M1 = M1.inverse().eval();

	MatrixXd M1_12row(2, 3);
	M1_12row = M1.block(0, 0, 2, 3);
	MatrixXd M2_12row(2, 3);
	M2_12row = M2.block(0, 0, 2, 3);

	Vector3d A, B;
	if (T.sgn_int > 0) {
		A = Vector2d(n1, n2).transpose()*beta(x6, y6)*M1_12row;
		B = Vector2d(n1, n2).transpose() * beta_(x6, y6)*M2_12row;
	}
	else {
		A = Vector2d(n1, n2).transpose()*beta_(x6, y6)*M1_12row;
		B = Vector2d(n1, n2).transpose() * beta(x6, y6)*M2_12row;
	}
	double num = (b1*n1 + b2*n2)*sign(T.sgn_int);

	//u4:   u1, u2, u3, c
	u4(0) = A(0);
	u4(1) = A(2)*u234(0) - B(0);
	u4(2) = A(2)*u234(1) - B(1);
	u4(3) = A(1)*a4 + A(2)*a5 - num;
	u4 = u4 / (-1 * (A(1) - B(2) + A(2)*u234(2)));

	// u5:  u1, u2, u3, c
	u5(0) = u234(2)*u4(0);
	u5(1) = u234(0) + u234(2)*u4(1);
	u5(2) = u234(1) + u234(2)*u4(2);
	u5(3) = u234(2)*u4(3);

	u4_ = u4;
	u5_ = u5;

	u4(3) += a4;
	u5(3) += a5;

}

void InhomogeneousThreeCavity::ustar2_int(Triangle &T, Vector4d &u5, Vector4d &u5_)
{
	double x1 = T.x1_int;
	double y1 = T.y1_int;
	double x2 = T.x2_int;
	double y2 = T.y2_int;
	double x3 = T.x3_int;
	double y3 = T.y3_int;
	double x5 = T.x5_int;
	double y5 = T.y5_int;
	double x6 = (x2 + x5) / 2;
	double y6 = (y2 + y5) / 2;

	double a2 = a(x2, y2);
	double a5 = a(x5, y5);
	Vector2d bvalue = b(x6, y6);
	double b1 = bvalue(0);
	double b2 = bvalue(1);
	double s = sqrt((x2 - x5)*(x2 - x5) + (y2 - y5) *(y2 - y5));
	double n1 = (y2 - y5) / s;
	double n2 = -(x2 - x5) / s;
	double num = (b1*n1 + b2*n2);

	Matrix3d  M1;
	M1 << x1, y1, 1,
		x2, y2, 1,
		x5, y5, 1;
	M1 = M1.inverse();

	Matrix3d  M2;
	M2 << x2, y2, 1,
		x3, y3, 1,
		x5, y5, 1;
	M2 = M2.inverse();

	MatrixXd M1_12row(2, 3);
	M1_12row = M1.block(0, 0, 2, 3);
	MatrixXd M2_12row(2, 3);
	M2_12row = M2.block(0, 0, 2, 3);

	Vector3d A, B;
	if (T.sgn_int > 0)
	{
		A = Vector2d(n1, n2).transpose() * beta(x6, y6)*M1_12row;
		B = Vector2d(n1, n2).transpose() * beta_(x6, y6)*M2_12row;

		u5(0) = A(0);
		u5(1) = A(1) - B(0);
		u5(2) = -B(1);
		u5(3) = B(0)*a2 + B(2)*a5 - num;
		u5 = u5 / (B(2) - A(2));

		u5_ = u5;
		u5_(3) -= a5;
		return;
	}
	if (T.sgn_int < 0)
	{
		A = Vector2d(n1, n2).transpose() * beta_(x6, y6)*M1_12row;
		B = Vector2d(n1, n2).transpose() * beta(x6, y6)*M2_12row;

		u5(0) = -A(0);
		u5(1) = -A(1) + B(0);
		u5(2) = B(1);
		u5(3) = A(1)*a2 + A(2)*a5 - num;
		u5 = u5 / (-B(2) + A(2));

		u5_ = u5;
		u5_(3) -= a5;
		return;
	}
}

void InhomogeneousThreeCavity::ustar_int1(Triangle &T, Vector4d &u4, Vector4d &u5, Vector4d &u4_, Vector4d &u5_)
{
	double x1 = T.x1_intex;
	double y1 = T.y1_intex;
	double x2 = T.x2_intex;
	double y2 = T.y2_intex;
	double x3 = T.x3_intex;
	double y3 = T.y3_intex;
	double x4 = T.x4_intex;
	double y4 = T.y4_intex;
	double x5 = T.x5_intex;
	double y5 = T.y5_intex;
	double x6 = (x4 + x5) / 2;
	double y6 = (y4 + y5) / 2;

	double a4 = a(x4, y4)*sign(T.sgn_intex);
	double a5 = a(x5, y5)*sign(T.sgn_intex);
	Vector2d bvalue = b(x6, y6);
	double b1 = bvalue(0);
	double b2 = bvalue(1);
	double s = sqrt((x4 - x5)*(x4 - x5) + (y4 - y5) *(y4 - y5));
	double n1 = (y4 - y5) / s*sign(T.sgn_intex);
	double n2 = -(x4 - x5) / s*sign(T.sgn_intex);

	Matrix3d  M2;
	M2 << x2, y2, 1,
		x3, y3, 1,
		x4, y4, 1;
	M2 = M2.inverse();

	Vector3d u234 = Vector3d(x5, y5, 1).transpose() * M2;

	Matrix3d  M1;
	M1 << x1, y1, 1,
		x4, y4, 1,
		x5, y5, 1;
	M1 = M1.inverse();

	MatrixXd M1_12row(2, 3);
	M1_12row = M1.block(0, 0, 2, 3);
	MatrixXd M2_12row(2, 3);
	M2_12row = M2.block(0, 0, 2, 3);

	Vector3d A, B;
	if (T.sgn_intex > 0) {
		A = Vector2d(n1, n2).transpose()*beta(x6, y6)*M1_12row;
		B = Vector2d(n1, n2).transpose() * beta_(x6, y6)*M2_12row;
	}
	else {
		A = Vector2d(n1, n2).transpose()*beta_(x6, y6)*M1_12row;
		B = Vector2d(n1, n2).transpose() * beta(x6, y6)*M2_12row;
	}
	double num = (b1*n1 + b2*n2)*sign(T.sgn_intex);

	//u4:   u1, u2, u3, c
	u4(0) = A(0);
	u4(1) = A(2)*u234(0) - B(0);
	u4(2) = A(2)*u234(1) - B(1);
	u4(3) = A(1)*a4 + A(2)*a5 - num;
	u4 = u4 / (-1 * (A(1) - B(2) + A(2)*u234(2)));

	// u5:  u1, u2, u3, c
	u5(0) = u234(2)*u4(0);
	u5(1) = u234(0) + u234(2)*u4(1);
	u5(2) = u234(1) + u234(2)*u4(2);
	u5(3) = u234(2)*u4(3);

	u4_ = u4;
	u5_ = u5;

	u4(3) += a4;
	u5(3) += a5;
}

void InhomogeneousThreeCavity::ustar2_int1(Triangle &T, Vector4d &u5, Vector4d &u5_)
{
	double x1 = T.x1_intex;
	double y1 = T.y1_intex;
	double x2 = T.x2_intex;
	double y2 = T.y2_intex;
	double x3 = T.x3_intex;
	double y3 = T.y3_intex;
	double x5 = T.x5_intex;
	double y5 = T.y5_intex;
	double x6 = (x2 + x5) / 2;
	double y6 = (y2 + y5) / 2;

	double a2 = a(x2, y2);
	double a5 = a(x5, y5);
	Vector2d bvalue = b(x6, y6);
	double b1 = bvalue(0);
	double b2 = bvalue(1);
	double s = sqrt((x2 - x5)*(x2 - x5) + (y2 - y5) *(y2 - y5));
	double n1 = (y2 - y5) / s;
	double n2 = -(x2 - x5) / s;
	double num = (b1*n1 + b2*n2);

	Matrix3d  M1;
	M1 << x1, y1, 1,
		x2, y2, 1,
		x5, y5, 1;
	M1 = M1.inverse();

	Matrix3d  M2;
	M2 << x2, y2, 1,
		x3, y3, 1,
		x5, y5, 1;
	M2 = M2.inverse();

	MatrixXd M1_12row(2, 3);
	M1_12row = M1.block(0, 0, 2, 3);
	MatrixXd M2_12row(2, 3);
	M2_12row = M2.block(0, 0, 2, 3);

	Vector3d A, B;
	if (T.sgn_intex > 0)
	{
		A = Vector2d(n1, n2).transpose() * beta(x6, y6)*M1_12row;
		B = Vector2d(n1, n2).transpose() * beta_(x6, y6)*M2_12row;

		u5(0) = A(0);
		u5(1) = A(1) - B(0);
		u5(2) = -B(1);
		u5(3) = B(0)*a2 + B(2)*a5 - num;
		u5 = u5 / (B(2) - A(2));

		u5_ = u5;
		u5_(3) -= a5;
		return;
	}
	if (T.sgn_intex < 0)
	{
		A = Vector2d(n1, n2).transpose() * beta_(x6, y6)*M1_12row;
		B = Vector2d(n1, n2).transpose() * beta(x6, y6)*M2_12row;

		u5(0) = -A(0);
		u5(1) = -A(1) + B(0);
		u5(2) = B(1);
		u5(3) = A(1)*a2 + A(2)*a5 - num;
		u5 = u5 / (-B(2) + A(2));

		u5_ = u5;
		u5_(3) -= a5;
		return;
	}
}


//***setGrid中需要使用的子函数***

complex<double> InhomogeneousThreeCavity::value(double x, double y)
{
	if (phi_int(x, y) > pow(10, -10) || abs(phi_int(x, y)) < pow(10, -10))
	{
		return u(x, y);
	}

	if (phi_int(x, y) < -pow(10, -10))
	{
		return u_(x, y);
	}
}


void InhomogeneousThreeCavity::reloadtri(Triangle &T2, Triangle &T)
{
	T2.sgn_int = T.sgn_int;
	T2.x_int = T.x_int;
	T2.y_int = T.y_int;
	if (abs(T2.sgn_int) != 10)
	{
		T2.r1_int = T.r1_int;
		T2.r2_int = T.r2_int;
		T2.x1_int = T.x1_int;
		T2.x2_int = T.x2_int;
		T2.x3_int = T.x3_int;
		T2.y1_int = T.y1_int;
		T2.y2_int = T.y2_int;
		T2.y3_int = T.y3_int;

		if (abs(T2.sgn_int) < 10)
		{
			T2.x4_int = T.x4_int;
			T2.y4_int = T.y4_int;
			T2.x5_int = T.x5_int;
			T2.y5_int = T.y5_int;
			T2.u4_int = T.u4_int;
			T2.u4_int_ = T.u4_int_;
			T2.u5_int = T.u5_int;
			T2.u5_int_ = T.u5_int_;
		}

		if (abs(T2.sgn_int) > 10)
		{
			T2.v_int = T.v_int;
			T2.x5_int = T.x5_int;
			T2.y5_int = T.y5_int;
			T2.u5_int = T.u5_int;
			T2.u5_int_ = T.u5_int_;
		}
	}
}


void InhomogeneousThreeCavity::reloadtriex(Triangle &T2, Triangle &T)
{
	T2.sgn_int = T.sgn_intex;
	T2.x_int = T.x_intex;
	T2.y_int = T.y_intex;
	if (abs(T2.sgn_int) != 10)
	{
		T2.r1_int = T.r1_intex;
		T2.r2_int = T.r2_intex;
		T2.x1_int = T.x1_intex;
		T2.x2_int = T.x2_intex;
		T2.x3_int = T.x3_intex;
		T2.y1_int = T.y1_intex;
		T2.y2_int = T.y2_intex;
		T2.y3_int = T.y3_intex;

		if (abs(T2.sgn_int) < 10)
		{
			T2.x4_int = T.x4_intex;
			T2.y4_int = T.y4_intex;
			T2.x5_int = T.x5_intex;
			T2.y5_int = T.y5_intex;
			T2.u4_int = T.u4_intex;
			T2.u4_int_ = T.u4_intex_;
			T2.u5_int = T.u5_intex;
			T2.u5_int_ = T.u5_intex_;
		}

		if (abs(T2.sgn_int) > 10)
		{
			T2.v_int = T.v_intex;
			T2.x5_int = T.x5_intex;
			T2.y5_int = T.y5_intex;
			T2.u5_int = T.u5_intex;
			T2.u5_int_ = T.u5_intex_;
		}
	}
}


Vector4cd InhomogeneousThreeCavity::solveweak4(Triangle &T, double v2, double v3)
{
	//phi = para.phi_int;    !!!!

	double x1 = T.x(0);
	double y1 = T.y(0);
	double x2 = T.x(1);
	double y2 = T.y(1);
	double x3 = T.x(2);
	double y3 = T.y(2);

	if (abs(T.sgn_int) == 10)
	{
		Vector3cd inte;
		Vector4cd out;
		if (T.sgn_int == 10)
		{
			inte = weak4(T, v2, v3, true);
			out = Vector4cd(inte(0), inte(1), inte(2), 0);
		}
		else
		{
			inte = weak4(T, v2, v3, false);
			complex<double> con = 0;
			if (abs(phi_int(x1, y1)) < 1e-10) {
				con = con - inte(0)*a(x1, y1);
			}
			if (abs(phi_int(x2, y2)) < 1e-10) {
				con = con - inte(1)*a(x2, y2);
			}
			if (abs(phi_int(x3, y3)) < 1e-10) {
				con = con - inte(2)*a(x3, y3);
			}
			out = Vector4cd(inte(0), inte(1), inte(2), con);
		}
		return out;
	}//end if : abs(T.sgn_int) == 10

	if (abs(T.sgn_int) == 1)
	{
		Vector4d u4_ = T.u4_int_;
		double x4 = T.x4_int;
		double y4 = T.y4_int;

		Triangle T1;
		T1.x = Vector3d(x4, x2, x3);
		T1.y = Vector3d(y4, y2, y3);

		Vector3cd out4;
		if (T.sgn_int == 1) {
			out4 = weak4(T1, v2, v3, false);
		}
		else {
			out4 = weak4(T1, v2, v3, true);
		}
		complex<double> c1, c2, c3, c4;
		c1 = out4(0)*u4_(0);
		c2 = out4(0)*u4_(1) + out4(1);
		c3 = out4(0)*u4_(2) + out4(2);
		c4 = out4(0)*u4_(3);
		Vector4cd out(c1, c2, c3, c4);
		return out;
	}//end if : abs(T.sgn_int) == 1

	if (abs(T.sgn_int) == 2)
	{
		Vector4d u4 = T.u4_int;
		Vector4d u5 = T.u5_int;
		Vector4d u4_ = T.u4_int_;
		double x4 = T.x4_int;
		double y4 = T.y4_int;
		double x5 = T.x5_int;
		double y5 = T.y5_int;

		Triangle T1, T2;
		T1.x = Vector3d(x5, x2, x4);
		T1.y = Vector3d(y5, y2, y4);

		T2.x = Vector3d(x1, x4, x3);
		T2.y = Vector3d(y1, y4, y3);

		double r1 = T.r1_int;
		double v4 = v2 + (v3 - v2)*r1;

		Vector3cd out4, out5;
		if (T.sgn_int == 2) {
			out4 = weak4(T1, v2, v4, true);
			out5 = weak4(T2, v4, v3, false);
		}
		else {
			out4 = weak4(T1, v2, v4, false);
			out5 = weak4(T2, v4, v3, true);
		}

		complex<double> c1, c2, c3, c4;
		c1 = out4(0)*u5(2) + out4(2)*u4(2) + out5(0) + out5(1)*u4_(2);
		c2 = out4(0)*u5(0) + out4(1) + out4(2)*u4(0) + out5(1)*u4_(0);
		c3 = out4(0)*u5(1) + out4(2)*u4(1) + out5(1)*u4_(1) + out5(2);
		c4 = out4(0)*u5(3) + out4(2)*u4(3) + out5(1)*u4_(3);
		Vector4cd out(c1, c2, c3, c4);
		return out;
	}//end if : abs(T.sgn_int) == 2

	if (abs(T.sgn_int) == 3)
	{
		Vector4d u4 = T.u4_int;
		Vector4d u5 = T.u5_int;
		Vector4d u4_ = T.u4_int_;
		Vector4d u5_ = T.u5_int_;

		double x4 = T.x4_int;
		double y4 = T.y4_int;
		double x5 = T.x5_int;
		double y5 = T.y5_int;

		Triangle T1, T2;
		T1.x = Vector3d(x4, x2, x5);
		T1.y = Vector3d(y4, y2, y5);

		T2.x = Vector3d(x4, x5, x3);
		T2.y = Vector3d(y4, y5, y3);

		double r2 = T.r2_int;
		double v5 = v3 + (v2 - v3)*r2;

		Vector3cd out4, out5;
		if (T.sgn_int == 3) {
			out4 = weak4(T1, v2, v5, false);
			out5 = weak4(T2, v5, v3, true);
		}
		else {
			out4 = weak4(T1, v2, v5, true);
			out5 = weak4(T2, v5, v3, false);
		}

		complex<double> c1, c2, c3, c4;
		c1 = out4(1)*u4_(2) + out4(3)*u5_(2) + out5(1)*u4(2) + out5(2)*u5(2);
		c2 = out4(1)*u4_(3) + out4(2) + out4(3)*u5_(3) + out5(1)*u4(3) + out5(2)*u5(3);
		c3 = out4(1)*u4_(1) + out4(3)*u5_(1) + out5(1)*u4(1) + out5(2)*u5(1) + out5(3);
		c4 = out4(1)*u4_(4) + out4(3)*u5_(4) + out5(1)*u4(4) + out5(2)*u5(4);

		Vector4cd out(c1, c2, c3, c4);
		return out;
	}//end if : abs(T.sgn_int) == 3

	if (abs(T.sgn_int) == 11)
	{
		double x5 = T.x5_int;
		double y5 = T.y5_int;

		Triangle T1;
		T1.x = Vector3d(x5, x2, x3);
		T1.y = Vector3d(y5, y2, y3);

		Vector3cd out4;
		Vector4d u5;
		if (T.sgn_int == 11) {
			out4 = weak4(T1, v2, v3, false);
			u5 = T.u5_int_;
		}
		else {
			out4 = weak4(T1, v2, v3, true);
			u5 = T.u5_int;
		}

		complex<double> c1, c2, c3, c4;
		c1 = out4(0)*u5(0);
		c2 = out4(0)*u5(1) + out4(1);
		c3 = out4(0)*u5(2) + out4(2);
		c4 = out4(0)*u5(3);

		if (T.sgn_int == 11) {
			c4 = c4 - out4(1)*a(x2, y2);
		}

		Vector4cd out(c1, c2, c3, c4);
		return out;
	}//end if : abs(T.sgn_int) == 11

	if (abs(T.sgn_int) == 12)
	{
		double x5 = T.x5_int;
		double y5 = T.y5_int;

		Triangle T1;
		T1.x = Vector3d(x5, x2, x3);
		T1.y = Vector3d(y5, y2, y3);

		Vector3cd out4;
		Vector4d u5;
		if (T.sgn_int == 12) {
			out4 = weak4(T1, v2, v3, true);
			u5 = T.u5_int;
		}
		else {
			out4 = weak4(T1, v2, v3, false);
			u5 = T.u5_int_;
		}

		complex<double> c1, c2, c3, c4;
		c1 = out4(0)*u5(2);
		c2 = out4(0)*u5(0) + out4(1);
		c3 = out4(0)*u5(1) + out4(2);
		c4 = out4(0)*u5(3);

		if (T.sgn_int == -12) {
			c4 = c4 - out4(2)*a(x3, y3);
		}
		Vector4cd out(c1, c2, c3, c4);
		return out;
	}//end if : abs(T.sgn_int) == 12

	if (abs(T.sgn_int) == 13)
	{
		Vector4d u5 = T.u5_int;
		Vector4d u5_ = T.u5_int_;
		double r2 = T.r2_int;
		double v5 = v3 + (v2 - v3)*r2;
		double x5 = T.x5_int;
		double y5 = T.y5_int;

		Triangle T1, T2;
		T1.x = Vector3d(x1, x2, x5);
		T1.y = Vector3d(y1, y2, y5);
		T2.x = Vector3d(x1, x5, x3);
		T2.y = Vector3d(y1, y5, y3);

		complex<double> c1, c2, c3, c4;
		Vector3cd out4, out5;
		if (T.sgn_int == 13)
		{
			out4 = weak4(T1, v2, v5, false);
			out5 = weak4(T2, v5, v3, true);
			c1 = out4(0) + out4(2)*u5_(1) + out5(0) + out5(1)*u5(1);
			c2 = out4(1) + out4(2)*u5_(2) + out5(1)*u5(2);
			c3 = out4(2)*u5_(0) + out5(1)*u5(0) + out5(2);
			c4 = out4(2)*u5_(3) + out5(1)*u5(3);
			c4 = c4 - out4(0)*a(x1, y1);
		}
		else
		{
			out4 = weak4(T1, v2, v5, true);
			out5 = weak4(T2, v5, v3, false);
			c1 = out4(0) + out4(2)*u5(1) + out5(0) + out5(1)*u5_(1);
			c2 = out4(1) + out4(2)*u5(2) + out5(1)*u5_(2);
			c3 = out4(2)*u5(0) + out5(1)*u5_(0) + out5(2);
			c4 = out4(2)*u5(3) + out5(1)*u5_(3);
			c4 = c4 - out5(0)*a(x1, y1);
		}
		Vector4cd out(c1, c2, c3, c4);
		return out;
	}//end if : abs(T.sgn_int) == 13
}


int InhomogeneousThreeCavity::findXinNbound_multiple(vector<vector<double>> nbound, double x, int start, int end)
{
	int size = nbound[0].size();
	int index = -1;
	bool hasFind = false;
	for (int i = start; i <= end; i++)
	{
		if (abs(nbound[1][i] - x) < 1e-10)
		{
			index = i-start;
			hasFind = true;
			break;
		}
	}
	if (hasFind)
		return index;
	else
		throw"not find";
}