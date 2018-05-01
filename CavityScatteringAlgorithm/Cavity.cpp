#include "Cavity.h"


Cavity::Cavity(unsigned int cavityType)
{
	this->cavityType = cavityType;

	//启动matlab引擎
	cout << endl << "正在启动matlab引擎,请稍后..." << endl;
	if (!(this->ep = engOpen(NULL)))
	{
		fprintf(stderr, "\n无法启动MATLAB引擎\n");
		return;
	}
}

Cavity::~Cavity()
{
	engClose(this->ep);
}


void Cavity::InitVirtualBorder(double top, double bottom, double left, double right)
{
	this->virtualBorderTop = top;
	this->virtualBorderBottom = bottom;
	this->virtualBorderLeft = left;
	this->virtualBorderRight = right;

	//Mark function InitVirtualBorder has been executed
	this->initalCheckKey = this->initalCheckKey | 1; 
}

void Cavity::InitMesh(int m, int n)
{
	// check if the first bit of the binary array initalCheckKey is 1
	if (this->initalCheckKey & 1)
	{
		// This shows that function InitVirtualBorder has been executed
		this->meshWidth = m;
		this->meshHeight = n;
		this->stepX = (virtualBorderRight - virtualBorderLeft) / m;
		this->stepY = (virtualBorderTop - virtualBorderBottom) / n;

		//Mark function InitMesh has been executed
		this->initalCheckKey = this->initalCheckKey | 2;
	}
	else
	{
		throw "Init Error!\n You should call InitVirtualBorder before use InitMesh.";
	}
}

/*
InitialCheck: Check if the cavity parameters are complete initialization
checkLog: record the the check result in string.
Return whether to pass the check(1:pass, 0:faile).
--------------------------------------------------------------------------
Detailed explanation of the checkLog
the 1th bit :Mark function InitVirtualBorder has been executed
the 2th bit :Mark function InitMesh has been executed
the 3th bit :Mark function InitElectromagneticParameter has been executed
the 4th bit :Mark function InitCavityShapeParameter has been executed
the 5th bit :Mark function  has been executed
*/
bool Cavity::InitialCheck(char *checkLog) {
	//this->initalCheckKey
	return true;
}


void Cavity::PlotAperture(string title, string xlabel, string ylabel, int sign)
{
	cout << endl << "正在调用matlab引擎进行绘图,请稍后..." << endl;

	///获取用于绘图的横纵坐标数据
	VectorXd plotX, plotY;
	drawAperture(plotX, plotY, nbound, sign);
	//cout << plotX << endl;
	//cout << plotY << endl;


	///获取程序当前工作路径
	char cur_path[FILENAME_MAX];
	getcwd(cur_path, FILENAME_MAX);
	cout << "当前路径：" << cur_path << endl;


	///构造修改当前工作路径的matlab命令
	// 例如：cd(\"F:\\计算数学\\CavityScatteringAlgorithm\\Test\");
	int len_path = -1; //获取路径的实际长度
	for (int i = 0; i < size(cur_path); i++)
	{
		if (cur_path[i] == '\0')
		{
			len_path = i;
			break;
		}
	}

	const char *cmd_start = "cd(\"";
	const char *cmd_end = "\");";
	int len_start = strlen(cmd_start);
	int len_end = strlen(cmd_end);

	//组合出整个命令
	//首先分配内存
	char *cdCommand = (char*)malloc(len_start + len_path + len_end + 1);
	if (!cdCommand)
		abort();
	//拼接
	memcpy(cdCommand, cmd_start, len_start);
	memcpy(cdCommand + len_start, cur_path, len_path);
	memcpy(cdCommand + len_start + len_path, cmd_end, len_end);
	cdCommand[len_start + len_path + len_end] = '\0';


	///调用matalb引擎画图
	mxArray *mx_xValue, *mx_apertureValue;

	//将C中数组转化为matlab数组mxArray
	mx_xValue = mxCreateDoubleMatrix(1, plotX.size(), mxREAL);
	//memcpy(mxGetPr(matrix), m, sizeof(double) * 3 * 3);
	double *xValuePr = mxGetPr(mx_xValue);
	for (int i = 0; i < plotX.size(); i++)
		*xValuePr++ = plotX(i);

	mx_apertureValue = mxCreateDoubleMatrix(1, plotX.size(), mxREAL);
	//memcpy(mxGetPr(matrix), m, sizeof(double) * 3 * 3);
	double *apertureValuePr = mxGetPr(mx_apertureValue);
	for (int i = 0; i < plotY.size(); i++)
		*apertureValuePr++ = plotY(i);


	///将需要的变量加入matlab引擎
	engPutVariable(ep, "xValue", mx_xValue);
	engPutVariable(ep, "apertureValue", mx_apertureValue);


	///调用matlab的cd命令
	//相当于调用类似如下命令：engEvalString(ep, "cd(\"F:\\计算数学\\CavityScatteringAlgorithm\\Test\");");
	engEvalString(ep, cdCommand);
	//cout << cdCommand << endl;

	//以弃用，修改起来太麻烦，不够灵活
	///调用编写的matlab函数plotAperture（plotAperture.m文件应存放于当前路径下）
	//engEvalString(ep, "plotAperture( xValue, apertureValue)");

	//直接调用matlab的内置函数绘图

	engEvalString(ep, "figure");
	engEvalString(ep, "plot(xValue, apertureValue, 'b');");

	engEvalString(ep, "xlabel('x', 'FontSize', 12);");

	if (sign == 0)//幅值
		engEvalString(ep, "ylabel('Magnitude', 'FontSize', 12);");
	else if (sign == 1)//实部
		engEvalString(ep, "ylabel('Real', 'FontSize', 12);");
	else if (sign == -1)//虚部
		engEvalString(ep, "ylabel('Imaginary', 'FontSize', 12);");
	else if (sign == 10)//相位
		engEvalString(ep, "ylabel('Phase', 'FontSize', 12);");
	else
	{
		;
	}

	engEvalString(ep, "set(gca, 'LineWidth', 1.5);");

	//set(gca, 'XLim', [para.box.left para.box.right]);
	//set(gca, 'XTick', -3:1 : 3);
	//set(gca, 'YLim', [0 1.5]);
	//set(gca, 'YTick', 0:0.5 : 1.5);
	//set(gca, 'YLim', [-20 50]);
	//set(gca, 'YTick', -20:10 : 50);


	//如无法成功绘图，需要重新注册matlab引擎
	//详见https://blog.csdn.net/xiaoqiang920/article/details/8949254
}





//---------------------------------protect---------------------------------

//------------------------------
//主要功能函数
//------------------------------
VectorXcd Cavity::solveX(SparseMatrix<complex<double>> &A, VectorXcd &B)
{
	//求解Ax = b
	VectorXcd x(B.size());
	//SolverClassName<SparseMatrix<double> > solver;

	BiCGSTAB<SparseMatrix <complex<double>>> solver;
	//IncompleteLUT<SparseMatrix <complex<double>> > solver;
	solver.compute(A);
	//solver.preconditioner=IncompleteLUT;
	if (solver.info() != Success) {
		throw "decomposition failed";
	}
	x = solver.solve(B);
	if (solver.info() != Success) {
		throw "solving failed";
	}
	return x;
}

VectorXcd Cavity::solveX_mx(mxArray *A_mx, mxArray *B_mx)
{
	///将需要的变量加入matlab引擎
	engPutVariable(this->ep, "A", A_mx);
	engPutVariable(this->ep, "B", B_mx);

	///调用matlab命令
	//将A转化为matlab中的稀疏矩阵(调用前不要忘记A需要转制，原因详见setA_mx函数最后)
	engEvalString(ep, "A = spconvert(A.');");  //此外，注意转置应使用A.'   （对于复矩阵A'是共轭转置） 
	// 解方程
	engEvalString(ep, "u_mx = A\\B;"); 

	//获取matlab求解结果
	mxArray *u_mx = engGetVariable(this->ep, "u_mx");

	//将matlab数组mxArray转化为C数组
	int u_size = mxGetNumberOfElements(u_mx);  //获得解向量的长度
	VectorXcd x(u_size);
	double *UvaluePr = mxGetPr(u_mx);
	double *UvaluePi = mxGetPi(u_mx);
	for (int i = 0; i < u_size; i++)
	{
		x(i) = complex<double>(*UvaluePr++, *UvaluePi++);
	}
	return x;
}



//------------------------------
//Auxiliary function 其它辅助函数
//------------------------------

void Cavity::checkVirtualBorder()
{

}

double Cavity::findzero(double x1, double y1, double x2, double y2)
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
	double left_phi = phi(x1, y1);    //左端点水平集
	double right_phi = phi(x2, y2);   //右端点水平集

	// 初始化
	double cut_old_x = (left_x + right_x)*0.5;       // 初始化起始分割点cut为(x1,y1)、(x2,y2)中点
	double cut_old_y = (left_y + right_y)*0.5;
	double cut_phi = phi(cut_old_x, cut_old_y);      // 起始分割点cut的水平集函数值

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
			right_phi = phi(right_x, right_y);          // 修改右端点水平集
		}
		else if (right_phi * cut_phi < 0)
		{
			// 交点位于cut_old——right_x
			cut_new_x = (cut_old_x + right_x)*0.5;     // 新的中点
			cut_new_y = (cut_old_y + right_y)*0.5;
			left_x = cut_old_x;                        // 修改右端点
			left_y = cut_old_y;
			left_phi = phi(left_x, left_y);            // 修改右端点水平集
		}
		else
		{
			throw "something wrong in findzero.m";
		}
		cut_old_x = cut_new_x;
		cut_old_y = cut_new_y;
		cut_phi = phi(cut_old_x, cut_old_y);
	}
	double d_old_1 = sqrt(pow(cut_old_x - x1, 2) + pow(cut_old_y - y1, 2));
	double d_1_2 = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));

	return d_old_1 / d_1_2;
}


void Cavity::initializeComplexVector(vector<complex<double>> &vector, complex<double> InitialValue, int vectorSize)
{
	if (vector.size() == 0)
		vector.resize(vectorSize);
	if (vector.size() != vectorSize)
		throw "Size conflict when initialize complex vector.";

	for (int i = 0; i < vector.size(); i++)
		vector[i] = InitialValue;
}

void Cavity::complexVectorAdd(vector<complex<double>> &vector1, vector<complex<double>> vector2)
{
	if(vector1.size()==0 || vector2.size()==0)
		throw "complex vector is empty";
	if(vector1.size()!=vector2.size())
		throw "Size conflict when Add complex vector.";

	for (int i = 0; i < vector1.size(); i++)
		vector1[i] += vector2[i];
}

double Cavity::calculateTriangleArea(Vector3d tri1, Vector3d tri2)
{
	//tri1三角形三个顶点的x坐标
	//tri2三角形三个顶点的y坐标（与tri1相互对应）

	double a, b, c;//三角形三边长
	a = sqrt(pow(tri1(1) - tri1(0), 2) + pow(tri2(1) - tri2(0), 2));
	b = sqrt(pow(tri1(2) - tri1(0), 2) + pow(tri2(2) - tri2(0), 2));
	c = sqrt(pow(tri1(2) - tri1(1), 2) + pow(tri2(2) - tri2(1), 2));

	double p;//半周长
	double s;//面积
	if (a + b > c && a + c > b && b + c > a) //判断是否可以构成三角形。
	{
		p = (a + b + c) / 2;//计算半周长
		s = sqrt(p*(p - a)*(p - b)*(p - c));//套用海伦公式，计算面积
		return s;
	}
	else
		throw "给定坐标不能构成三角形";
}

vector<complex<double>> Cavity::eigenVector2stdVector(VectorXcd eigenVector)
{
	int size = eigenVector.size();

	vector<complex<double>> stdVector(size);
	for (int i = 0; i < size; i++)
	{
		stdVector[i] = eigenVector(i);
	}
	return stdVector;
}

int Cavity::findXinNbound(vector<vector<double>> nbound, double x)
{
	int size = nbound[0].size();
	int index = -1;
	bool hasFind = false;
	for (int i = 0; i < size; i++)
	{
		if (abs(nbound[1][i] - x) < 1e-10)
		{
			index = i;
			hasFind = true;
			break;
		}
	}
	if (hasFind)
		return index;
	else
		throw"not find";
}

int Cavity::findXIndexinNbound(vector<vector<double>> nbound, double xIndex)
{
	int size = nbound[0].size();
	int index = -1;
	bool hasFind = false;
	for (int i = 0; i < size; i++)
	{
		if (abs(nbound[0][i] - xIndex) < 1e-10)
		{
			index = i;
			hasFind = true;
			break;
		}
	}
	if (hasFind)
		return index;
	else
		throw"not find";
}

SparseMatrix<complex<double>> Cavity::buildSparseMatrixA(double* pA, int MatrixSize, int ArowNum)
{
	//构造Ax=b中的稀疏矩阵A

	//-----------------------------------
	//方法1：使用insert函数逐元素插入
	//经实验，效率很低，因此弃用
	//-----------------------------------
	//SparseMatrix<complex<double>> A(MatrixSize, MatrixSize);
	//int ArowsNum = ArowNum;
	//int AcolsNum = 4;
	//int xIndex, yIndex;
	//double real, ima; //实部，虚部
	//for (int i = 0; i < ArowsNum; i++)
	//{
	//	xIndex = *(pA++);   //这个索引是从1计数的(因此使用时不要忘记-1)
	//	yIndex = *(pA++);
	//	real= *(pA++);
	//	ima = *(pA++);
	//	A.insert(xIndex-1, yIndex-1) = complex<double>(real, ima);
	//}
	//return A;



	//-----------------------------------
	//方法2：使用三元组构造稀疏矩阵
	//-----------------------------------
	int estimation_of_entries = ArowNum; //预计非零元素的个数
	vector<Triplet<complex<double>>> tripletList;   //复数三元组
	tripletList.reserve(estimation_of_entries);

	int xIndex, yIndex;
	double real, ima; //实部，虚部
	for (int i = 0; i < ArowNum; i++)
	{
		xIndex = *(pA++);   //这个索引是从1计数的(因此使用时不要忘记-1)
		yIndex = *(pA++);
		real= *(pA++);
		ima = *(pA++);
		complex<double> v_ij(real,ima);
		//cout << xIndex << "  " << yIndex << "  " << v_ij << endl;
		tripletList.push_back(Triplet<complex<double>>(xIndex-1, yIndex-1, v_ij));
	}

	SparseMatrix<complex<double>> A(MatrixSize, MatrixSize);
	A.setFromTriplets(tripletList.begin(), tripletList.end()); //根据三元组列表生成稀疏矩阵	
	return A;
}

int Cavity::sign(double num)
{
	int sign;
	if (num > 0)
		sign = 1;
	else if (num == 0)
		sign = 0;
	else
		sign = -1;
	return sign;
}

