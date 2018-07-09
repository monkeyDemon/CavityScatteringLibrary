#include "Cavity.h"


Cavity::Cavity(unsigned int cavityType)
{
	this->cavityType = cavityType;

	//����matlab����
	cout << endl << "��������matlab����,���Ժ�..." << endl;
	if (!(this->ep = engOpen(NULL)))
	{
		fprintf(stderr, "\n�޷�����MATLAB����\n");
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
	cout << endl << "���ڵ���matlab������л�ͼ,���Ժ�..." << endl;

	///��ȡ���ڻ�ͼ�ĺ�����������
	VectorXd plotX, plotY;
	drawAperture(plotX, plotY, nbound, sign);
	//cout << plotX << endl;
	//cout << plotY << endl;


	///��ȡ����ǰ����·��
	char cur_path[FILENAME_MAX];
	getcwd(cur_path, FILENAME_MAX);
	cout << "��ǰ·����" << cur_path << endl;


	///�����޸ĵ�ǰ����·����matlab����
	// ���磺cd(\"F:\\������ѧ\\CavityScatteringAlgorithm\\Test\");
	int len_path = -1; //��ȡ·����ʵ�ʳ���
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

	//��ϳ���������
	//���ȷ����ڴ�
	char *cdCommand = (char*)malloc(len_start + len_path + len_end + 1);
	if (!cdCommand)
		abort();
	//ƴ��
	memcpy(cdCommand, cmd_start, len_start);
	memcpy(cdCommand + len_start, cur_path, len_path);
	memcpy(cdCommand + len_start + len_path, cmd_end, len_end);
	cdCommand[len_start + len_path + len_end] = '\0';


	//����matalb���滭ͼ
	mxArray *mx_xValue, *mx_apertureValue;

	//��C������ת��Ϊmatlab����mxArray
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


	///����Ҫ�ı�������matlab����
	engPutVariable(ep, "xValue", mx_xValue);
	engPutVariable(ep, "apertureValue", mx_apertureValue);


	///����matlab��cd����
	//�൱�ڵ��������������engEvalString(ep, "cd(\"F:\\������ѧ\\CavityScatteringAlgorithm\\Test\");");
	engEvalString(ep, cdCommand);
	//cout << cdCommand << endl;

	//�����ã��޸�����̫�鷳���������
	///���ñ�д��matlab����plotAperture��plotAperture.m�ļ�Ӧ����ڵ�ǰ·���£�
	//engEvalString(ep, "plotAperture( xValue, apertureValue)");

	//ֱ�ӵ���matlab�����ú�����ͼ

	engEvalString(ep, "figure");
	engEvalString(ep, "plot(xValue, apertureValue, 'b');");

	engEvalString(ep, "xlabel('x', 'FontSize', 12);");

	if (sign == 0)//��ֵ
		engEvalString(ep, "ylabel('Magnitude', 'FontSize', 12);");
	else if (sign == 1)//ʵ��
		engEvalString(ep, "ylabel('Real', 'FontSize', 12);");
	else if (sign == -1)//�鲿
		engEvalString(ep, "ylabel('Imaginary', 'FontSize', 12);");
	else if (sign == 10)//��λ
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


	//���޷��ɹ���ͼ����Ҫ����ע��matlab����
	//���https://blog.csdn.net/xiaoqiang920/article/details/8949254
}


//��Ⲣ����RCS
//interval������RCS�ļ����Ŀǰ�ĺϷ�ȡֵ�У�0.125��0.25��0.5��1
void Cavity::SolveRCS(double interval)
{
	if (interval != 0.125 && interval != 0.25 && interval != 0.5 && interval != 1)
		throw("interval���ò��Ϸ�");

	//��ȡRCS������
	double temp = (int)(1 / interval);
	int rcs_length = 89 * temp + 1;

	//��ʼ��RCS����ⷶΧ
	VectorXd rcs_range(rcs_length);//rcs_range = 0:-0.25 : -89;
	temp = 0;
	for (int i = 0; i < rcs_length; i++)
	{
		rcs_range(i) = temp;
		temp -= interval;
	}

	//��ʼ��RCS�����
	VectorXd RCS(rcs_length);

	//��ʼ���㣡
	MyTimer myTimer(1);//��ʼ����ʱ��
	myTimer.Start("RCS");
	for (int index_rcs = 0; index_rcs < rcs_length; index_rcs++)
	{
		double aj = rcs_range(index_rcs);
		this->theta = aj*M_PI / 180;

		//��⵱ǰthetaʱ�ھ���Ľ�
		cout << "============================" << endl;
		cout << "theta = " << theta << endl;
		cout << "============================" << endl;
		this->Solve();
		cout << endl;

		// ���㵱ǰtheta��RCS
		int aperNum = this->nbound[0].size();
		complex<double> Ptheta = 0;
		double as, ac;
		complex<double> temp_i(0, 1);
		for (int k = 0; k < aperNum; k++)
		{
			double x_aperture = nbound[1][k];
			double ang = -theta + M_PI / 2;
			as = sin(ang);
			ac = cos(ang);
			Ptheta += this->solutionOfAperture(k)*(cos(ac*k0*x_aperture) + sin(ac*k0*x_aperture)*temp_i);
		}
		Ptheta = 0.5*k0*as*(1.0 / (aperNum + 1))*Ptheta;              // P(theta), defined in Zhao's PhD thesis, P26 (2.120)
		double gama_theta = 4.0*abs(Ptheta)*abs(Ptheta) / k0;   // gama(theta), defined in Zhao's PhD thesis, P26 (2.119)
		double ttt = log(gama_theta);
		RCS(index_rcs) = 10.0*log(gama_theta) / log(10.0);        // RCS
	}
	myTimer.EndAndPrint();


	// ����RCS
	cout << endl << "���ڵ���matlab������л�ͼ,���Ժ�..." << endl;
	//���޷��ɹ���ͼ����Ҫ����ע��matlab����
	//���https://blog.csdn.net/xiaoqiang920/article/details/8949254

	//��ʼ��RCS���Ƶ�X������
	VectorXd plotX(rcs_length);
	temp = 91;
	for (int i = 0; i < rcs_length; i++)
	{
		plotX(i) = temp;
		temp += interval;
	}
	//cout << plotX << endl;
	//cout << RCS << endl;


	///����matalb���滭ͼ
	mxArray *mx_xValue, *mx_RCS;

	//��C������ת��Ϊmatlab����mxArray
	mx_xValue = mxCreateDoubleMatrix(1, plotX.size(), mxREAL);
	//memcpy(mxGetPr(matrix), m, sizeof(double) * 3 * 3);
	double *xValuePr = mxGetPr(mx_xValue);
	for (int i = 0; i < plotX.size(); i++)
		*xValuePr++ = plotX(i);

	mx_RCS = mxCreateDoubleMatrix(1, RCS.size(), mxREAL);
	double *RCSValuePr = mxGetPr(mx_RCS);
	for (int i = 0; i < RCS.size(); i++)
		*RCSValuePr++ = RCS(i);


	///����Ҫ�ı�������matlab����
	engPutVariable(ep, "xValue", mx_xValue);
	engPutVariable(ep, "RCS", mx_RCS);

	//ֱ�ӵ���matlab�����ú�����ͼ
	engEvalString(ep, "figure");
	//engEvalString(ep, "plot(91 : 0.25 : 180, RCS, 'k-', 'LineWidth', 1.2);");
	engEvalString(ep, "plot(xValue, RCS, 'k-', 'LineWidth', 1.2);");
	//engEvalString(ep, "ylim([-70, 20]);");
	//engEvalString(ep, "set(gca, 'YTick', -70:10 : 20, 'LineWidth', 1.2);");
	engEvalString(ep, "xlabel('Observation angle(degree)', 'FontSize', 12);");
	engEvalString(ep, "ylabel('Backscatter RCS(dB)', 'FontSize', 12);");
}




//---------------------------------protect---------------------------------

//------------------------------
//��Ҫ���ܺ���
//------------------------------
VectorXcd Cavity::solveX(SparseMatrix<complex<double>> &A, VectorXcd &B)
{
	//���Ax = b
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
	///����Ҫ�ı�������matlab����
	engPutVariable(this->ep, "A", A_mx);
	engPutVariable(this->ep, "B", B_mx);

	///����matlab����
	//��Aת��Ϊmatlab�е�ϡ�����(����ǰ��Ҫ����A��Ҫת�ƣ�ԭ�����setA_mx�������)
	engEvalString(ep, "A = spconvert(A.');");  //���⣬ע��ת��Ӧʹ��A.'   �����ڸ�����A'�ǹ���ת�ã� 
	// �ⷽ��
	engEvalString(ep, "u_mx = A\\B;"); 

	//��ȡmatlab�����
	mxArray *u_mx = engGetVariable(this->ep, "u_mx");

	//��matlab����mxArrayת��ΪC����
	int u_size = mxGetNumberOfElements(u_mx);  //��ý������ĳ���
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
//Auxiliary function ������������
//------------------------------

void Cavity::checkVirtualBorder()
{

}

double Cavity::findzero(double x1, double y1, double x2, double y2)
{
	// (x1,y1)��(x2,y2)λ��interface�����࣬�ʼ���interface�����߶�(x1,y1)��(x2,y2)�ڵ�cut��
	// ���ڦ��ϵĵ�ˮƽ������ֵΪ0
	// findzero����ʹ�ö��ַ�Ѱ�ҵ�cut��ֱ�����㾫��Ҫ�� | phi(cut) | <eps
	// out���ؽ���Ķ�Ӧλ�ã���(x1,y1)->(x2,y2)�ֱ��Ӧ0->1

	double eps = 1e-10;

	double left_x = x1;
	double left_y = y1;
	double right_x = x2;
	double right_y = y2;
	double left_phi = phi(x1, y1);    //��˵�ˮƽ��
	double right_phi = phi(x2, y2);   //�Ҷ˵�ˮƽ��

	// ��ʼ��
	double cut_old_x = (left_x + right_x)*0.5;       // ��ʼ����ʼ�ָ��cutΪ(x1,y1)��(x2,y2)�е�
	double cut_old_y = (left_y + right_y)*0.5;
	double cut_phi = phi(cut_old_x, cut_old_y);      // ��ʼ�ָ��cut��ˮƽ������ֵ

	// ����ֱ�����㾫��Ҫ��
	double cut_new_x, cut_new_y;
	while (abs(cut_phi) >= eps)
	{
		//�����㾫��Ҫ�󣬼�������Ѱ��
		
		// �����ж������Ľ���λ����һ�ߣ�left����cut_old or cut_old����right
		if (left_phi * cut_phi < 0)
		{
			// ����λ�� left����cut_old
			cut_new_x = (cut_old_x + left_x)*0.5;       // �µ��е�
			cut_new_y = (cut_old_y + left_y)*0.5;
			right_x = cut_old_x;                        // �޸��Ҷ˵�
			right_y = cut_old_y;
			right_phi = phi(right_x, right_y);          // �޸��Ҷ˵�ˮƽ��
		}
		else if (right_phi * cut_phi < 0)
		{
			// ����λ��cut_old����right_x
			cut_new_x = (cut_old_x + right_x)*0.5;     // �µ��е�
			cut_new_y = (cut_old_y + right_y)*0.5;
			left_x = cut_old_x;                        // �޸��Ҷ˵�
			left_y = cut_old_y;
			left_phi = phi(left_x, left_y);            // �޸��Ҷ˵�ˮƽ��
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
	//tri1���������������x����
	//tri2���������������y���꣨��tri1�໥��Ӧ��

	double a, b, c;//���������߳�
	a = sqrt(pow(tri1(1) - tri1(0), 2) + pow(tri2(1) - tri2(0), 2));
	b = sqrt(pow(tri1(2) - tri1(0), 2) + pow(tri2(2) - tri2(0), 2));
	c = sqrt(pow(tri1(2) - tri1(1), 2) + pow(tri2(2) - tri2(1), 2));

	double p;//���ܳ�
	double s;//���
	if (a + b > c && a + c > b && b + c > a) //�ж��Ƿ���Թ��������Ρ�
	{
		p = (a + b + c) / 2;//������ܳ�
		s = sqrt(p*(p - a)*(p - b)*(p - c));//���ú��׹�ʽ���������
		return s;
	}
	else
		throw "�������겻�ܹ���������";
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
	//����Ax=b�е�ϡ�����A

	//-----------------------------------
	//����1��ʹ��insert������Ԫ�ز���
	//��ʵ�飬Ч�ʺܵͣ��������
	//-----------------------------------
	//SparseMatrix<complex<double>> A(MatrixSize, MatrixSize);
	//int ArowsNum = ArowNum;
	//int AcolsNum = 4;
	//int xIndex, yIndex;
	//double real, ima; //ʵ�����鲿
	//for (int i = 0; i < ArowsNum; i++)
	//{
	//	xIndex = *(pA++);   //��������Ǵ�1������(���ʹ��ʱ��Ҫ����-1)
	//	yIndex = *(pA++);
	//	real= *(pA++);
	//	ima = *(pA++);
	//	A.insert(xIndex-1, yIndex-1) = complex<double>(real, ima);
	//}
	//return A;



	//-----------------------------------
	//����2��ʹ����Ԫ�鹹��ϡ�����
	//-----------------------------------
	int estimation_of_entries = ArowNum; //Ԥ�Ʒ���Ԫ�صĸ���
	vector<Triplet<complex<double>>> tripletList;   //������Ԫ��
	tripletList.reserve(estimation_of_entries);

	int xIndex, yIndex;
	double real, ima; //ʵ�����鲿
	for (int i = 0; i < ArowNum; i++)
	{
		xIndex = *(pA++);   //��������Ǵ�1������(���ʹ��ʱ��Ҫ����-1)
		yIndex = *(pA++);
		real= *(pA++);
		ima = *(pA++);
		complex<double> v_ij(real,ima);
		//cout << xIndex << "  " << yIndex << "  " << v_ij << endl;
		tripletList.push_back(Triplet<complex<double>>(xIndex-1, yIndex-1, v_ij));
	}

	SparseMatrix<complex<double>> A(MatrixSize, MatrixSize);
	A.setFromTriplets(tripletList.begin(), tripletList.end()); //������Ԫ���б�����ϡ�����	
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

