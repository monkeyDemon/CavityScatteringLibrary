#pragma once

//#include<stdlib.h>
#include<iostream>
#include <complex>
#include <math.h>
#include <direct.h>

#include <Eigen/Dense>//Eigen���Դ�����
#include <Eigen/Sparse>

#include "engine.h" //���MATLAB����ͷ�ļ�

#include "Structure.h"
#include "ApertureIntegral.h"
#include "MyTimer.h"

using namespace Eigen;
using namespace std;

class _declspec(dllexport) Cavity
{
public:
	Cavity(unsigned int cavityType); //���캯��
	~Cavity();                                     //��������
	void InitVirtualBorder(double top, double bottom, double left, double right);
	void InitMesh(int m, int n);
	//void InitElectromagneticParameter() = 0;
	//void InitAperture(double apertureLeft, double apertureRight, double apertureY);
	//void InitCavityShapeParameter() = 0;
	bool InitialCheck(char *checkLog);
	virtual bool Solve()=0;
	virtual bool SolveAperture(string title, string xlabel, string ylabel, int sign) = 0;
	void SolveRCS(double interval);

	/*
	��������������
	*/
	virtual void PlotTriangleMesh(string title, string xlabel, string ylabel) = 0;

protected:
	Engine *ep;//matlab����

	unsigned int initalCheckKey = 0;   //ʹ�ô˶����ƴ����г�ʼ�����
	unsigned int cavityType;

	int nn;
	vector<int> nu;
	vector<vector<double>> nbound;

	// Solution of the Holmholtz equation on aperture
	VectorXcd solutionOfAperture;

	//�����
	double theta;
	//����
	double k0;

	// parameter of the virtual border which cover the cavity(initial by InitVirtualBorder)
	double virtualBorderTop;
	double virtualBorderBottom;
	double virtualBorderLeft;
	double virtualBorderRight;

	// parameter of the triangle grid(initial by InitMesh)
	int meshWidth;
	int meshHeight;
	double stepX;
	double stepY;


	//����������ͷ����йصĺ���
	virtual Matrix2d beta(double x, double y)=0;
	virtual complex<double> q(double x, double y)=0;
	virtual complex<double> f(double x, double y) = 0;
	virtual double phi(double x, double y) = 0;
	virtual complex<double> u(double x, double y) = 0;

	//��Ҫ���ܺ���
	virtual void setTri(TriangleMesh &U, TriangleMesh &L, int &nn, vector<int> &nu, vector<vector<double>> &nbound)=0;
	virtual VectorXcd compute_g(MatrixXcd &G, vector<vector<double>> &nbound)=0;
	virtual vector<vector<gridCell>> setGrid(TriangleMesh &U, TriangleMesh &L, vector<vector<double>> &nbound, vector<int> &nu)=0;
	virtual void weak(Triangle_Normal &T, Vector3d v, int topsign, vector<vector<double>> &nbound, vector <complex<double>> &out, vector <complex<double>> &out5m) = 0;
	virtual void weak(Triangle_All &T, Vector3d v, int topsign, vector<vector<double>> &nbound, vector <complex<double>> &out, vector <complex<double>> &out5m) = 0;
	virtual Vector3cd weak1(Triangle_Normal &T, Vector3d v)=0;
	virtual Vector3cd weak1(Triangle_All &T, Vector3d v) = 0;
	virtual Vector3cd weak3(Triangle_Normal &T, Vector3d v) = 0;
	virtual Vector3cd weak3(Triangle_All &T, Vector3d v) = 0;
	virtual Vector3cd weak4(Vector2d &p1, Vector2d &p2, Vector2d &p3, double v1, double v2, Vector2d &p6) = 0;
	virtual void weak5(double x1, double y1, double x2, double y2, double v1, double v2, vector<vector<double>> &nbound, vector<complex<double>> &out5m, complex<double> &out5) = 0;
	virtual VectorXcd setRightHand(TriangleMesh &U, TriangleMesh &L, vector<int> &nu)=0;
	virtual SparseMatrix<complex<double>> setA(double nn, vector<int> &nu, vector<vector<double>> &nbound, vector<vector<gridCell>> &grid) = 0;
	virtual mxArray* setA_mx(double nn, vector<int> &nu, vector<vector<double>> &nbound, vector<vector<gridCell>> &grid) = 0;
	virtual VectorXcd setB(vector<vector<gridCell>> &grid, VectorXcd &rh, int nn, vector<int> &nu) = 0;
	virtual mxArray* setB_mx(vector<vector<gridCell>> &grid, VectorXcd &rh, int nn, vector<int> &nu) = 0;
	VectorXcd solveX(SparseMatrix<complex<double>> &A, VectorXcd &B);
	VectorXcd solveX_mx(mxArray *A, mxArray *B);
	virtual void assign(VectorXcd u, TriangleMesh &U, TriangleMesh &L, vector<int> &nu) = 0;
	virtual VectorXcd getAperture(vector<vector<double>> &nbound, TriangleMesh &U, TriangleMesh &L) = 0;
	virtual void drawAperture(VectorXd &plotX, VectorXd &plotY, vector<vector<double>> &nbound, int sign) = 0;

	//Auxiliary function ������������
	virtual void checkVirtualBorder();
	double findzero(double x1, double y1, double x2, double y2);
	void initializeComplexVector(vector<complex<double>> &vector, complex<double> InitialValue, int vectorSize); //��ʼ����������
	void complexVectorAdd(vector<complex<double>> &vector1, vector<complex<double>> vector2); //���������ӷ�
	double calculateTriangleArea(Vector3d tri1, Vector3d tri2); //�������������
	vector<complex<double>> eigenVector2stdVector(VectorXcd eigenVector);
	int findXinNbound(vector<vector<double>> nbound, double x);
	int findXIndexinNbound(vector<vector<double>> nbound, double xIndex);
	SparseMatrix<complex<double>> buildSparseMatrixA(double* pA, int MatrixSize, int ArowNum);
	int sign(double num);

	/*
	���ƿھ��洦����ֵ��
	sign ��ʶ������ֵ������ͣ�
	sign = 0   ��ֵ
	sign = 1   ʵ��
	sign = -1  �鲿
	sign = 10 ��λ
	*/
	void PlotAperture(string title, string xlabel, string ylabel, int sign);
private:
	
};
