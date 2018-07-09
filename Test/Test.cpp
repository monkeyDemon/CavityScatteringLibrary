// Test.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include <iostream>
#include <time.h>

#define _USE_MATH_DEFINES //����ʹ�ô˺꣬math.h���޷��ҵ��еĶ���M_PI
#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Cavity.h>
#include <SingalCavity.h>
#include <SingalRectangleCavity.h>
#include <InhomogeneousSingalRectangleCavity.h>
#include <InhomogeneousSingalRectangleCavityAAMM1.h>
#include <InhomogeneousSingalRectangleCavityAAMM2.h>
#include <InhomogeneousThreeRectangleCircleCavity.h>
#include <InhomogeneousThreeRectangleRectangleCavity.h>
#include <ThreeRectangleCavity.h>

using namespace Eigen;
using namespace std;

typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix <double> SparseMatrixType;
typedef Eigen::SparseMatrix <complex<double>> SparseComplexMatrixType;

int main()
{
	
	#pragma region ��������ǻ��

	////��Ӧmatlab����F:\������ѧ\github-�����������ǻ������\A-fast-algorithm-for-the-Singal-Cavity-Scattering\PDE-Solver-SingleCavity - ����
	////ѡ��ǻ������
	//SingalRectangleCavity cavity(1);

	/////��ʼ������
	//double VirtualTop = 1;
	//double VirtualBottom = -1;
	//double VirtualLeft = -1;
	//double VirtualRight = 1;
	//cavity.InitVirtualBorder(VirtualTop, VirtualBottom, VirtualLeft, VirtualRight);

	//int m = 16;
	//int n = 16;
	//cavity.InitMesh(m, n);

	//double k0 = 2 * M_PI;
	//complex<double> epr(4, 1);
	//double theta = 0;
	//cavity.InitElectromagneticParameter(k0, epr, theta);

	//double apertureLeft = -0.5;
	//double apertureRight = 0.5;
	//double apertureY = 1;
	//cavity.InitAperture(apertureLeft, apertureRight, apertureY);

	//double cavityBottom = 0.75;
	//cavity.InitCavityShapeParameter(cavityBottom);

	//// ��������������
	//cavity.PlotTriangleMesh("title", "xlabel", "ylabel");

	////// ���
	////cavity.Solve();//���
	////
	////// �ھ����Ŀ��ӻ�
	////cavity.PlotAperture("title", "xlabel", "ylabel", 0);

	////// ���㲢����RCS
	////double interval = 0.25;
	////cavity.SolveRCS(interval);

	#pragma endregion



	#pragma region �����Ǿ��Ƚ��ʾ���ǻ��

	#pragma region AAMM1
	
	////ѡ��ǻ������
	//InhomogeneousSingalRectangleCavityAAMM1 cavity(1);

	/////��ʼ������
	//double VirtualTop = 0;
	//double VirtualBottom = -2;
	//double VirtualLeft = -0.5;
	//double VirtualRight = 1.5;
	//cavity.InitVirtualBorder(VirtualTop, VirtualBottom, VirtualLeft, VirtualRight);

	//int m = 24;
	//int n = 24;
	//cavity.InitMesh(m, n);

	//double k0 = 4 * M_PI;
	//complex<double> epr(4, 2);
	//complex<double> epr_int(2.39, 1.84);
	//double theta = 0;
	//cavity.InitElectromagneticParameter(k0, epr, epr_int, theta);

	//double apertureLeft = 0;
	//double apertureRight = 1;
	//double apertureY = 0;
	//cavity.InitAperture(apertureLeft, apertureRight, apertureY);

	//double cavityTop = 0;
	//double cavityBottom = -1;
	//double cavityLeft = 0;
	//double cavityRight = 1;
	//double cavityTop_ = -0.5;
	//double cavityBottom_ = -1;
	//double cavityLeft_ = 0;
	//double cavityRight_ = 1;
	//cavity.InitCavityShapeParameter(cavityTop, cavityBottom, cavityLeft, cavityRight, cavityTop_, cavityBottom_, cavityLeft_, cavityRight_);

	//// ��������������
	//cavity.PlotTriangleMesh("title", "xlabel", "ylabel");

	///////���
	////cavity.Solve();//���

	///////���ӻ�
	////cavity.PlotAperture("title", "xlabel", "ylabel", 0);

	////double interval = 0.25;
	////cavity.SolveRCS(interval);

	#pragma endregion

	#pragma region AAMM2
	
	////ѡ��ǻ������
	//InhomogeneousSingalRectangleCavityAAMM2 cavity(1);

	/////��ʼ������
	//double VirtualTop = 1;
	//double VirtualBottom = -1;
	//double VirtualLeft = -0.5;
	//double VirtualRight = 1.5;
	//cavity.InitVirtualBorder(VirtualTop, VirtualBottom, VirtualLeft, VirtualRight);

	//int m = 24;
	//int n = 24;
	//cavity.InitMesh(m, n);

	//double k0 = 4 * M_PI;
	//complex<double> epr(4, 2);
	//complex<double> epr_int(2.39, 1.84);
	//double theta = 0;
	//cavity.InitElectromagneticParameter(k0, epr, epr_int, theta);

	//double apertureLeft = 0;
	//double apertureRight = 1;
	//double apertureY = 1;
	//cavity.InitAperture(apertureLeft, apertureRight, apertureY);

	//double cavityTop = 1;
	//double cavityBottom = 0;
	//double cavityLeft = 0;
	//double cavityRight = 1;
	//double cavityTop_ = 0.125;
	//double cavityBottom_ = 0;
	//double cavityLeft_ = 0.25;
	//double cavityRight_ = 0.75;
	//cavity.InitCavityShapeParameter(cavityTop, cavityBottom, cavityLeft, cavityRight, cavityTop_, cavityBottom_, cavityLeft_, cavityRight_);

	//// ��������������
	//cavity.PlotTriangleMesh("title", "xlabel", "ylabel");

	///////���
	////cavity.Solve();//���

	///////���ӻ�
	////cavity.PlotAperture("title", "xlabel", "ylabel", 0);

	/////*double interval = 0.25;
	////cavity.SolveRCS(interval);*/

	#pragma endregion

	#pragma endregion

	#pragma region �����Ǿ��Ƚ���ǻ��

	#pragma region ����ǻ��+Բ�ηǾ��Ƚ���
	
	//InhomogeneousThreeRectangleCircleCavity cavity(1);
	/////��ʼ������
	//double VirtualTop = 0;
	//double VirtualBottom = -1.5;
	//double VirtualLeft = -3;
	//double VirtualRight = 3;
	//cavity.InitVirtualBorder(VirtualTop, VirtualBottom, VirtualLeft, VirtualRight);

	//int m = 48;
	//int n = 12;
	////int m = 192;
	////int n = 96;
	///*int m = 384;
	//int n = 192;*/
	////int m = 24;
	////int n = 6;
	//cavity.InitMesh(m, n);

	//double k0 = 64 * M_PI;
	//complex<double> epr1(1, 0);
	//complex<double> epr1_int(1, 0);
	//complex<double> epr2(1, 0);
	//complex<double> epr2_int(1, 0);
	//complex<double> epr3(1, 1);
	//complex<double> epr3_int(1, 1);
	//double theta = 0;
	//cavity.InitElectromagneticParameter(k0, epr1, epr1_int, epr2, epr2_int, epr3, epr3_int, theta);

	//double aperture1Left = -2.5;
	//double aperture1Right = -1.5;
	//double aperture2Left = -0.5;
	//double aperture2Right = 0.5;
	//double aperture3Left = 1.5;
	//double aperture3Right = 2.5;
	//double apertureY = 0;
	//cavity.InitAperture(apertureY, aperture1Left, aperture1Right, aperture2Left, aperture2Right, aperture3Left, aperture3Right);

	//cavityBox box1;
	//box1.cavityTop = 0;
	//box1.cavityBottom = -1;
	//box1.cavityLeft = -2.5;
	//box1.cavityRight = -1.5;
	//cavityBox box2;
	//box2.cavityTop = 0;
	//box2.cavityBottom = -1;
	//box2.cavityLeft = -0.5;
	//box2.cavityRight = 0.5;
	//cavityBox box3;
	//box3.cavityTop = 0;
	//box3.cavityBottom = -1;
	//box3.cavityLeft = 1.5;
	//box3.cavityRight = 2.5;
	//cavity.InitCavityShapeParameter(box1, box2, box3);

	//cavityCircle circle1;
	//circle1.circleCenter = Vector2d(-2, -1);
	//circle1.radius = 0.8;
	//cavityCircle circle2;
	//circle2.circleCenter = Vector2d(0, -1);
	//circle2.radius = 0.8; 
	//cavityCircle circle3;
	//circle3.circleCenter = Vector2d(2, -1);
	//circle3.radius = 0.8;
	//cavity.InitCavityInhomogeneousShapeParameter(circle1, circle2, circle3);

	//// ��������������
	//cavity.PlotTriangleMesh("title", "xlabel", "ylabel");

	///////���
	////cavity.Solve();//���

	///////���ӻ�
	////cavity.PlotAperture("title", "xlabel", "ylabel", 0);

	#pragma endregion

	#pragma region ����ǻ��+���ηǾ��Ƚ���

	InhomogeneousThreeRectangleRectangleCavity cavity(1);
	///��ʼ������
	double VirtualTop = 0;
	double VirtualBottom = -1.2;
	double VirtualLeft = -3;
	double VirtualRight = 3;
	cavity.InitVirtualBorder(VirtualTop, VirtualBottom, VirtualLeft, VirtualRight);

	//int m = 48;
	//int n = 12;
	int m = 96;
	int n = 24;
	//int m = 192;
	//int n = 48;
	/*int m = 384;
	int n = 96;*/
	cavity.InitMesh(m, n);

	double k0 = 4 * M_PI;
	complex<double> epr1(1, 0);
	complex<double> epr1_int(1, 0);
	complex<double> epr2(1, 0);
	complex<double> epr2_int(1, 0);
	complex<double> epr3(1, 1);
	complex<double> epr3_int(1, 1);
	double theta = 0;
	cavity.InitElectromagneticParameter(k0, epr1, epr1_int, epr2, epr2_int, epr3, epr3_int, theta);

	double aperture1Left = -2.5;
	double aperture1Right = -1.5;
	double aperture2Left = -0.5;
	double aperture2Right = 0.5;
	double aperture3Left = 1.5;
	double aperture3Right = 2.5;
	double apertureY = 0;
	cavity.InitAperture(apertureY, aperture1Left, aperture1Right, aperture2Left, aperture2Right, aperture3Left, aperture3Right);

	cavityBox box1;
	box1.cavityTop = 0;
	box1.cavityBottom = -1;
	box1.cavityLeft = -2.5;
	box1.cavityRight = -1.5;
	cavityBox box2;
	box2.cavityTop = 0;
	box2.cavityBottom = -1;
	box2.cavityLeft = -0.5;
	box2.cavityRight = 0.5;
	cavityBox box3;
	box3.cavityTop = 0;
	box3.cavityBottom = -1;
	box3.cavityLeft = 1.5;
	box3.cavityRight = 2.5;
	cavity.InitCavityShapeParameter(box1, box2, box3);

	cavityBox box1_int;
	box1_int.cavityTop = -0.8;
	box1_int.cavityBottom = -1;
	box1_int.cavityLeft = -2.25;
	box1_int.cavityRight = -1.75;
	cavityBox box2_int;
	box2_int.cavityTop = -0.8;
	box2_int.cavityBottom = -1;
	box2_int.cavityLeft = -0.25;
	box2_int.cavityRight = 0.25;
	cavityBox box3_int;
	box3_int.cavityTop = -0.8;
	box3_int.cavityBottom = -1;
	box3_int.cavityLeft = 1.75;
	box3_int.cavityRight = 2.25;
	cavity.InitCavityInhomogeneousShapeParameter(box1_int, box2_int, box3_int);

	// ��������������
	cavity.PlotTriangleMesh("title", "xlabel", "ylabel");

	/////���
	//cavity.Solve();//���

	/////���ӻ�
	//cavity.PlotAperture("title", "xlabel", "ylabel", 0);

	#pragma endregion


	#pragma endregion


	#pragma region �������Ƚ���ǻ��

	#pragma region ���ξ��Ƚ��ʶ�ǻ��
	///*
	//���ھ��Ƚ��ʶ�ǻ��ʹ�÷Ǿ��Ƚ��ʵļܹ�ʵ��
	//�ڲ��������߼���ǻ��������ֽ��ʣ����ǵ���ʱ��Ϊ��֤�����ֽ�����ͬ����
	//*/

	//ThreeRectangleCavity cavity(1);
	/////��ʼ������
	//double VirtualTop = 0;
	//double VirtualBottom = -1.5;
	//double VirtualLeft = -3;
	//double VirtualRight = 3;
	//cavity.InitVirtualBorder(VirtualTop, VirtualBottom, VirtualLeft, VirtualRight);

	//int m = 48;
	//int n = 24;
	////int m = 192;
	////int n = 96;
	///*int m = 384;
	//int n = 192;*/
	////int m = 24;
	////int n = 12;
	//cavity.InitMesh(m, n);

	//double k0 = 4 * M_PI;
	//complex<double> epr1(1, 0);
	//complex<double> epr1_int(1, 0);
	//complex<double> epr2(1, 0);
	//complex<double> epr2_int(1, 0);
	//complex<double> epr3(1, 1);
	//complex<double> epr3_int(1, 1);
	//double theta = 0;
	//cavity.InitElectromagneticParameter(k0, epr1, epr1_int, epr2, epr2_int, epr3, epr3_int, theta);

	//double aperture1Left = -2.5;
	//double aperture1Right = -1.5;
	//double aperture2Left = -0.5;
	//double aperture2Right = 0.5;
	//double aperture3Left = 1.5;
	//double aperture3Right = 2.5;
	//double apertureY = 0;
	//cavity.InitAperture(apertureY, aperture1Left, aperture1Right, aperture2Left, aperture2Right, aperture3Left, aperture3Right);

	//cavityBox box1;
	//box1.cavityTop = 0;
	//box1.cavityBottom = -1;
	//box1.cavityLeft = -2.5;
	//box1.cavityRight = -1.5;
	//cavityBox box2;
	//box2.cavityTop = 0;
	//box2.cavityBottom = -1;
	//box2.cavityLeft = -0.5;
	//box2.cavityRight = 0.5;
	//cavityBox box3;
	//box3.cavityTop = 0;
	//box3.cavityBottom = -1;
	//box3.cavityLeft = 1.5;
	//box3.cavityRight = 2.5;
	//cavity.InitCavityShapeParameter(box1, box2, box3);

	//cavityBox box1_int;
	//box1_int.cavityTop = -0.8;
	//box1_int.cavityBottom = -1;
	//box1_int.cavityLeft = -2.25;
	//box1_int.cavityRight = -1.75;
	//cavityBox box2_int;
	//box2_int.cavityTop = -0.8;
	//box2_int.cavityBottom = -1;
	//box2_int.cavityLeft = -0.25;
	//box2_int.cavityRight = 0.25;
	//cavityBox box3_int;
	//box3_int.cavityTop = -0.8;
	//box3_int.cavityBottom = -1;
	//box3_int.cavityLeft = 1.75;
	//box3_int.cavityRight = 2.25;
	//cavity.InitCavityInhomogeneousShapeParameter(box1_int, box2_int, box3_int);

	/////���
	//cavity.Solve();//���

	/////���ӻ�
	//cavity.PlotAperture("title", "xlabel", "ylabel", 0);

	#pragma endregion

	#pragma endregion


	system("PAUSE");
    return 0;
}

