// Test.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <time.h>

#define _USE_MATH_DEFINES //若不使用此宏，math.h中无法找到π的定义M_PI
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
#include <InhomogeneousThreeCircleRectangleCavity.h>
#include <InhomogeneousThreeCircleCircleCavity.h>
#include <ThreeRectangleCavity.h>

using namespace Eigen;
using namespace std;

typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix <double> SparseMatrixType;
typedef Eigen::SparseMatrix <complex<double>> SparseComplexMatrixType;

int main()
{
	
	#pragma region 单个矩形腔体

	////对应matlab程序：F:\计算数学\github-单个规则矩形腔体问题\A-fast-algorithm-for-the-Singal-Cavity-Scattering\PDE-Solver-SingleCavity - 副本
	////选择腔体类型
	//SingalRectangleCavity cavity(1);

	/////初始化参数
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

	//// 绘制三角形网格
	//cavity.PlotTriangleMesh("title", "xlabel", "ylabel");

	////求解-并绘制口径面处的解
	//cavity.SolveAperture("title", "xlabel", "ylabel", 0);

	////求解-并绘制整个腔体内的数值解
	//cavity.SolveCavity("title", "xlabel", "ylabel");

	////// 计算并绘制RCS
	////double interval = 0.25;
	////cavity.SolveRCS(interval);

	#pragma endregion

	#pragma region 单个非均匀介质矩形腔体

	#pragma region AAMM1
	
	////选择腔体类型
	//InhomogeneousSingalRectangleCavityAAMM1 cavity(1);

	/////初始化参数
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

	//// 绘制三角形网格
	//cavity.PlotTriangleMesh("title", "xlabel", "ylabel");

	////求解-并绘制口径面处的解
	//cavity.SolveAperture("title", "xlabel", "ylabel", 0);

	////求解-并绘制整个腔体内的数值解
	//cavity.SolveCavity("title", "xlabel", "ylabel");

	////double interval = 0.25;
	////cavity.SolveRCS(interval);

	#pragma endregion

	#pragma region AAMM2
	
	////选择腔体类型
	//InhomogeneousSingalRectangleCavityAAMM2 cavity(1);

	/////初始化参数
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

	//// 绘制三角形网格
	//cavity.PlotTriangleMesh("title", "xlabel", "ylabel");

	////求解-并绘制口径面处的解
	//cavity.SolveAperture("title", "xlabel", "ylabel", 0);

	////求解-并绘制整个腔体内的数值解
	//cavity.SolveCavity("title", "xlabel", "ylabel");

	/////*double interval = 0.25;
	////cavity.SolveRCS(interval);*/

	#pragma endregion

	#pragma endregion


	#pragma region 三个非均匀介质腔体

	#pragma region 矩形腔体+圆形非均匀介质
	
	//InhomogeneousThreeRectangleCircleCavity cavity(1);
	/////初始化参数
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

	//// 绘制三角形网格
	//cavity.PlotTriangleMesh("title", "xlabel", "ylabel");

	////求解-并绘制口径面处的解
	//cavity.SolveAperture("title", "xlabel", "ylabel", 0);

	////求解-并绘制整个腔体内的数值解
	//cavity.SolveCavity("title", "xlabel", "ylabel");
	
	////求解-并绘制RCS
	//double interval = 0.5;
	//cavity.SolveRCS(interval);
	#pragma endregion

	#pragma region 矩形腔体+矩形非均匀介质

	//InhomogeneousThreeRectangleRectangleCavity cavity(1);
	/////初始化参数
	//double VirtualTop = 0;
	//double VirtualBottom = -1.2;
	//double VirtualLeft = -3;
	//double VirtualRight = 3;
	//cavity.InitVirtualBorder(VirtualTop, VirtualBottom, VirtualLeft, VirtualRight);

	////int m = 48;
	////int n = 12;
	//int m = 96;
	//int n = 24;
	////int m = 192;
	////int n = 48;
	///*int m = 384;
	//int n = 96;*/
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

	//// 绘制三角形网格
	//cavity.PlotTriangleMesh("title", "xlabel", "ylabel");

	////求解-并绘制口径面处的解
	//cavity.SolveAperture("title", "xlabel", "ylabel", 0);

	////求解-并绘制整个腔体内的数值解
	//cavity.SolveCavity("title", "xlabel", "ylabel");
	
	////求解-并绘制RCS
	//double interval = 0.5;
	//cavity.SolveRCS(interval);

	#pragma endregion

	#pragma region 圆形腔体+矩形非均匀介质

	InhomogeneousThreeCircleRectangleCavity cavity(1);
	///初始化参数
	double VirtualTop = 0;
	double VirtualBottom = -0.6;
	double VirtualLeft = -3;
	double VirtualRight = 3;
	cavity.InitVirtualBorder(VirtualTop, VirtualBottom, VirtualLeft, VirtualRight);

	//int m = 48;
	//int n = 6;
	//int m = 192;
	//int n = 24;
	int m = 2400;
	int n = 300;
	//int m = 3360;
	//int n = 420;
	cavity.InitMesh(m, n);

	double k0 = 16 * M_PI;
	complex<double> epr1(1, 0);
	complex<double> epr1_int(1, 1);
	complex<double> epr2(1, 0);
	complex<double> epr2_int(1, 1);
	complex<double> epr3(1, 0);
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

	cavityCircle circle1;
	circle1.circleCenter = Vector2d(-2, 0);
	circle1.radius = 0.5;
	cavityCircle circle2;
	circle2.circleCenter = Vector2d(0, 0);
	circle2.radius = 0.5;
	cavityCircle circle3;
	circle3.circleCenter = Vector2d(2, 0);
	circle3.radius = 0.5;
	cavity.InitCavityShapeParameter(circle1, circle2, circle3);

	cavityBox box1;
	box1.cavityTop = -0.1;
	box1.cavityBottom = -0.3;
	box1.cavityLeft = -2.25;
	box1.cavityRight = -1.75;
	cavityBox box2;
	box2.cavityTop = -0.1;
	box2.cavityBottom = -0.2;
	box2.cavityLeft = -0.375;
	box2.cavityRight = 0.375;
	cavityBox box3;
	box3.cavityTop = -0.1;
	box3.cavityBottom = -0.3;
	box3.cavityLeft = 1.75;
	box3.cavityRight = 2.25;
	cavity.InitCavityInhomogeneousShapeParameter(box1, box2, box3);

	//// 绘制三角形网格
	//cavity.PlotTriangleMesh("title", "xlabel", "ylabel");

	////求解-并绘制口径面处的解
	//cavity.SolveAperture("title", "xlabel", "ylabel", 0);

	////求解-并绘制口径面处的解
	//cavity.SolveAperture("title", "xlabel", "ylabel", 0);

	////求解-并绘制整个腔体内的数值解
	//cavity.SolveCavity("title", "xlabel", "ylabel");

	//求解-并绘制RCS
	double interval = 0.5;
	cavity.SolveRCS(interval);

	#pragma endregion

	#pragma region 圆形腔体+圆形非均匀介质

	//InhomogeneousThreeCircleCircleCavity cavity(1);
	/////初始化参数
	//double VirtualTop = 0;
	//double VirtualBottom = -0.75;
	//double VirtualLeft = -3;
	//double VirtualRight = 3;
	//cavity.InitVirtualBorder(VirtualTop, VirtualBottom, VirtualLeft, VirtualRight);

	////int m = 48;
	////int n = 6;
	//int m = 192;
	//int n = 24;
	////int m = 192;
	////int n = 96;
	////int m = 2400;
	////int n = 300;
	////int m = 3360;
	////int n = 420;
	//cavity.InitMesh(m, n);

	//double k0 = 8 * M_PI;
	//complex<double> epr1(1, 0);
	//complex<double> epr1_int(1, 1);
	//complex<double> epr2(1, 0);
	//complex<double> epr2_int(4, 1);
	//complex<double> epr3(1, 0);
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

	//cavityCircle circle1;
	//circle1.circleCenter = Vector2d(-2, 0);
	//circle1.radius = 0.5;
	//cavityCircle circle2;
	//circle2.circleCenter = Vector2d(0, 0);
	//circle2.radius = 0.5;
	//cavityCircle circle3;
	//circle3.circleCenter = Vector2d(2, 0);
	//circle3.radius = 0.5;
	//cavity.InitCavityShapeParameter(circle1, circle2, circle3);

	//cavityCircle circle1_int;
	//circle1_int.circleCenter = Vector2d(-2, -1);
	//circle1_int.radius = 0.6;
	//cavityCircle circle2_int;
	//circle2_int.circleCenter = Vector2d(0, -1);
	//circle2_int.radius = 0.8;
	//cavityCircle circle3_int;
	//circle3_int.circleCenter = Vector2d(2, -1);
	//circle3_int.radius = 0.6;
	//cavity.InitCavityInhomogeneousShapeParameter(circle1_int, circle2_int, circle3_int);

	//// 绘制三角形网格
	//cavity.PlotTriangleMesh("title", "xlabel", "ylabel");

	////求解-并绘制口径面处的解
	//cavity.SolveAperture("title", "xlabel", "ylabel", 0);

	////求解-并绘制整个腔体内的数值解
	//cavity.SolveCavity("title", "xlabel", "ylabel");

	//// 计算并绘制RCS
	//double interval = 0.5;
	//cavity.SolveRCS(interval);

	#pragma endregion

	#pragma endregion


	#pragma region 三个均匀介质腔体

	#pragma region 矩形均匀介质多腔体
	///*
	//由于均匀介质多腔体使用非均匀介质的架构实现
	//内部暗含的逻辑是腔体具有两种介质，我们调用时人为保证这两种介质相同即可
	//*/

	//ThreeRectangleCavity cavity(1);
	/////初始化参数
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

	//// 绘制三角形网格
	//cavity.PlotTriangleMesh("title", "xlabel", "ylabel");

	////求解-并绘制口径面处的解
	//cavity.SolveAperture("title", "xlabel", "ylabel", 0);

	////求解-并绘制整个腔体内的数值解
	//cavity.SolveCavity("title", "xlabel", "ylabel");

	//// 计算并绘制RCS
	//double interval = 0.5;
	//cavity.SolveRCS(interval);

	#pragma endregion

	#pragma endregion


	system("PAUSE");
    return 0;
}

