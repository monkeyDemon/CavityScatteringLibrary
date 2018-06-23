#pragma once
#include "Cavity.h"

class _declspec(dllexport) SingalCavity : public Cavity
{
public:
	SingalCavity(unsigned int cavityType);
	//初始化电磁相关物理参数
	void InitElectromagneticParameter(double k0, complex<double> epr, double theta);
	//初始化口径面信息
	void InitAperture(double apertureLeft, double apertureRight, double apertureY);
	
	bool Solve();

protected:
	complex<double> epr;
	complex<double> k2;
	double apertureLeft;
	double apertureRight;
	double apertureY;

	MatrixXcd G_aperture;
	VectorXcd g_aperture;

	Matrix2d beta(double x, double y);
	complex<double> q(double x, double y);
	complex<double> f(double x, double y);
	complex<double> u(double x, double y);

	void setTri(TriangleMesh &U, TriangleMesh &L, int &nn, vector<int> &nu, vector<vector<double>> &nbound);
	VectorXcd compute_g(MatrixXcd &G, vector<vector<double>> &nbound);
	vector<vector<gridCell>> setGrid(TriangleMesh &U, TriangleMesh &L, vector<vector<double>> &nbound, vector<int> &nu);
	void weak(Triangle_Normal &T, Vector3d v, int topsign, vector<vector<double>> &nbound, vector <complex<double>> &out, vector <complex<double>> &out5m);
	void weak(Triangle_All &T, Vector3d v, int topsign, vector<vector<double>> &nbound, vector <complex<double>> &out, vector <complex<double>> &out5m);
	Vector3cd weak1(Triangle_Normal &T, Vector3d v);
	Vector3cd weak1(Triangle_All &T, Vector3d v);
	Vector3cd weak3(Triangle_Normal &T, Vector3d v);
	Vector3cd weak3(Triangle_All &T, Vector3d v);
	Vector3cd weak4(Vector2d &p1, Vector2d &p2, Vector2d &p3, double v1, double v2, Vector2d &p6);
	void weak5(double x1, double y1, double x2, double y2, double v1, double v2, vector<vector<double>> &nbound, vector<complex<double>> &out5m, complex<double> &out5);
	VectorXcd setRightHand(TriangleMesh &U, TriangleMesh &L, vector<int> &nu);
	SparseMatrix<complex<double>> setA(double nn, vector<int> &nu, vector<vector<double>> &nbound, vector<vector<gridCell>> &grid);
	mxArray* setA_mx(double nn, vector<int> &nu, vector<vector<double>> &nbound, vector<vector<gridCell>> &grid);
	VectorXcd setB(vector<vector<gridCell>> &grid, VectorXcd &rh, int nn, vector<int> &nu);
	mxArray* setB_mx(vector<vector<gridCell>> &grid, VectorXcd &rh, int nn, vector<int> &nu);
	void assign(VectorXcd u, TriangleMesh &U, TriangleMesh &L, vector<int> &nu);
	VectorXcd getAperture(vector<vector<double>> &nbound, TriangleMesh &U, TriangleMesh &L);
	void drawAperture(VectorXd &plotX, VectorXd &plotY, vector<vector<double>> &nbound, int sign);
private:
	//setTri中需要使用的子函数
	void setTriangleType(double x1, double x2, double x3, double y1, double y2, double y3, int &sgn);
	void tri(double x1, double x2, double x3, double y1, double y2, double y3, Triangle_Normal& l);
	void tri(double x1, double x2, double x3, double y1, double y2, double y3, Triangle_All& l);
	void solveTri(Triangle_All &l);
};
