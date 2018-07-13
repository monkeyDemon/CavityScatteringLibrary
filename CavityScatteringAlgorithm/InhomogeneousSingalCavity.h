#pragma once
#include "Cavity.h"

class  _declspec(dllexport) InhomogeneousSingalCavity : public Cavity
{
public:
	InhomogeneousSingalCavity(unsigned int cavityType);
	//初始化电磁相关物理参数
	void InitElectromagneticParameter(double k0, complex<double>epr, complex<double> epr_int, double theta);
	//初始化口径面信息
	void InitAperture(double apertureLeft, double apertureRight, double apertureY);
	
	bool Solve();
	bool SolveAperture(string title, string xlabel, string ylabel, int sign);

	void PlotTriangleMesh(string title, string xlabel, string ylabel);

protected:
	complex<double> epr, epr_int;
	complex<double> k2, k2_int;

	double apertureLeft;
	double apertureRight;
	double apertureY;


	MatrixXcd G_aperture;
	VectorXcd g_aperture;

	virtual Matrix2d beta(double x, double y);
	virtual Matrix2d beta_(double x, double y);
	virtual complex<double> q(double x, double y);
	virtual complex<double> q_(double x, double y);
	virtual complex<double> f(double x, double y);
	virtual complex<double> f_(double x, double y);
	virtual complex<double> u(double x, double y);
	virtual complex<double> u_(double x, double y);
	virtual double a(double x, double y);
	virtual Vector2d b(double x, double y);
	virtual double phi_int(double x, double y) =0;

	void setTri(TriangleMesh &U, TriangleMesh &L, int &nn, vector<int> &nu, vector<vector<double>> &nbound);
	VectorXcd compute_g(MatrixXcd &G, vector<vector<double>> &nbound);
	vector<vector<gridCell>> setGrid(TriangleMesh &U, TriangleMesh &L, vector<vector<double>> &nbound, vector<int> &nu);
	void weak(Triangle_Normal &T, Vector3d v, int topsign, vector<vector<double>> &nbound, vector <complex<double>> &out, vector <complex<double>> &out5m);
	void weak(Triangle_All &T, Vector3d v, int topsign, vector<vector<double>> &nbound, vector <complex<double>> &out, vector <complex<double>> &out5m);
	Vector4cd weak_int(Triangle_Normal &T, Vector3d v);
	Vector4cd weak_int(Triangle_All &T, Vector3d v);
	Vector3cd weak1(Triangle_Normal &T, Vector3d v);
	Vector3cd weak1(Triangle_All &T, Vector3d v);
	Vector4cd weak1_int(Triangle_Normal &T, Vector3d v,bool beta);
	Vector4cd weak1_int(Triangle_All &T, Vector3d v, bool beta);
	Vector3cd weak3(Triangle_Normal &T, Vector3d v);
	Vector3cd weak3(Triangle_All &T, Vector3d v);
	Vector4cd weak3_int(Triangle_Normal &T, Vector3d v, bool q);
	Vector4cd weak3_int(Triangle_All &T, Vector3d v, bool q);
	Vector3cd weak4(Vector2d &p1, Vector2d &p2, Vector2d &p3, double v1, double v2, Vector2d &p6);
	Vector3cd weak4(Triangle_All &T, double v2, double v3, bool beta);
	double weak4_int(Vector2d p1, Vector2d p2, Vector2d b, double v1, double v2);
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
	void setTriangleType(double x1, double x2, double x3, double y1, double y2, double y3, int &sgn, int &sgn_int, int &sgn_intex);
	int setTriangleType_int(double x1, double x2, double x3, double y1, double y2, double y3);
	int setTriangleType_intex(double x1, double x2, double x3, double y1, double y2, double y3);
	void tri(double x1, double x2, double x3, double y1, double y2, double y3, Triangle_All& l);
	void tri_int1(double x1, double  x2, double  x3, double y1, double y2, double y3, Triangle_All &l);
	void tri(double x1, double x2, double x3, double y1, double y2, double y3, Triangle_Normal& l);
	//void tri(double a1, double a2, double a3, double x1, double x2, double x3, double y1, double y2, double y3, Triangle& l);
	//void tri_int(double a1, double a2, double a3, double x1, double  x2, double  x3, double y1, double y2, double y3, Triangle &l);
	//void tri_int1(double a1, double a2, double a3, double x1, double  x2, double  x3, double y1, double y2, double y3, Triangle &l);
	void solveTri(Triangle_All &l);
	void solveTri_int(Triangle_All &l);
	void solveTri_int1(Triangle_All &l);
	double findzero_int(double x1, double y1, double x2, double y2);
	void ustar_int(Triangle_All &T, Vector4d &u4, Vector4d &u5, Vector4d &u4_, Vector4d &u5_);
	void ustar2_int(Triangle_All &T, Vector4d &u5, Vector4d &u5_);
	void ustar_int1(Triangle_All &T, Vector4d &u4, Vector4d &u5, Vector4d &u4_, Vector4d &u5_);
	void ustar2_int1(Triangle_All &T, Vector4d &u5, Vector4d &u5_);

	//setGrid中需要使用的子函数
	complex<double> value(double x, double y);
	void reloadtri(Triangle_All &T2, Triangle_All &T);
	void reloadtriex(Triangle_All &T3, Triangle_All &T);
	Vector4cd solveweak4(Triangle_All &T, double v2, double v3);
	
	//setRightHand中需要使用的子函数
	complex<double> fInt_f(Triangle_Normal &T, Vector3d v);
	complex<double> fInt_f(Triangle_All &T, Vector3d v);
	complex<double> fInt_int(Triangle_Normal &T, Vector3d v);
	complex<double> fInt_int(Triangle_All &T, Vector3d v);
	complex<double> tr(Triangle_Normal &T, Vector3d v, bool f);
	complex<double> tr(Triangle_All &T, Vector3d v, bool f);

	//assign中需要使用的子函数
	void assign_int(Triangle_All &T);
	void assign_intex(Triangle_All &T);
};

