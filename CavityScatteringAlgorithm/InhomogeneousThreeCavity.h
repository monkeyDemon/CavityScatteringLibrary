#pragma once
#include "Cavity.h"

class  _declspec(dllexport) InhomogeneousThreeCavity : public Cavity
{
public:
	InhomogeneousThreeCavity(unsigned int cavityType);
	//��ʼ���������������
	void InitElectromagneticParameter(double k0, complex<double>epr1, complex<double> epr1_int, complex<double>epr2, complex<double> epr2_int, complex<double>epr3, complex<double> epr3_int, double theta);
	//��ʼ���ھ�����Ϣ
	void InitAperture(double apertureY, double aperture1Left, double aperture1Right, double aperture2Left, double aperture2Right, double aperture3Left, double aperture3Right);

	bool Solve();

protected:
	double k0;
	double theta;
	complex<double> epr1, epr1_int;
	complex<double> epr2, epr2_int;
	complex<double> epr3, epr3_int;
	complex<double> k1_2, k1_2int;
	complex<double> k2_2, k2_2int;
	complex<double> k3_2, k3_2int;

	int cavity_number;

	double apertureY;

	double aperture1Left;
	double aperture1Right;
	double aperture2Left;
	double aperture2Right;
	double aperture3Left;
	double aperture3Right;

	double separator1_2;
	double separator2_3;


	MatrixXcd G1_aperture;
	MatrixXcd G2_aperture;
	MatrixXcd G3_aperture;
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
	virtual double phi_int(double x, double y) = 0;

	void setTri(TriangleMesh &U, TriangleMesh &L, int &nn, vector<int> &nu, vector<vector<double>> &nbound);
	VectorXcd compute_g(MatrixXcd &G, vector<vector<double>> &nbound);
	vector<vector<gridCell>> setGrid(TriangleMesh &U, TriangleMesh &L, vector<vector<double>> &nbound, vector<int> &nu);
	void weak(Triangle &T, Vector3d v, int topsign, vector<vector<double>> &nbound, vector <complex<double>> &out, vector <complex<double>> &out5m);
	Vector4cd weak_int(Triangle &T, Vector3d v);
	Vector3cd weak1(Triangle &T, Vector3d v);
	Vector4cd weak1_int(Triangle &T, Vector3d v, bool beta);
	Vector3cd weak3(Triangle &T, Vector3d v);
	Vector4cd weak3_int(Triangle &T, Vector3d v, bool q);
	Vector3cd weak4(Vector2d &p1, Vector2d &p2, Vector2d &p3, double v1, double v2, Vector2d &p6);
	Vector3cd weak4(Triangle &T, double v2, double v3, bool beta);
	double weak4_int(Vector2d p1, Vector2d p2, Vector2d b, double v1, double v2);
	void weak5(double x1, double y1, double x2, double y2, double v1, double v2, vector<vector<double>> &nbound, vector<complex<double>> &out5m, complex<double> &out5);
	VectorXcd setRightHand(TriangleMesh &U, TriangleMesh &L, vector<int> &nu);
	SparseMatrix<complex<double>> setA(double nn, vector<int> &nu, vector<vector<double>> &nbound, vector<vector<gridCell>> &grid);
	mxArray* setA_mx(double nn, vector<int> &nu, vector<vector<double>> &nbound, vector<vector<gridCell>> &grid) ;
	VectorXcd setB(vector<vector<gridCell>> &grid, VectorXcd &rh, int nn, vector<int> &nu);
	mxArray* setB_mx(vector<vector<gridCell>> &grid, VectorXcd &rh, int nn, vector<int> &nu);
	void assign(VectorXcd u, TriangleMesh &U, TriangleMesh &L, vector<int> &nu);
	VectorXcd getAperture(vector<vector<double>> &nbound, TriangleMesh &U, TriangleMesh &L);
	void drawAperture(VectorXd &plotX, VectorXd &plotY, vector<vector<double>> &nbound, int sign);

private:
	//setTri����Ҫʹ�õ��Ӻ���
	Triangle tri(double a1, double a2, double a3, double x1, double x2, double x3, double y1, double y2, double y3);
	void tri_int(double a1, double a2, double a3, double x1, double  x2, double  x3, double y1, double y2, double y3, Triangle &l);
	void tri_int1(double a1, double a2, double a3, double x1, double  x2, double  x3, double y1, double y2, double y3, Triangle &l);
	void solveTri(Triangle &l);
	void solveTri_int(Triangle &l);
	void solveTri_int1(Triangle &l);
	double findzero_int(double x1, double y1, double x2, double y2);
	void ustar_int(Triangle &T, Vector4d &u4, Vector4d &u5, Vector4d &u4_, Vector4d &u5_);
	void ustar2_int(Triangle &T, Vector4d &u5, Vector4d &u5_);
	void ustar_int1(Triangle &T, Vector4d &u4, Vector4d &u5, Vector4d &u4_, Vector4d &u5_);
	void ustar2_int1(Triangle &T, Vector4d &u5, Vector4d &u5_);

	//setGrid����Ҫʹ�õ��Ӻ���
	complex<double> value(double x, double y);
	void reloadtri(Triangle &T2, Triangle &T);
	void reloadtriex(Triangle &T3, Triangle &T);
	Vector4cd solveweak4(Triangle &T, double v2, double v3);
	int findXinNbound_multiple(vector<vector<double>> nbound, double x, int start, int end);

	//setRightHand����Ҫʹ�õ��Ӻ���
	complex<double> fInt_f(Triangle &T, Vector3d v);
	complex<double> fInt_int(Triangle &T, Vector3d v);
	complex<double> tr(Triangle &T, Vector3d v, bool f);

	//assign����Ҫʹ�õ��Ӻ���
	void assign_int(Triangle &T);
	void assign_intex(Triangle &T);
};