#pragma once
#include<complex>
#include<math.h>
#include<vector>
#include<Eigen\Dense>

using namespace std;
using namespace Eigen;

struct cavityBox
{
	double cavityLeft;
	double cavityRight;
	double cavityTop;
	double cavityBottom;
};


struct cavityCircle
{
	Vector2d circleCenter;  //圆心
	double radius;               //半径
};



struct gridCell
{
	complex<double> _c;
	complex<double> _w;
	complex<double> _s;
	complex<double> _n;
	complex<double> _e;
	complex<double> _nw;
	complex<double> _se;
	complex<double> _const;
	vector<complex<double>> tm;
};

struct Triangle
{
	Vector3d x;
	Vector3d y;
	int sgn;
	double r1;
	double r2;
	Vector2d A;
	Vector2d B;
	Vector2d C;
	Vector3i v;
	double x1;
	double y1;
	double x2;
	double y2;
	double x3;
	double y3;
	double x4;
	double y4;
	double x5;
	double y5;
	Vector2d p5;
	Vector3cd z;
	complex<double> z1;
	complex<double> z2;
	complex<double> z3;
	complex<double> z4;
	complex<double> z5;

	//非均匀介质时需要处理的属性
	Vector3d x_int;
	Vector3d y_int;
	int sgn_int;
	double r1_int;
	double r2_int;
	Vector2d A_int;
	Vector2d B_int;
	Vector2d C_int;
	Vector3i v_int;
	double x1_int;
	double y1_int;
	double x2_int;
	double y2_int;
	double x3_int;
	double y3_int;
	double x4_int;
	double y4_int;
	double x5_int;
	double y5_int;
	Vector2d p5_int;
	Vector4d u4_int;
	Vector4d u4_int_;
	Vector4d u5_int;
	Vector4d u5_int_;
	Vector3cd z_int;
	complex<double> z1_int;
	complex<double> z2_int;
	complex<double> z2_int_;
	complex<double> z3_int;
	complex<double> z4_int;
	complex<double> z5_int;
	complex<double> z4_int_;
	complex<double> z5_int_;

	//
	Vector3d x_intex;
	Vector3d y_intex;
	int sgn_intex;
	double r1_intex;
	double r2_intex;
	Vector2d A_intex;
	Vector2d B_intex;
	Vector2d C_intex;
	Vector3i v_intex;
	double x1_intex;
	double y1_intex;
	double x2_intex;
	double y2_intex;
	double x3_intex;
	double y3_intex;
	double x4_intex;
	double y4_intex;
	double x5_intex;
	double y5_intex;
	Vector2d p5_intex;
	Vector4d u4_intex;
	Vector4d u4_intex_;
	Vector4d u5_intex;
	Vector4d u5_intex_;
	Vector3cd z_intex;
	complex<double> z1_intex;
	complex<double> z2_intex;
	complex<double> z2_intex_;
	complex<double> z3_intex;
	complex<double> z4_intex;
	complex<double> z5_intex;
	complex<double> z4_intex_;
	complex<double> z5_intex_;

	int n;

	Triangle()
	{
		x = Vector3d(3);
		y = Vector3d(3);
		sgn = INT_MIN;
		r1 = DBL_MIN;
		r2 = DBL_MIN;
		A = Vector2d(2);
		B = Vector2d(2);
		C = Vector2d(2);
		v = Vector3i(3);
		x1 = DBL_MIN;
		y1 = DBL_MIN;
		x2 = DBL_MIN;
		y2 = DBL_MIN;
		x3 = DBL_MIN;
		y3 = DBL_MIN;
		x4 = DBL_MIN;
		y4 = DBL_MIN;
		x5 = DBL_MIN;
		y5 = DBL_MIN;
		p5 = Vector2d(2);
		z = Vector3cd(3);
		z1 = complex<double>(DBL_MIN, DBL_MIN);
		z2 = complex<double>(DBL_MIN, DBL_MIN);
		z3 = complex<double>(DBL_MIN, DBL_MIN);
		z4 = complex<double>(DBL_MIN, DBL_MIN);
		z5 = complex<double>(DBL_MIN, DBL_MIN);

		//非均匀介质时需要处理的属性
		x_int = Vector3d(3);
		y_int = Vector3d(3);
		sgn_int = INT_MIN;
		r1_int = DBL_MIN;
		r2_int = DBL_MIN;
		A_int = Vector2d(2);
		B_int = Vector2d(2);
		C_int = Vector2d(2);
		v_int = Vector3i(3);
		x1_int = DBL_MIN;
		y1_int = DBL_MIN;
		x2_int = DBL_MIN;
		y2_int = DBL_MIN;
		x3_int = DBL_MIN;
		y3_int = DBL_MIN;
		x4_int = DBL_MIN;
		y4_int = DBL_MIN;
		x5_int = DBL_MIN;
		y5_int = DBL_MIN;
		p5_int = Vector2d(2);
		u4_int = Vector4d(4);
		u4_int_ = Vector4d(4);
		u5_int = Vector4d(4);
		u5_int_ = Vector4d(4);
		z_int = Vector3cd(3);
		z1_int = complex<double>(DBL_MIN, DBL_MIN);
		z2_int = complex<double>(DBL_MIN, DBL_MIN);
		z2_int_ = complex<double>(DBL_MIN, DBL_MIN);
		z3_int = complex<double>(DBL_MIN, DBL_MIN);
		z4_int = complex<double>(DBL_MIN, DBL_MIN);
		z5_int = complex<double>(DBL_MIN, DBL_MIN);
		z4_int_ = complex<double>(DBL_MIN, DBL_MIN);
		z5_int_ = complex<double>(DBL_MIN, DBL_MIN);

		x_intex = Vector3d(3);
		y_intex = Vector3d(3);
		sgn_intex = INT_MIN;
		r1_intex = DBL_MIN;
		r2_intex = DBL_MIN;
		A_intex = Vector2d(2);
		B_intex = Vector2d(2);
		C_intex = Vector2d(2);
		v_intex = Vector3i(3);
		x1_intex = DBL_MIN;
		y1_intex = DBL_MIN;
		x2_intex = DBL_MIN;
		y2_intex = DBL_MIN;
		x3_intex = DBL_MIN;
		y3_intex = DBL_MIN;
		x4_intex = DBL_MIN;
		y4_intex = DBL_MIN;
		x5_intex = DBL_MIN;
		y5_intex = DBL_MIN;
		p5_intex = Vector2d(2);
		u4_intex = Vector4d(4);
		u4_intex_ = Vector4d(4);
		u5_intex = Vector4d(4);
		u5_intex_ = Vector4d(4);
		z_intex = Vector3cd(3);
		z1_intex = complex<double>(DBL_MIN, DBL_MIN);
		z2_intex = complex<double>(DBL_MIN, DBL_MIN);
		z2_intex_ = complex<double>(DBL_MIN, DBL_MIN);
		z3_intex = complex<double>(DBL_MIN, DBL_MIN);
		z4_intex = complex<double>(DBL_MIN, DBL_MIN);
		z5_intex = complex<double>(DBL_MIN, DBL_MIN);
		z4_intex_ = complex<double>(DBL_MIN, DBL_MIN);
		z5_intex_ = complex<double>(DBL_MIN, DBL_MIN);

		n = INT_MIN;
	}
};


class TriangleMesh
{
public:
	TriangleMesh(int meshWidth, int meshHeight)
	{
		tri.resize(meshWidth);
		for (int i = 0; i < tri.size(); i++)
			tri[i].resize(meshHeight);
	}

	Triangle Get(int i, int j)
	{
		return tri[i][j];
	}
	Triangle* GetPt(int i, int j)
	{
		return &tri[i][j];
	}
	void Set(int i, int j, Triangle &curTri)
	{
		tri[i][j] = curTri;
	}
private:
	vector<vector<Triangle>> tri;
};