#include "InhomogeneousSingalRectangleCavity.h"

InhomogeneousSingalRectangleCavity::InhomogeneousSingalRectangleCavity(unsigned int cavityType) :InhomogeneousSingalCavity(cavityType)
{
	int test = 0;
}

void InhomogeneousSingalRectangleCavity::InitCavityShapeParameter(double cT, double cB, double cL, double cR, double cT_, double cB_, double cL_, double cR_)
{
	this->cavityTop = cT;
	this->cavityBottom = cB;
	this->cavityLeft = cL;
	this->cavityRight = cR;

	this->cavityTop_int = cT_;
	this->cavityBottom_int = cB_;
	this->cavityLeft_int = cL_;
	this->cavityRight_int = cR_;

	//Mark function InitCavityShapeParameter has been executed
	this->initalCheckKey = this->initalCheckKey | 16; // 16 means 0001 0000
}

double InhomogeneousSingalRectangleCavity::phi(double x, double y)
{
	double x_left = this->cavityLeft;
	double x_right = this->cavityRight;
	double y_bottom = this->cavityBottom;
	double y_top = this->cavityTop;
	double out = 0;

	double dist_x, dist_y;
	if (x > x_right || x < x_left || y<y_bottom || y>y_top)
	{
		dist_x = 0;
		if (x > x_right)
			dist_x = x - x_right;
		if (x < x_left)
			dist_x = x_left - x;

		dist_y = 0;
		if (y > y_top)
			dist_y = y - y_top;
		if (y < y_bottom)
			dist_y = y_bottom - y;

		out = sqrt(dist_x*dist_x + dist_y*dist_y);
	}
	else if (x<x_right && x>x_left && y > y_bottom && y < y_top)
	{
		dist_x = x_right - x;
		dist_x = (x - x_left) < dist_x ? (x - x_left) : dist_x;
		dist_y = y_top - y;
		dist_y = (y - y_bottom) < dist_y ? (y - y_bottom) : dist_y;
		out = -1 * (dist_x < dist_y ? dist_x : dist_y);
	}
	else if (y == y_top && x<x_right && x>x_left)
	{
		dist_x = x_right - x;
		dist_x = (x - x_left) < dist_x ? (x - x_left) : dist_x;
		out = -1 * dist_x;
	}
	else
		out = 0;
	return out;
}


double InhomogeneousSingalRectangleCavity::phi_int(double x, double y)
{
	double x_left = this->cavityLeft_int;
	double x_right = this->cavityRight_int;
	double y_bottom = this->cavityBottom_int;
	double y_top = this->cavityTop_int;
	double out = 0;

	double dist_x, dist_y;
	if (x > x_right || x < x_left || y<y_bottom || y>y_top)
	{
		dist_x = 0;
		if (x > x_right)
			dist_x = x - x_right;
		if (x < x_left)
			dist_x = x_left - x;

		dist_y = 0;
		if (y > y_top)
			dist_y = y - y_top;
		if (y < y_bottom)
			dist_y = y_bottom - y;

		out = sqrt(dist_x*dist_x + dist_y*dist_y);
	}
	else if (x<x_right && x>x_left && y > y_bottom && y < y_top)
	{
		dist_x = x_right - x;
		dist_x = (x - x_left) < dist_x ? (x - x_left) : dist_x;
		dist_y = y_top - y;
		dist_y = (y - y_bottom) < dist_y ? (y - y_bottom) : dist_y;
		out = -1 * (dist_x < dist_y ? dist_x : dist_y);
	}
	else
		out = 0;
	return out;
}
