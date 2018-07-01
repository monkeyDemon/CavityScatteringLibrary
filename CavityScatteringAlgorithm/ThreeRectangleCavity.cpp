#include "ThreeRectangleCavity.h"

ThreeRectangleCavity::ThreeRectangleCavity(unsigned int cavityType) :InhomogeneousThreeRectangleCavity(cavityType)
{
	int test = 0;
}

void ThreeRectangleCavity::InitCavityInhomogeneousShapeParameter(cavityBox box1, cavityBox box2, cavityBox box3)
{
	this->cavity1Box_int = box1;
	this->cavity2Box_int = box2;
	this->cavity3Box_int = box3;

	//Mark function InitCavityShapeParameter has been executed
	this->initalCheckKey = this->initalCheckKey | 32; // 32 means 0010 0000
}


double ThreeRectangleCavity::phi_int(double x, double y)
{
	double x_left;
	double x_right;
	double y_bottom;
	double y_top;
	double out = 0;

	if (x < this->separator1_2)
	{
		x_left = this->cavity1Box_int.cavityLeft;
		x_right = this->cavity1Box_int.cavityRight;
		y_bottom = this->cavity1Box_int.cavityBottom;
		y_top = this->cavity1Box_int.cavityTop;
	}
	else if (separator1_2 <= x && x <= separator2_3)
	{
		x_left = this->cavity2Box_int.cavityLeft;
		x_right = this->cavity2Box_int.cavityRight;
		y_bottom = this->cavity2Box_int.cavityBottom;
		y_top = this->cavity2Box_int.cavityTop;
	}
	else
	{
		x_left = this->cavity3Box_int.cavityLeft;
		x_right = this->cavity3Box_int.cavityRight;
		y_bottom = this->cavity3Box_int.cavityBottom;
		y_top = this->cavity3Box_int.cavityTop;
	}


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