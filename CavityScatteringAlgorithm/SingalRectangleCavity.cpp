#include "SingalRectangleCavity.h"


SingalRectangleCavity::SingalRectangleCavity(unsigned int cavityType) :SingalCavity(cavityType)
{
	int test = 0;
}

void SingalRectangleCavity::InitCavityShapeParameter(double cavityBottom)
{
	//set the parameter of cavity

	//set the bottom of the cavity
	//here the cavity's shape is rectangle, the parameter 'cavityBottom' together with 'apertureLeft',
	//'apertureRight' and 'apertureY' can be used to describe the border of cavity
	this->cavityBottom = cavityBottom;

	//Mark function InitCavityShapeParameter has been executed
	this->initalCheckKey = this->initalCheckKey | 16; // 16 means 0001 0000
}



double SingalRectangleCavity::phi(double x, double y)
{
	double x_left = this->apertureLeft;
	double x_right = this->apertureRight;
	double y_bottom = this->cavityBottom;
	double y_top = this->apertureY;
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