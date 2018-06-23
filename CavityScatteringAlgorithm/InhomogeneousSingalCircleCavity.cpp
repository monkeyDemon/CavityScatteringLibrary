#include "InhomogeneousSingalCircleCavity.h"

InhomogeneousSingalCircleCavity::InhomogeneousSingalCircleCavity(unsigned int cavityType) :InhomogeneousSingalCavity(cavityType)
{
	int test = 0;
}

void InhomogeneousSingalCircleCavity::InitCavityShapeParameter(Vector2d circleCenter, double radius, Vector2d circleCenter_int, double radius_int)
{
	this->circleCenter = circleCenter;
	this->radius = radius;
	this->circleCenter_int = circleCenter_int;
	this->radius_int = radius_int;

	//Mark function InitCavityShapeParameter has been executed
	this->initalCheckKey = this->initalCheckKey | 16; // 16 means 0001 0000
}

double InhomogeneousSingalCircleCavity::phi(double x, double y)
{
	throw("暂未实现");
	double out = 0;
	return out;
}


double InhomogeneousSingalCircleCavity::phi_int(double x, double y)
{
	throw("暂未实现");
	double out = 0;
	return out;
}

