#include "InhomogeneousThreeCircleCavity.h"

InhomogeneousThreeCircleCavity::InhomogeneousThreeCircleCavity(unsigned int cavityType) :InhomogeneousThreeCavity(cavityType)
{
	int test = 0;
}

void InhomogeneousThreeCircleCavity::InitCavityShapeParameter(cavityCircle circle1, cavityCircle circle2, cavityCircle circle3)
{
	this->cavity1Circle = circle1;
	this->cavity2Circle = circle2;
	this->cavity3Circle = circle3;


	//Mark function InitCavityShapeParameter has been executed
	this->initalCheckKey = this->initalCheckKey | 16; // 16 means 0001 0000
}

double InhomogeneousThreeCircleCavity::phi(double x, double y)
{
	double centerX, centerY;
	double radius;

	if (x < separator1_2)
	{
		// belongs to the first domain, the left cavity is in it
		centerX = cavity1Circle.circleCenter[0];
		centerY = cavity1Circle.circleCenter[1];
		radius = cavity1Circle.radius;
	}
	else if (separator1_2 <= x && x <= separator2_3)
	{
		// belongs to the second domain, the middle cavity is in it
		centerX = cavity2Circle.circleCenter[0];
		centerY = cavity2Circle.circleCenter[1];
		radius = cavity2Circle.radius;
	}
	else if (x > separator2_3)
	{
		//belongs to the third domain, the right cavity is in it
		centerX = cavity3Circle.circleCenter[0];
		centerY = cavity3Circle.circleCenter[1];
		radius = cavity3Circle.radius;
	}
	double out = pow(x - centerX, 2) + pow(y - centerY, 2) - pow(radius, 2);
	return out;
}