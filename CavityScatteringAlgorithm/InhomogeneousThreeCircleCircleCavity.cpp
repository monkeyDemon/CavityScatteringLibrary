#include "InhomogeneousThreeCircleCircleCavity.h"

InhomogeneousThreeCircleCircleCavity::InhomogeneousThreeCircleCircleCavity(unsigned int cavityType) :InhomogeneousThreeCircleCavity(cavityType)
{
	int test = 0;
}

void InhomogeneousThreeCircleCircleCavity::InitCavityInhomogeneousShapeParameter(cavityCircle circle1, cavityCircle circle2, cavityCircle circle3)
{
	this->cavity1Circle_int = circle1;
	this->cavity2Circle_int = circle2;
	this->cavity3Circle_int = circle3;


	//Mark function InitCavityShapeParameter has been executed
	this->initalCheckKey = this->initalCheckKey | 32; // 32 means 0010 0000
}


double InhomogeneousThreeCircleCircleCavity::phi_int(double x, double y)
{
	double centerX, centerY;
	double radius;

	if (x < separator1_2)
	{
		// belongs to the first domain, the left cavity is in it
		centerX = cavity1Circle_int.circleCenter[0];
		centerY = cavity1Circle_int.circleCenter[1];
		radius = cavity1Circle_int.radius;
	}
	else if (separator1_2 <= x && x <= separator2_3)
	{
		// belongs to the second domain, the middle cavity is in it
		centerX = cavity2Circle_int.circleCenter[0];
		centerY = cavity2Circle_int.circleCenter[1];
		radius = cavity2Circle_int.radius;
	}
	else if (x > separator2_3)
	{
		//belongs to the third domain, the right cavity is in it
		centerX = cavity3Circle_int.circleCenter[0];
		centerY = cavity3Circle_int.circleCenter[1];
		radius = cavity3Circle_int.radius;
	}
	double out = pow(x - centerX, 2) + pow(y - centerY, 2) - pow(radius, 2);
	return out;
}