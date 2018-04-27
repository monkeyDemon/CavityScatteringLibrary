#pragma once
#include "InhomogeneousSingalCavity.h"

class _declspec(dllexport) InhomogeneousSingalRectangleCavity : public InhomogeneousSingalCavity
{
public:
	InhomogeneousSingalRectangleCavity(unsigned int cavityType);

	//初始化腔体形状相关参数
	void InitCavityShapeParameter(double cT, double cB, double cL,double cR, double cT_, double cB_, double cL_, double cR_);


protected:
	double cavityLeft;
	double cavityRight;
	double cavityTop;
	double cavityBottom;

	double cavityLeft_int;
	double cavityRight_int;
	double cavityTop_int;
	double cavityBottom_int;

	double phi(double x, double y);
	double phi_int(double x, double y);

};
