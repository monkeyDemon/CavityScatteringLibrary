#pragma once
#include "InhomogeneousThreeCircleCavity.h"

class _declspec(dllexport) InhomogeneousThreeCircleCircleCavity : public InhomogeneousThreeCircleCavity
{
public:
	InhomogeneousThreeCircleCircleCavity(unsigned int cavityType);

	//初始化腔体形状相关参数
	void InitCavityInhomogeneousShapeParameter(cavityCircle circle1, cavityCircle circle2, cavityCircle circle3);


protected:
	cavityCircle cavity1Circle_int;
	cavityCircle cavity2Circle_int;
	cavityCircle cavity3Circle_int;

	double phi_int(double x, double y);

};
