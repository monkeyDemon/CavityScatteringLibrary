#pragma once
#include "InhomogeneousThreeCircleCavity.h"

class _declspec(dllexport) InhomogeneousThreeCircleRectangleCavity : public InhomogeneousThreeCircleCavity
{
public:
	InhomogeneousThreeCircleRectangleCavity(unsigned int cavityType);

	//初始化腔体形状相关参数
	void InitCavityInhomogeneousShapeParameter(cavityBox box1, cavityBox box2, cavityBox box3);


protected:
	cavityBox cavity1Box;
	cavityBox cavity2Box;
	cavityBox cavity3Box;

	double phi_int(double x, double y);

};