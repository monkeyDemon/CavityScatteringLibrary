#pragma once
#include "InhomogeneousThreeRectangleCavity.h"

class _declspec(dllexport) InhomogeneousThreeRectangleRectangleCavity : public InhomogeneousThreeRectangleCavity
{
public:
	InhomogeneousThreeRectangleRectangleCavity(unsigned int cavityType);

	//初始化腔体形状相关参数
	void InitCavityInhomogeneousShapeParameter(cavityBox box1, cavityBox box2, cavityBox box3);


protected:
	cavityBox cavity1Box_int;
	cavityBox cavity2Box_int;
	cavityBox cavity3Box_int;

	double phi_int(double x, double y);

};