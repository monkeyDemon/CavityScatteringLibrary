#pragma once
#include "InhomogeneousThreeCavity.h"

class _declspec(dllexport) InhomogeneousThreeRectangleCavity : public InhomogeneousThreeCavity
{
public:
	InhomogeneousThreeRectangleCavity(unsigned int cavityType);

	//初始化腔体形状相关参数
	void InitCavityShapeParameter(cavityBox box1, cavityBox box2,cavityBox box3);


protected:
	cavityBox cavity1Box;
	cavityBox cavity2Box;
	cavityBox cavity3Box;

	double phi(double x, double y);
	//double phi_int(double x, double y);

};