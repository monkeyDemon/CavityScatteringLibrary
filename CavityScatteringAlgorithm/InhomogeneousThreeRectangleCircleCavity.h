#pragma once
#include "InhomogeneousThreeRectangleCavity.h"

class _declspec(dllexport) InhomogeneousThreeRectangleCircleCavity : public InhomogeneousThreeRectangleCavity
{
public:
	InhomogeneousThreeRectangleCircleCavity(unsigned int cavityType);

	//��ʼ��ǻ����״��ز���
	void InitCavityInhomogeneousShapeParameter(cavityCircle circle1, cavityCircle circle2, cavityCircle circle3);


protected:
	cavityCircle cavity1Circle;
	cavityCircle cavity2Circle;
	cavityCircle cavity3Circle;

	double phi_int(double x, double y);

};
