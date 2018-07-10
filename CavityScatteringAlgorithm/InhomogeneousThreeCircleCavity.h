#pragma once
#include "InhomogeneousThreeCavity.h"

class _declspec(dllexport) InhomogeneousThreeCircleCavity : public InhomogeneousThreeCavity
{
public:
	InhomogeneousThreeCircleCavity(unsigned int cavityType);

	//��ʼ��ǻ����״��ز���
	void InitCavityShapeParameter(cavityCircle circle1, cavityCircle circle2, cavityCircle circle3);

protected:
	cavityCircle cavity1Circle;
	cavityCircle cavity2Circle;
	cavityCircle cavity3Circle;

	double phi(double x, double y);

};
