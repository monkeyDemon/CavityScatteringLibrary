#pragma once
#include "InhomogeneousSingalCavity.h"

class _declspec(dllexport) InhomogeneousSingalCircleCavity : public InhomogeneousSingalCavity
{
public:
	InhomogeneousSingalCircleCavity(unsigned int cavityType);

	//��ʼ��ǻ����״��ز���
	void InitCavityShapeParameter(Vector2d circleCenter, double radius, Vector2d circleCenter_int, double radius_int);


protected:
	Vector2d circleCenter;
	double radius;
	Vector2d circleCenter_int;
	double radius_int;


	double phi(double x, double y);
	double phi_int(double x, double y);

};
