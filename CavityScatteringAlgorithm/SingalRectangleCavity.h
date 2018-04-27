#pragma once
#include "SingalCavity.h"

class _declspec(dllexport) SingalRectangleCavity : public SingalCavity
{
public:
	SingalRectangleCavity(unsigned int cavityType);

	//��ʼ��ǻ����״����
	void InitCavityShapeParameter(double cavityBottom);

protected:
	double cavityBottom;

	double phi(double x, double y);

	//�����Ӧ�ľ�����������о�ȷ�⣬��Ҫ��д������
	//complex<double> u(double x, double y);
};