#pragma once
#include "SingalCavity.h"

class _declspec(dllexport) SingalRectangleCavity : public SingalCavity
{
public:
	SingalRectangleCavity(unsigned int cavityType);

	//初始化腔体形状参数
	void InitCavityShapeParameter(double cavityBottom);

protected:
	double cavityBottom;

	double phi(double x, double y);

	//子类对应的具体算例如果有精确解，需要覆写本方法
	//complex<double> u(double x, double y);
};