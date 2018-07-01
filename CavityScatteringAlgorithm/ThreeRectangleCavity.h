#pragma once
#include "InhomogeneousThreeRectangleCavity.h"

/*
������ͨ���Ƚ��ʵĶ�ǻ����Կ�����
�Ǿ��Ƚ��ʵĶ�ǻ���������������������ͬ�Ľ��ʵĶ�ǻ�壩
��˽���ͨ���Ƚ��ʶ�ǻ�������������Ƕ���Ǿ��Ƚ��ʵĻ�����
����ThreeRectangleCavity��ʵ����ȫ�հ�InhomogeneousThreeRectangleRectangleCavity
����ʱ��ֻ����Ϊȷ�����ֽ�����ͬ����ֵ��ͬ��k��epr������

����Ȼ�����߼���ThreeRectangleCavity��̳�Cavity�࣬���ֽṹ��һ���ԣ�������Ե�������ưɣ�
*/
class _declspec(dllexport) ThreeRectangleCavity : public InhomogeneousThreeRectangleCavity
{
public:
	ThreeRectangleCavity(unsigned int cavityType);

	//��ʼ��ǻ����״��ز���
	void InitCavityInhomogeneousShapeParameter(cavityBox box1, cavityBox box2, cavityBox box3);


protected:
	cavityBox cavity1Box_int;
	cavityBox cavity2Box_int;
	cavityBox cavity3Box_int;

	double phi_int(double x, double y);

};
