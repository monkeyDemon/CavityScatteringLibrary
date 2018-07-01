#pragma once
#include "InhomogeneousThreeRectangleCavity.h"

/*
由于普通均匀介质的多腔体可以看做是
非均匀介质的多腔体的特例（即具有两种相同的介质的多腔体）
因此将普通均匀介质多腔体的算例代码内嵌到非均匀介质的基类下
这里ThreeRectangleCavity的实现完全照搬InhomogeneousThreeRectangleRectangleCavity
调用时，只需人为确保两种介质相同（赋值相同的k和epr）即可

（当然最优逻辑是ThreeRectangleCavity类继承Cavity类，保持结构的一致性，静待有缘人来完善吧）
*/
class _declspec(dllexport) ThreeRectangleCavity : public InhomogeneousThreeRectangleCavity
{
public:
	ThreeRectangleCavity(unsigned int cavityType);

	//初始化腔体形状相关参数
	void InitCavityInhomogeneousShapeParameter(cavityBox box1, cavityBox box2, cavityBox box3);


protected:
	cavityBox cavity1Box_int;
	cavityBox cavity2Box_int;
	cavityBox cavity3Box_int;

	double phi_int(double x, double y);

};
