#pragma once
#include <time.h>
#include <iostream>

class MyTimer
{
public:
	MyTimer(int timerType)
	{
		this->timerType = timerType;
		//timerType 用于设定计时器的类别，（暂未实现，目前只有一种类别）
		//不同的计时器实现所对应的平台及精度不同
		//目前实现的计时器能够跨平台使用，精度较低但够用
	}
	void Start(char* operation)
	{
		this->operation = operation;
		this->start = clock();
	}
	void End()
	{
		this->end = clock();
	}
	void EndAndPrint()
	{
		this->end = clock();
		double dur = (double)(end - start);
		printf("%s use Time:%f\n", operation,(dur / CLOCKS_PER_SEC));
	}
private:
	int timerType;
	char * operation;
	clock_t start, end;
};