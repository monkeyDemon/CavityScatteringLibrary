#pragma once
#include <time.h>
#include <iostream>

class MyTimer
{
public:
	MyTimer(int timerType)
	{
		this->timerType = timerType;
		//timerType �����趨��ʱ������𣬣���δʵ�֣�Ŀǰֻ��һ�����
		//��ͬ�ļ�ʱ��ʵ������Ӧ��ƽ̨�����Ȳ�ͬ
		//Ŀǰʵ�ֵļ�ʱ���ܹ���ƽ̨ʹ�ã����Ƚϵ͵�����
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