#ifndef SCMETHOD_H
#define SCMETHOD_H

#include <cstddef>
#include <vector>
#include <math.h>

using std::vector;

const double epsilon=1e-9;
const int NAG_INF=-100000;
const int INVALID=-100000;

// (th, freq) 法的 alpha， beta
const float a=0.0f;
const float b=1.0f;
// 计算小段内方差，方差太小，认为停止
const float varTh=0.5;
// 是否做低通滤波，做
const bool doLpf=true;
// 滤波窗口长度，用不到，实际取 FPS/3
const int numtaps=10;
// 统一步数补偿
const int compensation=01;
//达不到scTh停下，则认为不是真正行走，计步归零
const int scTh=5;

const size_t BUFSZ=10;
const int S2NS=1000;

const int METHOD_DYTH=0,
	METHOD_DYPEAK=1,
	METHOD_DYZCROSS=2;

#pragma region (th, freq) 法
template<typename T>
size_t dyThresholdOffline(vector<T> dataSet);

size_t dyThresholdOnline(double value, double ts, bool doLpf=false);

size_t getDythSteps();
void resetDythCounter();

#pragma endregion //--------------(th, freq) 法


#pragma region dyPeak 法

size_t dyPeakOnline(double value, double ts, bool doLpf=false);
size_t getDypeakSteps();
void resetDypeakCounter();

template<typename T>
double getPeakMean(const T &buf, size_t begin, size_t end){
	double pmean=NAG_INF;

	int pc=0;
	double psum=0;
	for(size_t i=begin; i<end; i++){
		if(i==0 || i==end-1)
			continue;
		double v=buf[i];
		double forwardSlope=buf[i+1]-v;
		double backwardSlope=v-buf[i-1];
		if(forwardSlope<0 && backwardSlope>0){
			pc++;
			psum+=v;
		}
	}

	if(pc!=0)
		pmean=psum*1.0/pc;
	return pmean;
}//getPeakMean


#pragma endregion //--------------dyPeak 法


#pragma region dyZcross 法
size_t dyZcrossOnline(double value, double ts, bool doLpf=false);
size_t getDyzcSteps();
void resetDyzcCounter();




#pragma endregion //--------------dyZcross 法


inline size_t pushData(double value, double ts, int whichMethod=METHOD_DYTH){
	size_t steps=NAG_INF;
	switch (whichMethod)
	{
	case METHOD_DYTH:
		steps=dyThresholdOnline(value, ts, true);
		break;
	case METHOD_DYPEAK:
		steps=dyPeakOnline(value, ts, true);
		break;
	case METHOD_DYZCROSS:
		steps=dyZcrossOnline(value, ts, true);
		break;
	default:
		break;
	}
	return steps;
}//pushData

inline size_t pushData(double ax, double ay, double az, double ts, int whichMethod=METHOD_DYTH){
	double totAcc=sqrt(ax*ax+ay*ay+az*az);
	return pushData(totAcc, ts, whichMethod);
}//pushData

inline size_t getSteps(int whichMethod=METHOD_DYTH){
	size_t steps;
	switch (whichMethod)
	{
	case METHOD_DYTH:
		steps=getDythSteps();
		break;
	case METHOD_DYPEAK:
		steps=getDypeakSteps();
		break;
	case METHOD_DYZCROSS:
		steps=getDyzcSteps();
		break;
	default:
		break;
	}
	return steps;
}//getSteps

inline void resetCounter(int whichMethod=METHOD_DYTH){
	switch (whichMethod)
	{
	case METHOD_DYTH:
		resetDythCounter();
		break;
	case METHOD_DYPEAK:
		resetDypeakCounter();
		break;
	case METHOD_DYZCROSS:
		resetDyzcCounter();
		break;
	default:
		break;
	}
}//resetCounter

#endif//SCMETHOD_H
