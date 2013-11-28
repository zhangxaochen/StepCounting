#ifndef SCMETHOD_H
#define SCMETHOD_H

#include <vector>
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

template<typename T>
size_t getStepByPmean(const T &buf, size_t begin, size_t end, double pmean, size_t startp){
	size_t sc=0;
	for(size_t i=begin; i<end; i++){
		if(i==0 || i==end-1)
			continue;
		double v=buf[i];
		double forwardSlope=buf[i+1]-v;
		double backwardSlope=v-buf[i-1];
		if(forwardSlope<0 && backwardSlope>0 &&
			v>coeff*pmean && v>baseTh)
		{
			if(startp+i-lastZcrossIdx>fps/4){
				//cout<<"startp: "<<startp<<"\t"<<startp+i<<endl;
				sc++;
				lastZcrossIdx=startp+i;
			}
		}
	}
	return sc;
}//getStepByPmean

#pragma endregion //--------------dyPeak 法


#pragma region dyZcross 法
size_t dyZcrossOnline(double value, double ts, bool doLpf=false);
size_t getDyzcSteps();
void resetDyzcCounter();

template<typename T>
double getZeroLine(const T &buf, size_t begin, size_t end){
	int pc=0;
	int vc=0;
	double lastPeak=NAG_INF;
	double lastValley=NAG_INF;

	double psum=0;
	double vsum=0;

	for(size_t i=begin; i<end; i++){
		if(i==0 || i==end-1)
			continue;
		double v=buf[i];
		double forwardSlope=buf[i+1]-v;
		double backwardSlope=v-buf[i-1];

		//一段中 pc, vc 至多算2个，防止静止时期的毛刺干扰（尽管lpf，毛刺仍存在）
		//abs(lastValley-v)>0.1 用于过滤掉中间没滤掉的褶皱， 下同。。
		//peak:
		if(forwardSlope<0 && backwardSlope>0 && pc<2){
			if(abs(lastValley-v)<0.1){
				vsum-=lastValley;
				vc--;
				lastValley=NAG_INF;
			}
			else
			{
				lastPeak=v;
				pc++;
				psum+=v;
			}
		}
		//valley:
		if(forwardSlope>0 && backwardSlope<0 && vc<2){
			if(abs(lastPeak-v)<0.1){
				psum-=lastPeak;
				pc--;
				lastPeak=NAG_INF;
			}
			else
			{
				lastValley=v;
				vc++;
				vsum+=v;
			}
		}
	}//for
	//#若波峰或波谷有一个不存在，以端点值为波峰or波谷
	if(pc==0){
		pc=2;
		psum=buf[begin]+buf[end-1];
	}
	//#其实波峰or波谷只可能一个不存在：
	if(vc==0){
		vc=2;
		vsum=buf[begin]+buf[end-1];
	}

	double pvmean=(psum*1.0/pc+vsum*1.0/vc)/2;
	return pvmean;
}//getZeroLine

/**#统计过零点次数
 **@param startp 本段buf起点在全局中的位置
 **
 **/
template<typename T>
size_t getStepByZline(const T &buf, size_t begin, size_t end, double zline, size_t startp){
	//cout<<"in getStepByZline"<<endl;
	size_t sc=0;
	for(size_t i=begin; i<end; i++){
		if(i==0)
			continue;
		double v=buf[i];
		double vprev=buf[i-1];
		//if((v-zline)*(vprev-zline)<0){
		//只看下降沿：
		if(vprev-zline>0 && v-zline<0){
			//cout<<"down,,,,,"<<vprev<<"\t"<<zline<<"\t"<<v<<endl;
			if(startp+i-lastZcrossIdx>2){
				//cout<<"startp: "<<startp<<"\t"<<startp+i<<endl;
				sc++;
				lastZcrossIdx=startp+i;
			}
		}
	}
	return sc;
}//getStepByZline


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
