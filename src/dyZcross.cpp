#include "utils.h"
#include "scMethods.h"
#include "circular.h"
#include <vector>
using namespace std;
typedef circular_buffer<double> cbuf_type;

static bool start=false;

static size_t steps=0;
static size_t i=0;
static size_t fps=NAG_INF;
//相当于offline的numtaps
static size_t winsz=NAG_INF;

static const int CAPACITY=250;

static cbuf_type cbuf_ts(CAPACITY);
//之后需要 reserve 到 2*winsz：
static cbuf_type cbuf_v(CAPACITY);
static cbuf_type cbuf_v_lpf(CAPACITY);
static const cbuf_type *cbufPtr=NULL;

static vector<double> *hamWin;
static const int overlap=1;
static size_t lastZcrossIdx=0;

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


size_t dyZcrossOnline(double value, double ts, bool doLpf){
	/*cout<<"i: "<<i<<endl;*/
	/*BUFSZ 用来算 fps, winsz 用来算 hamWin, fps 用作此方法实际缓存大小 ×→ 仍用 winsz*/

	if(i<BUFSZ)
		cbuf_ts.push_back(ts);
	else if(fps==INVALID){
		fps=int(S2NS*(BUFSZ-1)*1.0/(cbuf_ts[BUFSZ-1]-cbuf_ts[0]));
		winsz=int(fps/3);
		hamWin=new vector<double>(getHammingWin(winsz) );
	}

	cbuf_v.push_back(value);
	if(fps==INVALID)
		cbuf_v_lpf.push_back(value);
	else
		cbuf_v_lpf.push_back(getDataLpf(cbuf_v, *hamWin));

	//real 真正要处理的缓存大小
	size_t rbufsz=fps+0;
	//rbufsz=winsz;

	if(fps==INVALID || cbuf_v.size()<rbufsz*2){
		i++;
		return steps;
	}
	/*#若 i>=rbufsz*2:*/
	if(doLpf)
		cbufPtr=&cbuf_v_lpf;
	else
		cbufPtr=&cbuf_v;
	const cbuf_type &buf=*cbufPtr;
	
	double va=max(maxMinVar(buf, 0, rbufsz), maxMinVar(buf, rbufsz-1, rbufsz*2-1) );
	double zline=getZeroLine(buf, 0, rbufsz);
	if(va<varTh){
		for(size_t i=0; i<winsz; i++){
			cbuf_v.pop_front();
			cbuf_v_lpf.pop_front();
		}
	}
	else{
		/*cout<<"=============steps+="<<endl;*/
		steps+=getStepByZline(buf, 0, rbufsz, zline, i+1-rbufsz*2);
		for(size_t i=0; i<rbufsz-overlap; i++){
			cbuf_v.pop_front();
			cbuf_v_lpf.pop_front();
		}
	}
	i++;
	return steps;
}//dyZcrossOnline


void resetDyzcCounter(){
	lastZcrossIdx=0;
	start=false;

	i=0;
	steps=0;
	fps=NAG_INF;
	winsz=NAG_INF;

	cbuf_ts.clear();
	cbuf_v.clear();
	cbuf_v_lpf.clear();

}//resetDyzcCounter

size_t getDyzcSteps(){
	return steps;
}//getDyzcSteps