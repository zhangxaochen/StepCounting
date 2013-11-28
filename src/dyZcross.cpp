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