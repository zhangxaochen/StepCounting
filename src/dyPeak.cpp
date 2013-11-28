#include "utils.h"
#include "scMethods.h"
#include "circular.h"
#include <vector>
using namespace std;
typedef circular_buffer<double> cbuf_type;

static const int overlap=1;
static const double coeff=0.8;
static const double baseTh=2.;
static const int CAPACITY=100;

static bool start=false;
static size_t steps=0;
static size_t i=0;
static size_t fps=NAG_INF;
//相当于offline的numtaps
static size_t winsz=NAG_INF;
static size_t lastZcrossIdx=0;



static cbuf_type cbuf_ts(CAPACITY);
//之后需要 reserve 到 2*winsz：
static cbuf_type cbuf_v(CAPACITY);
static cbuf_type cbuf_v_lpf(CAPACITY);
static const cbuf_type *cbufPtr=NULL;
static vector<double> *hamWin;


size_t dyPeakOnline(double value, double ts, bool doLpf){
	/*cout<<"i: "<<i<<endl;*/
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
	if(va<varTh){
		//del buf[:]
		for(size_t i=0; i<rbufsz; i++){
			/*cout<<"del buf[:], "<<i<<endl;*/
			cbuf_v.pop_front();
			cbuf_v_lpf.pop_front();
		}
	}
	//即 va>=varTh 时：
	else{
		double pmean=getPeakMean(buf, 0, rbufsz);
		/*cout<<"~~~~~~"<<i<<endl;*/
		steps+=getStepByPmean(buf, 0, rbufsz, pmean, i+1-rbufsz*2);

		for(size_t i=0; i<rbufsz-overlap; i++){
			cbuf_v.pop_front();
			cbuf_v_lpf.pop_front();
		}
	}
	return steps;
}//dyPeakOnline

void resetDypeakCounter(){
	lastZcrossIdx=0;
	start=false;

	i=0;
	steps=0;
	fps=NAG_INF;
	winsz=NAG_INF;

	cbuf_ts.clear();
	cbuf_v.clear();
	cbuf_v_lpf.clear();
}//resetDypeakCounter

size_t getDypeakSteps(){
	return steps;
}//getDypeakSteps