#include "utils.h"
#include "scMethods.h"
#include "circular.h"
#include <vector>
using namespace std;
typedef circular_buffer<double> cbuf_type;

// ������ֵ
const double baseTh=1;

static size_t k=0;
static double maxv=NAG_INF;
static size_t i=0;
static size_t steps=0;
static bool start=false;

static double fps=NAG_INF;
//�൱��offline��numtaps
static size_t winsz=NAG_INF;

// vector<double> buf;
static cbuf_type cbuf_ts(BUFSZ);
//֮����Ҫ reserve �� 2*winsz��
static cbuf_type cbuf_v(BUFSZ);
static cbuf_type cbuf_v_lpf(BUFSZ);
static const cbuf_type *cbufPtr=NULL;

static vector<double> *hamWin;


template size_t dyThresholdOffline<double>(vector<double> dataSet);
template<typename T>
size_t dyThresholdOffline(vector<T> dataSet){
	size_t k=0;
	double max=dataSet[k];
	size_t i=k+1;
	// 	�����������εļƲ�������Ϊ���� scTh �Ա�
	size_t steps=0;
	// 	�������ݵļƲ�������
	size_t totSteps=0;
	//�費��Ҫstart��ʶ��
	bool start=false;

	for(i=0; i<dataSet.size(); i++){
		double v=dataSet[i];

		//i ��1��ʼ���� dataSet.size()-numtaps ��ֹ
		if(i==0 || i==dataSet.size()-1 
			|| i+numtaps>=dataSet.size())
			continue;
		// 		if (i==999)
		// 			cout<<i<<" maxMinVar(dataSet, i, i+numtaps): "<<maxMinVar(dataSet, i, i+numtaps)<<endl;

		if (maxMinVar(dataSet, i, i+numtaps)<varTh){
			if(start){
				start=false;
				// 				if(steps>=scTh)
				// 					totSteps+=steps;
				// 				steps=0;
			}
			continue;
		}
		else
		{
			if(!start)
				start=true;
		}

		double th=a*1.f/(i-k)+b;
		if (v<dataSet[i+1] && max-v>=th && max-v>baseTh){
			steps++;
			k=i;
			max=v;
			// 			cout<<"i: "<<i<<endl;
		}
		else{
			if(v>max){
				max=v;
				k=i;
			}
		}
	}
	return steps+compensation;
}//dyThresholdOffline

size_t dyThresholdOnline(double value, double ts, bool doLpf){
	//cout<<"========i, value, ts: "<<i<<", "<<value<<", "<<ts<<endl;
	//���bufsz���˲����ڳ��ȣ� ֻ�����ڼ��� fps
	if(cbuf_ts.size()<BUFSZ){
		/*cout<<"cbuf_ts.size()<BUFSZ"<<endl;*/
		//��ֹ��һ֡���⣺
		if(cbuf_ts.size()==1 && abs(cbuf_ts.front()-ts)>S2NS){
			cout<<"@@@@@@@cbuf_ts.pop_front()"<<endl;
			cbuf_v.pop_front();
			cbuf_v_lpf.pop_front();

			cbuf_ts.pop_front();
		}
		cbuf_v.push_back(value);
		cbuf_v_lpf.push_back(value);
		cbuf_ts.push_back(ts);
		return 0;
	}
	if(fps==NAG_INF){
		fps=S2NS*BUFSZ*1.0/(cbuf_ts.back()-cbuf_ts[0]);
		winsz=int(fps/3);
		/*cout<<"fps, winsz: "<<fps<<", "<<winsz<<endl;*/
		hamWin=new vector<double>(getHammingWin(winsz));
	}

	if(cbuf_v.capacity()<winsz*2){
		cbuf_v.reserve(winsz*2);
		cbuf_v_lpf.reserve(winsz*2);
	}
	if(cbuf_v.size()<winsz*2){
		cbuf_v.push_back(value);
		//cbuf_v_lpf.push_back(value);
		cbuf_v_lpf.push_back(getDataLpf(cbuf_v, *hamWin));
		i=winsz-2;
		return 0;
	}

	cbuf_v.push_back(value);
	// 	cbuf_v_lpf.push_back(value);
	cbuf_v_lpf.push_back(getDataLpf(cbuf_v, *hamWin));

	if(doLpf)
		cbufPtr=&cbuf_v_lpf;
	else
		cbufPtr=&cbuf_v;

	double v=(*cbufPtr)[winsz-1];
	/*cout<<"v: "<<v<<endl;*/
	i++;

	if(maxMinVar((*cbufPtr), 0, winsz)<varTh &&
		maxMinVar((*cbufPtr), winsz-1, (*cbufPtr).size())<varTh
		){
			return steps;
	}

	double th=a*1.f/(i-k)+b;
	//cout<<"th: "<<th<<endl;

	if(v<(*cbufPtr)[winsz] && maxv-v>=th && maxv-v>baseTh){
		// 		if(i>100 && i<120)
		// 			cout<<"~~~~~~"<<i<<", "<<maxv<<", "<<v<<", "<<maxv-v<<endl;
		steps++;
		k=i;
		maxv=v;
		/*cout<<"~~~~~~~~~~"<<i<<endl;*/
	}
	else{
		//cout<<"---v, maxv: "<<v<<", "<<maxv<<endl;
		if(v>maxv){
			//cout<<"v, maxv: "<<v<<", "<<maxv<<endl;
			maxv=v;
			k=i;
		}
	}

	return steps;
}//dyThresholdOnline

void resetDythCounter(){
	k=0;
	maxv=NAG_INF;
	i=0;
	steps=0;
	start=false;
	fps=NAG_INF;
	//�൱��offline��numtaps
	winsz=NAG_INF;

	cbuf_ts.clear();
	cbuf_v.clear();
	cbuf_v_lpf.clear();
}//resetDythCounter

size_t getDythSteps(){
	return steps;
}//getDythSteps
