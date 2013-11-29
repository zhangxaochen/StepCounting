
#include "utils.h"
#include "scMethods.h"

#include <iostream>
#include <algorithm>
using namespace std;

#include <vector>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>



int main(){
	const char *fname="./data/ZC30fast_a9_0.xml";
	//fname="./data/ZC30fast_a9_7.xml";
	fname="./data/ZCrun120_a5_0.xml";
	vector<vector<double>> list=getAxyzBF(fname);
	vector<double> tslist=list[0];
	vector<double> axyzlpf=getDataLpf(list[1], ham10);

	//--------------------dyTh
	cout<<"----------dyThreshold:"<<endl;
	int steps=dyThresholdOffline(axyzlpf);
	cout<<"offline: "<<steps<<endl;
	
	int s2=0;
	/*doLpf==false, 传入的是已低通滤波过的数据*/
	for(size_t j=00; j<tslist.size(); j++)
		s2=dyThresholdOnline(axyzlpf[j], tslist[j]);
	cout<<"online: "<<s2<<endl;	//1*实际步数（统计值）

	/*doLpf==true, 传入原始数据*/
	for(size_t j=00; j<tslist.size(); j++)
		s2=dyThresholdOnline(list[1][j], tslist[j], true);
	cout<<"online: "<<s2<<endl;	//2*实际步数（统计值）, 因为没有 reset 的话计步接上次累积
	resetDythCounter();	//重置计步器

	int olds2=0;
	for(size_t j=00; j<tslist.size(); j++){
		s2=dyThresholdOnline(list[1][j], tslist[j], true);
		if(s2!=olds2){
			//DO SOMETHING...
			//e.g., 触发事件， 或发射信号(Qt)
			olds2=s2;
		}
	}
	cout<<"online: "<<s2<<endl;	//1*实际步数（统计值），因为 reset 过了
	//也可以用 getSteps 得到步数
	cout<<"online, getSteps: "<<getDythSteps()<<endl;
	resetDythCounter();
	
	const char *dataDir="./data";
	const char *gtFname="./data/groundtruth_SC.txt";
	countAndPrint(dataDir, gtFname);

	cout<<"----------dyPeak:"<<endl;
	for(size_t j=00; j<tslist.size(); j++)
		s2=dyPeakOnline(list[1][j], tslist[j], true);
	cout<<"online: "<<s2<<endl;
	resetDypeakCounter();

	countAndPrint(dataDir, gtFname, METHOD_DYPEAK);

	cout<<"----------dyZcross:"<<endl;
	for(size_t j=00; j<tslist.size(); j++)
		s2=dyZcrossOnline(list[1][j], tslist[j], true);
	cout<<"online: "<<s2<<endl;
	resetDyzcCounter();

	countAndPrint(dataDir, gtFname, METHOD_DYZCROSS);
	
	return 0;
}//main

