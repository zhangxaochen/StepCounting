
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
	/*doLpf==false, ��������ѵ�ͨ�˲���������*/
	for(size_t j=00; j<tslist.size(); j++)
		s2=dyThresholdOnline(axyzlpf[j], tslist[j]);
	cout<<"online: "<<s2<<endl;	//1*ʵ�ʲ�����ͳ��ֵ��

	/*doLpf==true, ����ԭʼ����*/
	for(size_t j=00; j<tslist.size(); j++)
		s2=dyThresholdOnline(list[1][j], tslist[j], true);
	cout<<"online: "<<s2<<endl;	//2*ʵ�ʲ�����ͳ��ֵ��, ��Ϊû�� reset �Ļ��Ʋ����ϴ��ۻ�
	resetDythCounter();	//���üƲ���

	int olds2=0;
	for(size_t j=00; j<tslist.size(); j++){
		s2=dyThresholdOnline(list[1][j], tslist[j], true);
		if(s2!=olds2){
			//DO SOMETHING...
			//e.g., �����¼��� �����ź�(Qt)
			olds2=s2;
		}
	}
	cout<<"online: "<<s2<<endl;	//1*ʵ�ʲ�����ͳ��ֵ������Ϊ reset ����
	//Ҳ������ getSteps �õ�����
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

