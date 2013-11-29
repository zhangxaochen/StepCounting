#ifndef UTILS_H
#define UTILS_H
#define _SCL_SECURE_NO_WARNINGS

#include "scMethods.h"
#include <cstddef>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <pugixml.hpp>

#define _USE_MATH_DEFINES
#include <math.h>
// #include <cmath>


namespace fs=boost::filesystem;

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::map;
using std::ifstream;
using boost::mpl::pair;
using namespace boost::accumulators;
using namespace pugi;

typedef map<string, int>gtmap_type;

const double PI = 2*acos(0.0);

//python 里得到的， 取 -numtaps=10, cutoff=40, nyq=800
const double ham10[]={0.01509692, 0.03662365, 0.09209563, 0.15669359, 0.19949021, 0.19949021, 0.15669359, 0.09209563, 0.03662365, 0.01509692};

// old-style
const char *const kRoot="CaptureSession",
	*const kNodes="Nodes",
	*const kNode="Node",
	*const kData="Data",

	*const kTs="ts",
	*const kTimestamp="timestamp",

	*const kAx="Ax",
	*const kAy="Ay",
	*const kAz="Az",
	*const kGx="Gx",
	*const kGy="Gy",
	*const kGz="Gz",
	*const kMx="Mx",
	*const kMy="My",
	*const kMz="Mz",
	*const kRx="Rx",
	*const kRy="Ry",
	*const kRz="Rz",
	*const kRw="Rw";

/**根据key返回传感器数据 vector
 **@param nNode pugiXml 解析得到的 old-style 的 nodeNode 节点
 **@param key 需要哪个轴数据
 **@return vector<double> 返回对应轴数据序列
 **/
vector<double> getDataList(const xml_node &nNode, const char *key);
vector<double> getDataList(const char *fname, const char *key);

vector<vector<double> > getAxyzBF(const xml_node &nNode);

vector<vector<double> > getAxyzBF(const char *fname);

template<typename T>
vector<double> getDataLpf(const vector<T> &data, const vector<double> &hamWin){
	vector<double> res;
	for(size_t n=0; n<data.size(); n++){
		double y=0;
		for(size_t i=0; i<hamWin.size(); i++){
			if(n<i)
				break;
			y+=hamWin[i]*data[n-i];
		}
		res.push_back(y);
	}
	return res;
}//getDataLpf

template<typename T, size_t N>
vector<double> getDataLpf(const vector<T> &data, const double (&hamWin)[N]){
	vector<double> vecWin(hamWin, hamWin+N);
	return getDataLpf(data, vecWin);
}//getDataLpf

/**对 data 最后一帧进行 lpf
 **
 **/
template<typename T>
double getDataLpf(const T &data, const vector<double> &hamWin){
	size_t n=data.size()-1;
	double y=0;
	for(size_t i=0; i<hamWin.size(); i++){
		if(n<i)
			break;
		y+=hamWin[i]*data[n-i];
	}
	return y;
}//getDataLpf

template<typename T>
double getDataLpf2(const T &data, const vector<double> &hamWin){
	size_t n=data.size()-1;
	if(n<hamWin.size()-1)
		return data[n];
	
	double y=0;
	for(size_t i=0; i<hamWin.size(); i++)
		y+=hamWin[i]*data[n-i];
	return y;
}//getDataLpf

/** get the population variance of src[begin:end]
 ** @param src The input src vector
 ** @param begin Begin index in src
 ** @param end End index in src, src[end] is not included
 ** @return The population variance 
 **/
template<typename T>
double getVarp(const vector<T> &src, size_t begin, size_t end){
	assert(begin>=0 && end>0);
	accumulator_set<T, stats<tag::lazy_variance> > acc_var;
	for(size_t i=begin; i<end; i++){
		acc_var(src[i]);
	}
	double var=variance(acc_var);
	return var;
}//getVarp

/** get the population variance of src
 ** @param src The input src vector
 ** @return The population variance 
 **/
template<typename T>
double getVarp(const vector<T> &src){
	return getVarp(src, 0, src.size());
}//getVarp

template<typename T, size_t N>
double getVarp(const T (&src)[N],  size_t begin, size_t end){
	assert(begin>=0 && end>0);
	accumulator_set<T, stats<tag::lazy_variance> > acc_var;
	for(size_t i=begin; i<end; i++){
		acc_var(src[i]);
	}
	return variance(acc_var);
}//getVarp

template<typename T, size_t N>
double getVarp(const T (&src)[N]){
	return getVarp(src, 0, N);
}//getVarp

/**对每一帧向前数 numtaps 个，求 var
 **@param data 源
 **@param numtaps 向前窗口大小
 **/
template<typename T>
vector<double> getVarPrev(const vector<T> &data, size_t numtaps){
	vector<double> res;
	for(size_t i=0; i<data.size(); i++){
		double v=0;
		if(i<numtaps-1)
			v=getVarp(data, 0, i+1);
		else
			v=getVarp(data, i-(numtaps-1), i+1);
		res.push_back(v);
	}
	return res;
}//getVarPrev

// template<typename T>
// double maxMinVar(const vector<T> &data, size_t begin, size_t end){
// 	assert(begin>=0 && end<=data.size());
// 	T max, min;
// 	max=min=data[begin];
// 	for(size_t i=begin; i<end; i++){
// 		T v=data[i];
// 		if (v>max)
// 			max=v;
// 		else if (v<min)
// 			min=v;
// 	}
// 	double t[]={min, max};
// 	return getVarp(t);
// }//maxMinVar

template<typename T>
double maxMinVar(const T &data, size_t begin, size_t end){
	assert(begin>=0 && end<=data.size());
	double max, min;
	max=min=data[begin];
	for(size_t i=begin; i<end; i++){
		double v=data[i];
		if (v>max)
			max=v;
		else if (v<min)
			min=v;
	}
	double t[]={min, max};
	return getVarp(t);
}//maxMinVar

// template<typename T>
// double maxMinVar(const vector<T> &data){
// 	return maxMinVar(data, 0, data.size());
// }//maxMinVar

template<typename T>
double maxMinVar(const T &data){
	return maxMinVar(data, 0, data.size());
}//maxMinVar

/**获取长度n的海明窗
 **@param n 海明窗长度
 **/
vector<double> getHammingWin(size_t n);

vector<fs::path> listFiles(const fs::path &dir, const string &ext);

gtmap_type getGtMap(const string &gtFname);

void countAndPrint(const char *dataDir, const char* gtFname, int whichMethod=METHOD_DYTH);


#endif //UTILS_H
