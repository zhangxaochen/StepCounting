#include "utils.h"


//�Լ����� var�� sd�ο���
// http://www.softwareandfinance.com/CPP/MeanVarianceStdDevi.html

/**����key���ش��������� vector
 **@param nNode pugiXml �����õ��� old-style �� nodeNode �ڵ�
 **@param key ��Ҫ�ĸ�������
 **@return vector<double> ���ض�Ӧ����������
 **/
vector<double> getDataList(const xml_node &nNode, const char *key){
	vector<double> res;
	for(auto &dnode: nNode.children()){
		double v=dnode.attribute(key).as_double();
		res.push_back(v);
	}
	return res;
}//getDataList

vector<double> getDataList(const char *fname, const char *key){
	xml_document doc;
	xml_parse_result res=doc.load_file(fname);
	assert(res.status==xml_parse_status::status_ok);
	xml_node &nNode=doc.child(kRoot).child(kNodes).child(kNode);
	return getDataList(nNode, key);
}//getDataList

vector<vector<double>> getAxyzBF(const xml_node &nNode){
	vector<vector<double>> res;
	res.resize(2);

	for(auto &dnode: nNode.children()){
		double ts=dnode.attribute(kTs).as_double();
		if(abs(ts-0.0)<epsilon)
			ts=dnode.attribute(kTimestamp).as_double();
		
		double ax=dnode.attribute(kAx).as_double();
		double ay=dnode.attribute(kAy).as_double();
		double az=dnode.attribute(kAz).as_double();
		double v=sqrt(ax*ax+ay*ay+az*az);
		res[0].push_back(ts);
		res[1].push_back(v);
	}
	return res;
}//getAxyzBF

vector<vector<double>> getAxyzBF(const char *fname){
	xml_document doc;
	xml_parse_result res=doc.load_file(fname);
	assert(res.status==xml_parse_status::status_ok);
	xml_node &nNode=doc.child(kRoot).child(kNodes).child(kNode);
	return getAxyzBF(nNode);
}//getAxyzBF

vector<double> getHammingWin(size_t n){
	vector<double> res(n);
	double sum=0;
	for(size_t i=0; i<n; i++){
		res[i]=0.54-0.46*cos(2*M_PI*i/(n-1));
		sum+=res[i];
	}
	for(auto &i:res){
		i/=sum;
	}
	return res;
}//getHammingWin

vector<fs::path> listFiles(const fs::path &dir, const string &ext){
	vector<fs::path> res;
	if(fs::exists(dir) && fs::is_directory(dir)){
		fs::directory_iterator it(dir);
		fs::directory_iterator endit;
		while(it!=endit){
			//cout<<it->path().extension()<<endl;
			
			if(fs::is_regular_file(*it) && 
				it->path().extension().string().find(ext)!=string::npos)
				res.push_back(it->path().filename());
			++it;
		}
	}
	return res;
}//listFiles

map<string, int> getGtMap(const string &gtFname){
	map<string, int> res;

	ifstream ifs(gtFname);
	//vector<string> lines;
	string key;
	int steps;
	while(ifs>>key>>steps){
		//res.insert(pair<string, int>(key, steps));
		res[key]=steps;
	}
	return res;
}//getGtMap

/**��ӡ�Ʋ������������
 **@param dataDir ���������ļ���·��
 **@param gtFname ground truth �ļ�·��
 **/
void countAndPrint(const char *dataDir, const char* gtFname, int whichMethod){
	gtmap_type gtmap=getGtMap(gtFname);
	vector<fs::path> fnames=listFiles(dataDir, ".xml");
	double errSum=0;

	vector<double> errors;
	vector<double> gts;
	for(auto fname:fnames){
		string path=string(dataDir)+"/"+fname.string();
		vector<vector<double>> alist=getAxyzBF(path.c_str());
		vector<double> tslist=alist[0],
			axyzlist=alist[1];
		for(size_t i=0; i<tslist.size(); i++){
			pushData(axyzlist[i], tslist[i], whichMethod);
		}
		int steps=getSteps(whichMethod);
		resetCounter(whichMethod);

// 		string strFname=fname.string();
// 		int pos=strFname.find('_');
// 		//ʵ��id��
// 		string eid=strFname.substr(0, pos-1);
		vector<string> tmp;
		string eid=boost::split(tmp, fname.string(), boost::is_any_of("_"))[0];
		//cout<<"eid: "<<eid<<endl;
		
		int gtv=gtmap[eid];
		double err=abs(gtv-steps)*100.0/gtv;
		errSum+=err;
		cout<<fname.string()<<"\t"<<gtv<<"\t"<<steps<<"\t"<<err<<"%"<<endl;

		errors.push_back(err);
		gts.push_back(gtv);
	}
	//double meanErr=errSum/fnames.size();
	double meanErr=0;
	double gtsum=0;
	for(size_t i=0; i<gts.size(); i++){
		gtsum+=gts[i];
	}
	for(size_t i=0; i<gts.size(); i++){
		meanErr+=errors[i]*gts[i]/gtsum;
	}
	cout<<"mean error: "<<meanErr<<endl;
}//countAndPrint
