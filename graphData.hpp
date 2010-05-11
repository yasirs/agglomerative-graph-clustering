#ifndef GRAPHDATA_HPP
#define GRAPHDATA_HPP
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

void my_Tokenize(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiter);

class graphData{
	public:
		int numV;
		char gtype;
		float aveP;
		typedef std::map<int, float> destList;
		std::map<int, destList> edgeList;
		bool readWeighted(const char* filename);
		bool readBinary(const char* filename);
		int degree(int i);
		std::set<int> neighbors(int i);
};

bool graphData::readBinary(const char* fn) {
	gtype = 'b';
	numV = 0;
	std::string strline;
	std::ifstream file;
	std::vector<std::string> tok;
	int u,v;
	float sum = 0;
	file.open(fn,std::ios::in);
	if (not file.is_open()) return 0;
	while (!file.eof()) {
		getline(file,strline);
		tok.clear();
		my_Tokenize(strline,tok," \t");
		std::istringstream(tok[0]) >> u;
		std::istringstream(tok[1]) >> v;
		if (numV<(u+1)) numV = u+1;
		if (numV<(v+1)) numV = v+1;
		if (edgeList.find(u)==edgeList.end()) {
			edgeList[u] = destList();
		}
		edgeList[u][v] = 1.0f;
		if (edgeList.find(v)==edgeList.end()) {
			edgeList[v] = destList();
		}
		edgeList[v][u] = 1.0f;
		sum += 2;
	}
	aveP = sum/(numV * numV);
	return 1;
}

void my_Tokenize(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiters = "")  // via http://www.geocities.com/eric6930/cplus.html unknown author presumed Eric.
{
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
}





#endif

