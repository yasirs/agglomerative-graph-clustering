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
		float Etot;
		typedef std::map<int, float> destList;
		std::map<int, destList*> edgeList;
		bool readWeighted(const char* filename);
		bool readBinary(const char* filename);
		int degree(int i);
		std::set<int> neighbors(int i);
		~graphData() {
			std::map<int, destList*>::iterator it;
			for (it=edgeList.begin(); it!=edgeList.end(); ++it) {
				(*it).second->clear();
				delete (*it).second;
			}
			edgeList.clear();
		}
};


bool graphData::readWeighted(const char* fn) {
	gtype = 'w';
	numV = 0;
	Etot = 0.0f;
	std::string strline;
	std::ifstream file;
	std::vector<std::string> tok;
	destList* pdl;
	int u,v;
	float sum = 0;
	float weight;
	file.open(fn,std::ios::in);
	if (not file.is_open()) return 0;
	while (!file.eof()) {
		getline(file,strline);
		tok.clear();
		my_Tokenize(strline,tok," \t");
		if (tok.size()>1) {
			std::istringstream(tok[0]) >> u;
			std::istringstream(tok[1]) >> v;
			if (tok.size()>2) {
				std::istringstream(tok[2]) >> weight;
			} else weight = 0;
			if (numV<(u+1)) numV = u+1;
			if (numV<(v+1)) numV = v+1;
			if (edgeList.find(u)==edgeList.end()) {
				pdl = new destList;
				edgeList[u] = pdl;
			}
			(*edgeList[u])[v] = weight;
			if (edgeList.find(v)==edgeList.end()) {
				pdl = new destList;
				edgeList[v] = pdl;
			}
			(*edgeList[v])[u] = weight;
			sum += 2; // NOTE : this shouldn't really matter as the sum isn't used for
					// weighted graphs
			Etot += weight;
		}
	}
	aveP = sum/(numV * numV);
	return 1;
}




bool graphData::readBinary(const char* fn) {
	gtype = 'b';
	numV = 0;
	Etot = 0.0f;
	std::string strline;
	std::ifstream file;
	std::vector<std::string> tok;
	destList* pdl;
	int u,v;
	float sum = 0;
	file.open(fn,std::ios::in);
	if (not file.is_open()) return 0;
	while (!file.eof()) {
		getline(file,strline);
		tok.clear();
		my_Tokenize(strline,tok," \t");
		if (tok.size()>1) {
			std::istringstream(tok[0]) >> u;
			std::istringstream(tok[1]) >> v;
			if (numV<(u+1)) numV = u+1;
			if (numV<(v+1)) numV = v+1;
			if (edgeList.find(u)==edgeList.end()) {
				pdl = new destList;
				edgeList[u] = pdl;
			}
			(*edgeList[u])[v] = 1.0f;
			if (edgeList.find(v)==edgeList.end()) {
				pdl = new destList;
				edgeList[v] = pdl;
			}
			(*edgeList[v])[u] = 1.0f;
			sum += 2;
			Etot += 1.0f;
		}
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

