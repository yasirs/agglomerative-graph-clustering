#ifndef GRAPHDATA_HPP
#define GRAPHDATA_HPP
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>

void my_Tokenize(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiter);

std::string getstvec(std::vector<std::string> a, int i) {
	return a[i];
}

void putstvec(std::vector<std::string> a, int i, const char* st) {
	a[i] = st;
}

class graphData{
	public:
		int numV;
		char gtype;
		float aveP;
		float Etot;
		typedef std::map<int, float> destList;
		std::map<int, std::string> int2Name;
		std::map<std::string, int> name2Int;
		std::map<int, destList*> edgeList;
		float get_uv(int u, int v);
		void set_uv(int u, int v, float w);
		bool has_uv(int u, int v);
		int Add_uv(int u, int v, float w); 
		bool readWeighted(const char* filename);
		bool readBinary(const char* filename);
		void writeBoth(const char* filename);
		void writeSingle(const char* filename);
		int degree(int i);
		//std::set<int> neighbors(int i);
		~graphData() {
			std::map<int, destList*>::iterator it;
			for (it=edgeList.begin(); it!=edgeList.end(); ++it) {
				(*it).second->clear();
				delete (*it).second;
			}
			edgeList.clear();
		}
};


void graphData::writeSingle(const char* fn) {
	int u,v;
	float weight;
	std::ofstream file;
	file.open(fn,std::ios::out);
	std::map<int, destList*>::iterator outIt;
	destList::iterator inIt;
	for (outIt = edgeList.begin(); outIt != edgeList.end(); ++outIt) {
		u = outIt->first;
		for (inIt = edgeList[u]->begin(); inIt != edgeList[u]->end(); ++inIt) {
			v = inIt->first;
			weight = inIt->second;
			if (u<=v) {
				if (gtype=='w') {
					file << int2Name[u] <<'\t' << int2Name[v] << '\t' << weight << '\n';
				}
				else if (gtype=='b') {
					file << int2Name[u] <<'\t' << int2Name[v] << '\n';
				}
				else {
					std::cerr << gtype << " type of graph not recognized in writing graphs\n";
					throw 1;
				}
			}
		}
	}
		
};



void graphData::writeBoth(const char* fn) {
	int u,v;
	float weight;
	std::ofstream file;
	file.open(fn,std::ios::out);
	std::map<int, destList*>::iterator outIt;
	destList::iterator inIt;
	for (outIt = edgeList.begin(); outIt != edgeList.end(); ++outIt) {
		u = outIt->first;
		for (inIt = edgeList[u]->begin(); inIt != edgeList[u]->end(); ++inIt) {
			v = inIt->first;
			weight = inIt->second;
			if (gtype=='w') {
				file << int2Name[u] <<'\t' << int2Name[v] << '\t' << weight << '\n';
			}
			else if (gtype=='b') {
				file << int2Name[u] <<'\t' << int2Name[v] << '\n';
			}
			else {
				std::cerr << gtype << " type of graph not recognized in writing graphs\n";
				throw 1;
			}
		}
	}
		
};


bool graphData::has_uv(int u, int v) {
	if (edgeList.find(u)!=edgeList.end()) {
		if (edgeList[u]->find(v) != edgeList[u]->end()) {
			return 1;
		}
	}
	return 0;
};


float graphData::get_uv(int u, int v) {
	if (edgeList.find(u)!=edgeList.end()) {
		if (edgeList[u]->find(v) != edgeList[u]->end()) {
			return (*edgeList[u])[v];
		}
	}
	return 0;
};



bool graphData::readWeighted(const char* fn) {
	gtype = 'w';
	Etot = 0.0f;
	std::string strline;
	std::ifstream file;
	std::vector<std::string> tok;
	destList* pdl;
	int u,v;
	float sum = 0;
	float weight;
	std::string St1, St2;
	std::istringstream temp;
	file.open(fn,std::ios::in);
	if (! file.is_open()) return 0;
	while (!file.eof()) {
		getline(file,strline);
		tok.clear();
		my_Tokenize(strline,tok," \t");
		if (tok.size()>1) {
			if (name2Int.find(tok[0])==name2Int.end()) {
				u = name2Int.size();
				name2Int[tok[0]] = u;
				int2Name[u] = tok[0];
			} else {
				u = name2Int[tok[0]];
			}
			if (name2Int.find(tok[1])==name2Int.end()) {
				v = name2Int.size();
				name2Int[tok[1]] = v;
				int2Name[v] = tok[1];
			} else {
				v = name2Int[tok[1]];
			}
			if (tok.size()>2) {
				std::istringstream(tok[2]) >> weight;
			} else weight = 1;
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
	numV = name2Int.size();
	aveP = sum/(numV * numV); // shouldn't matter because aveP is to be used for binary
	return 1;
}

void graphData::set_uv(int u, int v, float w) {
	if (edgeList.find(u)==edgeList.end()) {
		edgeList[u] = new destList;
	}
	(*edgeList[u])[v] = w;
};

int graphData::Add_uv(int u, int v, float w) {
	if (edgeList.find(u)==edgeList.end()) {
		edgeList[u] = new destList;
		(*edgeList[u])[v] = w;
		return 0;
	} else {
		if ((*edgeList[u]).find(v) != (*edgeList[u]).end()) {
			(*edgeList[u])[v] += w;
			return 1;
		} else {
			(*edgeList[u])[v] = w;
			return 0;
		}
	}
};


	


bool graphData::readBinary(const char* fn) {
	gtype = 'b';
	Etot = 0.0f;
	std::string strline;
	std::ifstream file;
	std::vector<std::string> tok;
	destList* pdl;
	std::istringstream temp;
	int u,v;
	float sum = 0;
	file.open(fn,std::ios::in);
	if (! file.is_open()) return 0;
	while (!file.eof()) {
		getline(file,strline);
		tok.clear();
		my_Tokenize(strline,tok," \t");
		if (tok.size()>1) {
			if (name2Int.find(tok[0])==name2Int.end()) {
				u = name2Int.size();
				name2Int[tok[0]] = u;
				int2Name[u] = tok[0];
			} else {
				u = name2Int[tok[0]];
			}
			if (name2Int.find(tok[1])==name2Int.end()) {
				v = name2Int.size();
				name2Int[tok[1]] = v;
				int2Name[v] = tok[1];
			} else {
				v = name2Int[tok[1]];
			}
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
	numV = name2Int.size();
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

