#ifndef LABELDATA_HPP
#define LABELDATA_HPP
#include "graphData.hpp"
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <tr1/unordered_map>

void my_Tokenize(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiter);

int getcm(std::set<int>* s1, std::set<int>* s2) {
	std::set<int>::iterator x1, x2, l1, l2;
	x1 = s1->begin(); x2 = s2->begin();
	l1 = s1->end(); l2 = s2->end();
	int c=0;
	while ((x1 != l1)&&(x2 != l2)) {
		if (*x1>*x2) {
			x2++;
		} else if (*x2 > *x1) {
			x1++;
		}
		else {
			c++;
			x1++;
			x2++;
		}
	}
	return c;
}


class labelData{
	public:
		typedef std::tr1::unordered_map<int, float*> destList;
		std::map<int, std::string> int2Name;
		std::map<std::string, int> name2Int;
		std::map<int, destList*> edgeList;
		std::vector<std::string> headers;
		int levels;
		labelData(graphData* Gorginal, const char* filename, int reg);
		void putSoFar(graphData* GsoFar, int l);
		void write(const char* filename);
		void populateLocal(const char* filename,int lput);
		~labelData() {
			std::map<int, destList*>::iterator it;
			destList::iterator Init;
			for (it=edgeList.begin(); it!=edgeList.end(); ++it) {
				for (Init = it->second->begin(); Init != it->second->end(); ++Init) {
					delete[] Init->second;
				}
				(*it).second->clear();
				delete (*it).second;
			}
			edgeList.clear();
		}
};

void labelData::populateLocal(const char* filename, int lput) {
	std::ifstream file;
	std::string strline;
	std::vector<std::string> tok;
	int u,v;
	file.open(filename,std::ios::in);
	typedef std::set<int> dlist;
	std::map<int, dlist> elist;
	if (! file.is_open()) throw(0);
	while (!file.eof()) {
		getline(file,strline);
		tok.clear();
		my_Tokenize(strline,tok," \t");
		if (tok.size()>1) {
			if (name2Int.find(tok[0])==name2Int.end()) {
				continue;
			} else {
				u = name2Int[tok[0]];
			}
			if (name2Int.find(tok[1])==name2Int.end()) {
				continue;
			} else {
				v = name2Int[tok[1]];
			}
			elist[u].insert(v);
			elist[v].insert(u);
		}
	}
	float du, dv, cm;
	std::set<int>* destu;
	std::set<int>* destv;
	std::map<int, destList*>::iterator outIt;
	destList::iterator inIt;
	for (outIt = edgeList.begin(); outIt != edgeList.end(); ++outIt) {
		u = outIt->first;
		destu = &elist[u];
		du = destu->size();
		for (inIt = outIt->second->begin(); inIt != outIt->second->end(); ++inIt) {
			v = inIt->first;
			destv = &elist[v];
			dv = destv->size();
			cm = getcm(destu, destv);
			inIt->second[lput] = du*dv;
			inIt->second[lput+1] = cm;
			inIt->second[lput+2] = ( (float) cm)/(du+dv-cm);
		}
	}
	headers[lput] = std::string("dprod");
	headers[lput+1] = std::string("cneighb");
	headers[lput+2] = std::string("jaccard");
}

	

void labelData::putSoFar(graphData* GsoFar, int levthis) {
	int u,v;
	std::map<int, destList*>::iterator outIt;
	destList::iterator inIt;
	for (outIt = edgeList.begin(); outIt != edgeList.end(); ++outIt) {
		u = outIt->first;
		for (inIt = outIt->second->begin(); inIt != outIt->second->end(); ++inIt) {
			v = inIt->first;
			inIt->second[levthis] = GsoFar->get_uv(u,v);
		}
	}
	std::stringstream temp;
	temp << levthis;
	std::string lstr;
	temp >> lstr;
	headers[levthis] = std::string("soFar") + lstr;
};

void labelData::write(const char* fn) {
	int u,v;
	std::ofstream file;
	file.open(fn,std::ios::out);
	std::map<int, destList*>::iterator outIt;
	destList::iterator inIt;
	std::string uname, vname;
	file << "u\tv";
	for (int i=0;i<= levels; i++) {
		file << '\t' << headers[i];
	}
	file << '\n';
	for (outIt = edgeList.begin(); outIt != edgeList.end(); ++outIt) {
		u = outIt->first;
		uname = int2Name[u];
		for (inIt = outIt->second->begin(); inIt != outIt->second->end(); ++inIt) {
			v = inIt->first;
			vname = int2Name[v];
			file << uname << '\t' << vname;
			for (int i=0; i <= levels; i++) {
				file << '\t' << inIt->second[i];
			}
			file << '\n';
		}
	}
	file.close();
};



labelData::labelData(graphData* Goriginal, const char* filename, int reg) {
	this->levels = reg;
	headers.resize(reg+1);
	headers[0] = std::string("label");
	std::string strline;
	std::ifstream file;
	std::vector<std::string> tok;
	destList* pdl;
	int u,v;
	float weight;
	int linesread=0;
	std::string St1, St2;
	file.open(filename,std::ios::in);
	if (! file.is_open()) throw( 0);
	this->int2Name = Goriginal->int2Name;
	this->name2Int = Goriginal->name2Int;
	while (! file.eof()) {
		getline(file,strline);
		tok.clear();
		my_Tokenize(strline,tok," \t");
		if (tok.size()>1) {
			linesread++;
			if (Goriginal->name2Int.find(tok[0])==Goriginal->name2Int.end()) {
				// do nothing
				continue;
				throw( 0);
			} else {
				u = Goriginal->name2Int[tok[0]];
			}
			if (Goriginal->name2Int.find(tok[1])==Goriginal->name2Int.end()) {
				// do nothing
				continue;
				throw( 0);
			} else {
				v = Goriginal->name2Int[tok[1]];
			}
			std::stringstream temp;
			temp << tok[2];
			temp >> weight;
			if (edgeList.find(u)==edgeList.end()) {
				pdl = new destList;
				edgeList[u] = pdl;
			}
			float* fp = new float[levels+1];
			fp[0] = weight;
			(*edgeList[u])[v] = fp;
		}
	}
};








#endif

